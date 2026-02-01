from Bio import Entrez, SeqIO
import time
from collections import defaultdict
from datetime import timedelta
import os
import io
import re

OUTPUT_DIR = "FINAL_OUTPUT"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ================= CONFIG =================
Entrez.email = "a71364@ualg.pt"

INPUT_FILE = "Galaxy14_PHMMER.txt"
OUTPUT_DIR = "FINAL_OUTPUT"

SLEEP = 0.34
PROGRESS_STEP = 10
# =========================================


FAMILIES = {
    "CALCR": ["calcitonin", "calcr"],
    "CRHR": ["corticotropin", "crhr"],
    "PTHR": ["parathyroid", "pthr"],
    "GCGR": ["glucagon", "gcgr"],
    "VIP_SCTR": [
        "vasoactive intestinal", "vipr",
        "secretin", "sctr"
    ],
}



def is_np(acc):
    return acc.startswith("NP_")


def get_gene_id(acc):
    try:
        h = Entrez.elink(dbfrom="protein", db="gene", id=acc)
        r = Entrez.read(h)
        h.close()
        links = r[0].get("LinkSetDb", [])
        if links:
            return links[0]["Link"][0]["Id"]
    except:
        pass
    return None


def get_protein_sequence(acc):
    try:
        h = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        seq = SeqIO.read(io.StringIO(h.read()), "fasta")
        h.close()
        return str(seq.seq)
    except:
        return "NA"


def get_species_abbrev(description_line):
    try:
        match = re.findall(r"\[([^\]]+)\]", description_line)
        if not match:
            return "NA"

        species = match[-1].strip()
        parts = [p for p in species.split() if p]
        if len(parts) < 2:
            return "NA"

        first = parts[0]
        last = parts[-1]
        abbrev = f"{first[0].upper()}{last[:2].lower()}"
        return abbrev
    except:
        return "NA"


def format_sequence_with_header(acc, sequence, description_line):
    abbrev = get_species_abbrev(description_line)
    if sequence == "NA":
        return f">{abbrev} {acc} NA"
    return f">{abbrev} {acc} {sequence}"


def get_gene_summary(gene_id):
    try:
        h = Entrez.esummary(db="gene", id=gene_id, retmode="xml")
        summary = Entrez.read(h)
        h.close()

        docset = summary.get("DocumentSummarySet", {})
        docs = docset.get("DocumentSummary", [])
        if not docs:
            return {}

        doc = docs[0]
        genomic = (doc.get("GenomicInfo") or [{}])[0]

        exon_count = genomic.get("ExonCount")
        chromosome = doc.get("Chromosome") or genomic.get("ChrLoc")
        chr_accver = genomic.get("ChrAccVer")
        chr_start = genomic.get("ChrStart")
        chr_stop = genomic.get("ChrStop")

        return {
            "exon_count": exon_count,
            "chromosome": chromosome,
            "chr_accver": chr_accver,
            "chr_start": chr_start,
            "chr_stop": chr_stop,
        }
    except:
        return {}


# ================= PIPELINE =================
os.makedirs(OUTPUT_DIR, exist_ok=True)

records = []
unique_accessions = []
seen = set()

# ---------- READ INPUT ----------
with open(INPUT_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        parts = line.split()
        acc = parts[0]
        tlen = int(parts[2])

        records.append({
            "line": line.strip(),
            "accession": acc,
            "tlen": tlen
        })

        if acc not in seen:
            unique_accessions.append(acc)
            seen.add(acc)

total = len(unique_accessions)

gene_cache = {}
seq_cache = {}
gene_info_cache = {}

start = time.time()
next_progress = PROGRESS_STEP

print(f"Total unique accessions: {total}\n")

# ---------- FETCH GENE IDS + SEQUENCES ----------
for i, acc in enumerate(unique_accessions, start=1):
    gene_cache[acc] = get_gene_id(acc)
    seq_cache[acc] = get_protein_sequence(acc)

    percent = (i / total) * 100
    elapsed = time.time() - start

    if percent >= next_progress or i == total:
        eta = (elapsed / i) * (total - i)
        print(
            f"[{int(percent)}%] {i}/{total} | "
            f"elapsed: {timedelta(seconds=int(elapsed))} | "
            f"ETA: {timedelta(seconds=int(eta))}"
        )
        next_progress += PROGRESS_STEP

    time.sleep(SLEEP)

# ---------- FETCH GENE SUMMARY (EXON + CHR POSITION) ----------
unique_gene_ids = sorted({gid for gid in gene_cache.values() if gid})
gene_total = len(unique_gene_ids)
start = time.time()
next_progress = PROGRESS_STEP

print(f"\nTotal unique GeneIDs: {gene_total}\n")

for i, gid in enumerate(unique_gene_ids, start=1):
    gene_info_cache[gid] = get_gene_summary(gid)

    percent = (i / gene_total) * 100 if gene_total else 100
    elapsed = time.time() - start

    if percent >= next_progress or i == gene_total:
        eta = (elapsed / i) * (gene_total - i) if i else 0
        print(
            f"[Gene {int(percent)}%] {i}/{gene_total} | "
            f"elapsed: {timedelta(seconds=int(elapsed))} | "
            f"ETA: {timedelta(seconds=int(eta))}"
        )
        next_progress += PROGRESS_STEP

    time.sleep(SLEEP)

# ---------- ATTACH DATA ----------
for r in records:
    acc = r["accession"]
    r["gene_id"] = gene_cache.get(acc)
    r["sequence"] = seq_cache.get(acc)
    r["sequence_formatted"] = format_sequence_with_header(
        acc,
        r["sequence"],
        r["line"]
    )
    info = gene_info_cache.get(r["gene_id"], {})
    r["exon_count"] = info.get("exon_count") or "NA"
    r["chromosome"] = info.get("chromosome") or "NA"
    r["chr_accver"] = info.get("chr_accver") or "NA"
    r["chr_start"] = info.get("chr_start") or "NA"
    r["chr_stop"] = info.get("chr_stop") or "NA"

# ---------- GROUP BY GENEID ----------
by_gene = defaultdict(list)
for r in records:
    if r["gene_id"]:
        by_gene[r["gene_id"]].append(r)

# ---------- FILTER (largest, NP > XP) ----------
filtered = []
for gid, group in by_gene.items():
    group.sort(
        key=lambda x: (x["tlen"], is_np(x["accession"])),
        reverse=True
    )
    filtered.append(group[0])

# ---------- SPLIT OUTPUT ----------
family_files = {}
for fam in FAMILIES:
    family_files[fam] = open(f"{OUTPUT_DIR}/{fam}.tsv", "w")

family_files["UNCLASSIFIED"] = open(
    f"{OUTPUT_DIR}/UNCLASSIFIED.tsv", "w"
)

header = (
    "accession\ttlen\tGeneID\texon_count\tchromosome\tchr_accver\t"
    "chr_start\tchr_stop\tprotein_sequence\tdescription\n"
)
for fh in family_files.values():
    fh.write(header)

for r in filtered:
    line_lower = r["line"].lower()
    assigned = False

    for fam, keys in FAMILIES.items():
        if any(k in line_lower for k in keys):
            family_files[fam].write(
                f"{r['accession']}\t{r['tlen']}\t{r['gene_id']}\t{r['exon_count']}\t"
                f"{r['chromosome']}\t{r['chr_accver']}\t{r['chr_start']}\t"
                f"{r['chr_stop']}\t{r['sequence_formatted']}\t{r['line']}\n"
            )
            assigned = True
            break

    if not assigned:
        family_files["UNCLASSIFIED"].write(
            f"{r['accession']}\t{r['tlen']}\t{r['gene_id']}\t{r['exon_count']}\t"
            f"{r['chromosome']}\t{r['chr_accver']}\t{r['chr_start']}\t"
            f"{r['chr_stop']}\t{r['sequence_formatted']}\t{r['line']}\n"
        )

for fh in family_files.values():
    fh.close()

print("\nPIPELINE COMPLETE ✅")
print(f"Results written to: {OUTPUT_DIR}/")