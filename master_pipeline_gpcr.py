from Bio import Entrez, SeqIO
import time
from collections import defaultdict
from datetime import timedelta
import os
import io

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

# ---------- ATTACH DATA ----------
for r in records:
    acc = r["accession"]
    r["gene_id"] = gene_cache.get(acc)
    r["sequence"] = seq_cache.get(acc)

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

header = "accession\ttlen\tGeneID\tprotein_sequence\tdescription\n"
for fh in family_files.values():
    fh.write(header)

for r in filtered:
    line_lower = r["line"].lower()
    assigned = False

    for fam, keys in FAMILIES.items():
        if any(k in line_lower for k in keys):
            family_files[fam].write(
                f"{r['accession']}\t{r['tlen']}\t{r['gene_id']}\t{r['sequence']}\t{r['line']}\n"
            )
            assigned = True
            break

    if not assigned:
        family_files["UNCLASSIFIED"].write(
            f"{r['accession']}\t{r['tlen']}\t{r['gene_id']}\t{r['sequence']}\t{r['line']}\n"
        )

for fh in family_files.values():
    fh.close()

print("\nPIPELINE COMPLETE ✅")
print(f"Results written to: {OUTPUT_DIR}/")
