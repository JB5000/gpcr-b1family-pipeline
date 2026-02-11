#!/usr/bin/env python3
import argparse
import io
import json
import os
import time
from datetime import datetime, timedelta

from Bio import Entrez, SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract NCBI metadata for PHMMER accessions and store JSON cache"
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Directory with PHMMER .txt/.tsv files (default: phmmer_inputs)",
    )
    parser.add_argument(
        "--cache-out",
        default=None,
        help="Output JSON cache path (default: outputs/extraction_cache_<timestamp>.json)",
    )
    parser.add_argument(
        "--email",
        default="a71364@ualg.pt",
        help="Entrez email",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.34,
        help="Seconds between Entrez requests",
    )
    parser.add_argument(
        "--progress-step",
        type=int,
        default=10,
        help="Progress update step in percent",
    )
    return parser.parse_args()


def sanitize_name(value):
    safe = "".join(ch if ch.isalnum() else "_" for ch in str(value)).strip("_")
    return safe or "output"


def find_input_files(input_dir):
    if input_dir is None:
        input_dir = os.path.join(os.path.dirname(__file__), "phmmer_inputs")
    if not os.path.isdir(input_dir):
        raise SystemExit(f"Input dir not found: {input_dir}")

    files = [
        os.path.join(input_dir, f)
        for f in sorted(os.listdir(input_dir))
        if f.lower().endswith((".txt", ".tsv"))
    ]
    if not files:
        raise SystemExit("No PHMMER input files found.")
    return files


def parse_phmmer_file(path):
    rows = []
    with open(path) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            acc = parts[0]
            try:
                tlen = int(parts[2])
            except ValueError:
                tlen = None
            rows.append(
                {
                    "line_number": line_number,
                    "line": line.strip(),
                    "accession": acc,
                    "tlen": tlen,
                }
            )
    return rows


def get_gene_id(accession):
    try:
        h = Entrez.elink(dbfrom="protein", db="gene", id=accession)
        r = Entrez.read(h)
        h.close()
        links = r[0].get("LinkSetDb", [])
        if links:
            return links[0]["Link"][0]["Id"]
    except Exception:
        pass
    return None


def get_protein_sequence(accession):
    try:
        h = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        seq = SeqIO.read(io.StringIO(h.read()), "fasta")
        h.close()
        return str(seq.seq)
    except Exception:
        return "NA"


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

        return {
            "exon_count": genomic.get("ExonCount") or "NA",
            "chromosome": doc.get("Chromosome") or genomic.get("ChrLoc") or "NA",
            "chr_accver": genomic.get("ChrAccVer") or "NA",
            "chr_start": genomic.get("ChrStart") or "NA",
            "chr_stop": genomic.get("ChrStop") or "NA",
        }
    except Exception:
        return {}


def print_progress(prefix, i, total, start, step_state, step):
    if total <= 0:
        return step_state
    percent = int((i / total) * 100)
    if percent >= step_state or i == total:
        elapsed = int(time.time() - start)
        eta = int((elapsed / i) * (total - i)) if i else 0
        print(
            f"[{prefix} {percent}%] {i}/{total} | "
            f"elapsed: {timedelta(seconds=elapsed)} | "
            f"ETA: {timedelta(seconds=eta)}"
        )
        return step_state + step
    return step_state


def resolve_cache_out(path):
    if path:
        return path
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    return os.path.join(os.getcwd(), "outputs", f"extraction_cache_{run_id}.json")


def main():
    args = parse_args()
    Entrez.email = args.email

    files = find_input_files(args.input_dir)
    file_records = {}
    unique_accessions = []
    seen = set()

    for path in files:
        key = sanitize_name(os.path.splitext(os.path.basename(path))[0])
        rows = parse_phmmer_file(path)
        file_records[key] = {
            "input_file": path,
            "safe_base": key,
            "records": rows,
        }
        for row in rows:
            acc = row["accession"]
            if acc not in seen:
                seen.add(acc)
                unique_accessions.append(acc)

    print(f"Parsed input files: {len(files)}")
    print(f"Unique accessions: {len(unique_accessions)}")

    accession_data = {}
    start = time.time()
    next_step = args.progress_step

    for i, acc in enumerate(unique_accessions, start=1):
        gene_id = get_gene_id(acc)
        sequence = get_protein_sequence(acc)
        accession_data[acc] = {
            "gene_id": gene_id or "NA",
            "sequence": sequence or "NA",
        }
        next_step = print_progress("ACC", i, len(unique_accessions), start, next_step, args.progress_step)
        time.sleep(args.sleep)

    unique_gene_ids = sorted({v["gene_id"] for v in accession_data.values() if v.get("gene_id") not in (None, "NA")})
    print(f"Unique GeneIDs: {len(unique_gene_ids)}")

    gene_summary = {}
    start = time.time()
    next_step = args.progress_step

    for i, gid in enumerate(unique_gene_ids, start=1):
        gene_summary[gid] = get_gene_summary(gid)
        next_step = print_progress("GENE", i, len(unique_gene_ids), start, next_step, args.progress_step)
        time.sleep(args.sleep)

    out = {
        "generated_at": datetime.now().isoformat(),
        "entrez_email": args.email,
        "input_dir": args.input_dir,
        "file_records": file_records,
        "accession_data": accession_data,
        "gene_summary": gene_summary,
    }

    out_path = resolve_cache_out(args.cache_out)
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f)

    print(f"Cache saved: {out_path}")


if __name__ == "__main__":
    main()
