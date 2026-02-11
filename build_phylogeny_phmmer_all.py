#!/usr/bin/env python3
import argparse
import json
import os
import re
import time
import subprocess
import shutil
from Bio import Entrez, SeqIO
import io

Entrez.email = "a71364@ualg.pt"

SLEEP = 0.34


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a fast phylogeny from all PHMMER input hits"
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Folder with PHMMER input files",
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="Output folder (default: outputs/phmmer_all_run_<timestamp>)",
    )
    parser.add_argument(
        "--cache",
        default=None,
        help="JSON cache for accession->sequence",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="Fetch missing sequences from NCBI (slower)",
    )
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=500,
        help="Maximum number of accessions to include (default: 500)",
    )
    return parser.parse_args()


def sanitize_name(value):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return safe or "output"


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
        return f"{first[0].upper()}{last[:2].lower()}"
    except:
        return "NA"


def find_inputs(input_dir):
    if input_dir is None:
        base = os.path.join(os.path.dirname(__file__), "phmmer_inputs")
        input_dir = base
    if not os.path.isdir(input_dir):
        raise SystemExit(f"Input dir not found: {input_dir}")
    files = []
    for name in sorted(os.listdir(input_dir)):
        if name.lower().endswith((".txt", ".tsv")):
            files.append(os.path.join(input_dir, name))
    if not files:
        raise SystemExit("No PHMMER input files found.")
    return files


def collect_accessions(files):
    accs = []
    seen = set()
    for path in files:
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                acc = parts[0]
                if acc not in seen:
                    seen.add(acc)
                    accs.append(acc)
    return accs


def load_cache(path):
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path) as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except:
        return {}


def save_cache(path, cache):
    if not path:
        return
    with open(path, "w") as f:
        json.dump(cache, f)


def fetch_sequence(acc):
    try:
        h = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        seq = SeqIO.read(io.StringIO(h.read()), "fasta")
        h.close()
        return seq.description, str(seq.seq)
    except:
        return None, None


def resolve_output_root(output_root):
    if output_root:
        return output_root
    run_id = time.strftime("%Y%m%d_%H%M%S")
    return os.path.join(os.getcwd(), "outputs", f"phmmer_all_run_{run_id}")


def resolve_fasttree_cmd():
    for candidate in ["fasttree", "FastTree"]:
        if shutil.which(candidate):
            return candidate
    raise FileNotFoundError("fasttree not found in PATH")


def main():
    args = parse_args()
    files = find_inputs(args.input_dir)
    accs = collect_accessions(files)
    if args.max_seqs and len(accs) > args.max_seqs:
        accs = accs[: args.max_seqs]

    output_root = resolve_output_root(args.output_root)
    os.makedirs(output_root, exist_ok=True)

    cache_path = args.cache or os.path.join(output_root, "phmmer_all_cache.json")
    cache = load_cache(cache_path)

    fasta_path = os.path.join(output_root, "PHMMER_ALL.fasta")
    aln_path = os.path.join(output_root, "PHMMER_ALL.aln.fasta")
    tree_path = os.path.join(output_root, "PHMMER_ALL.fasttree.nwk")

    with open(fasta_path, "w") as f_out:
        total = len(accs)
        for i, acc in enumerate(accs, start=1):
            cached = cache.get(acc)
            description = None
            sequence = None
            if isinstance(cached, dict):
                description = cached.get("description")
                sequence = cached.get("sequence")
            elif isinstance(cached, str):
                sequence = cached

            if not sequence and args.fetch:
                description, sequence = fetch_sequence(acc)
                cache[acc] = {
                    "description": description or "",
                    "sequence": sequence or "",
                }
                time.sleep(SLEEP)

            if not sequence:
                continue

            abbrev = get_species_abbrev(description or "")
            label = f"{abbrev}_{acc}"
            f_out.write(f">{label}\n{sequence}\n")

            if i % 200 == 0 or i == total:
                print(f"Fetched {i}/{total}")

    save_cache(cache_path, cache)

    with open(aln_path, "w") as out_f:
        subprocess.run([
            "mafft",
            "--6merpair",
            "--retree",
            "1",
            "--maxiterate",
            "0",
            fasta_path,
        ], check=True, stdout=out_f)

    fasttree_cmd = resolve_fasttree_cmd()
    with open(tree_path, "w") as out_f:
        subprocess.run([fasttree_cmd, "-lg", "-nosupport", aln_path], check=True, stdout=out_f)

    print(f"All-hits phylogeny outputs saved to: {output_root}")
    print(f"Tree file: {tree_path}")


if __name__ == "__main__":
    main()
