#!/usr/bin/env python3

import argparse
import os
import subprocess
import csv
import matplotlib.pyplot as plt

# ----------------------------
# Argumentos
# ----------------------------
parser = argparse.ArgumentParser(
    description="Local automated BLASTP pipeline (genome vs GPCR isoforms)"
)

parser.add_argument("--genome_fasta", required=True, help="Proteome FASTA (genome proteins)")
parser.add_argument("--query_fasta", required=True, help="FASTA with query sequences (GPCRs)")
parser.add_argument("--outdir", default="blast_output", help="Output directory")
parser.add_argument("--evalue", default="1e-10", help="E-value cutoff")
parser.add_argument("--bitscore_threshold", type=float, default=300.0, help="Bitscore threshold")
parser.add_argument("--plot", action="store_true", help="Generate plot")

args = parser.parse_args()

# ----------------------------
# Setup
# ----------------------------
os.makedirs(args.outdir, exist_ok=True)
db_name = os.path.join(args.outdir, "blast_db")
results_dir = os.path.join(args.outdir, "results")
queries_dir = os.path.join(args.outdir, "queries")

os.makedirs(results_dir, exist_ok=True)
os.makedirs(queries_dir, exist_ok=True)

# ----------------------------
# 1. Criar BLAST database
# ----------------------------
print("[1/5] Creating BLAST database...")

subprocess.run([
    "makeblastdb",
    "-in", args.genome_fasta,
    "-dbtype", "prot",
    "-out", db_name
], check=True)

# ----------------------------
# 2. Separar FASTA de queries
# ----------------------------
print("[2/5] Splitting query FASTA...")

query_files = []
current_file = None

with open(args.query_fasta) as f:
    for line in f:
        if line.startswith(">"):
            if current_file:
                current_file.close()
            name = line[1:].split()[0]
            filepath = os.path.join(queries_dir, f"{name}.faa")
            query_files.append(filepath)
            current_file = open(filepath, "w")
        if current_file:
            current_file.write(line)

if current_file:
    current_file.close()

# ----------------------------
# 3. BLASTP
# ----------------------------
print("[3/5] Running BLASTP...")

blast_outputs = []

for q in query_files:
    out = os.path.join(
        results_dir,
        os.path.basename(q).replace(".faa", ".blast")
    )
    blast_outputs.append(out)

    subprocess.run([
        "blastp",
        "-query", q,
        "-db", db_name,
        "-evalue", args.evalue,
        "-outfmt", "6",
        "-out", out
    ], check=True)

# ----------------------------
# 4. Parse BLAST results
# ----------------------------
print("[4/5] Parsing results...")

summary_csv = os.path.join(args.outdir, "summary.csv")

isoforms = []
best_bitscores = []
num_strong_hits = []

with open(summary_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Isoform", "Best_bitscore", "Strong_hits"])

    for bf in blast_outputs:
        isoform = os.path.basename(bf).replace(".blast", "")
        best = 0.0
        strong = 0

        with open(bf) as f:
            for line in f:
                cols = line.split()
                bitscore = float(cols[11])

                if bitscore > best:
                    best = bitscore
                if bitscore >= args.bitscore_threshold:
                    strong += 1

        isoforms.append(isoform)
        best_bitscores.append(best)
        num_strong_hits.append(strong)

        writer.writerow([isoform, best, strong])

# ----------------------------
# 5. Plot
# ----------------------------
if args.plot:
    print("[5/5] Generating plot...")

    plt.figure(figsize=(14, 6))
    plt.bar(isoforms, best_bitscores)
    plt.axhline(
        args.bitscore_threshold,
        color="red",
        linestyle="--",
        linewidth=1,
        label="Strong GPCR threshold"
    )

    plt.xlabel("GPCR isoforms")
    plt.ylabel("Best BLAST bitscore")
    plt.title("GPCR conservation across genome")
    plt.xticks(rotation=90)
    plt.legend()
    plt.tight_layout()

    plt.savefig(os.path.join(args.outdir, "gpcr_conservation.png"), dpi=300)
    plt.show()

print("Pipeline finished successfully.")
