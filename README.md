# GPCR B1 Family Pipeline

Reproducible bioinformatics pipeline to identify, classify, and phylogenetically analyze **Class B1 GPCRs** from PHMMER outputs across vertebrate species.

This repository is structured as a practical, end-to-end workflow suitable for comparative genomics and receptor evolution studies.

## Why This Project Matters

Class B1 GPCRs (e.g., `GCGR`, `GLP1R`, `CRHR1/2`, `PTH1R/2`) are key signaling receptors with broad physiological and translational relevance. This pipeline automates:

- accession parsing from PHMMER hits,
- sequence and gene-metadata retrieval from NCBI,
- receptor subfamily classification via curated keyword rules,
- per-species and combined Excel reports,
- multi-species alignment and phylogeny reconstruction,
- publication-ready tree rendering in PDF.

## Core Features

- Staged execution (`retrieve`, `classify`, `tree`) plus one-command runs (`full`, `from-cache`)
- NCBI caching to avoid repeated API calls
- Family/subfamily-aware classification logic for B1 GPCR repertoire
- Species-level and global reporting (`.xlsx`)
- MAFFT + FastTree workflow with configurable alignment modes
- Color-coded phylogenetic trees by family or subfamily

## Repository Structure

```text
gpcr-b1family-pipeline/
├── pipeline_gpcr_staged.py        # Main orchestrator (recommended entry point)
├── extract_ncbi_cache.py          # Step 1: retrieve NCBI metadata + sequences into JSON cache
├── build_excels_from_cache.py     # Step 2: classify + generate Excel reports
├── build_phylogeny.py             # Step 3: FASTA/alignment/tree generation
├── render_tree_pdf.py             # Tree-to-PDF renderer with color modes
├── master_pipeline_gpcr.py        # Legacy monolithic script
├── build_phylogeny_phmmer_all.py  # Alternate phylogeny helper
├── phmmer_inputs/                 # Input PHMMER result files (.txt/.tsv)
├── requirements.txt
└── README.md
```

## Workflow

1. **Retrieve**
Fetch GeneID, sequence, and basic genomic metadata for all unique accessions found in PHMMER files.

2. **Classify + Report**
Classify selected hits into B1 GPCR subfamilies and groups (`CALCR`, `CRHR`, `GCGR`, `SCTR`, `PTHR`), then generate:
- per-species workbooks
- one combined workbook

3. **Phylogeny**
Aggregate classified sequences, align with MAFFT, infer tree with FastTree, and render annotated PDFs.

## Installation

### 1. Python dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. External tools
Install and expose these binaries in your `PATH`:

- `mafft`
- `fasttree` (or `FastTree`)

## Usage

### A) Full pipeline (recommended)

```bash
python pipeline_gpcr_staged.py full \
  --cache-out outputs/cache_reconstructed_from_excels.json \
  --email your.email@domain.com \
  --mafft-mode auto \
  --prefix ALL_CLASSIFIED_AUTOHQ
```

### B) Start from an existing cache (faster reruns)

```bash
python pipeline_gpcr_staged.py from-cache \
  --cache outputs/cache_reconstructed_from_excels.json \
  --output-root outputs/run_from_cache \
  --mafft-mode auto \
  --prefix ALL_CLASSIFIED_AUTOHQ
```

### C) Run steps independently

```bash
# Step 1
python pipeline_gpcr_staged.py retrieve \
  --input-dir phmmer_inputs \
  --cache-out outputs/cache.json \
  --email your.email@domain.com

# Step 2
python pipeline_gpcr_staged.py classify \
  --cache outputs/cache.json \
  --output-root outputs/run_manual

# Step 3
python pipeline_gpcr_staged.py tree \
  --run-root outputs/run_manual \
  --mafft-mode auto \
  --prefix ALL_CLASSIFIED_AUTOHQ
```

## Outputs

Expected artifacts inside each run folder:

- `*_FamilyB1.xlsx` (per species)
- `All_species_combined.xlsx`
- `phylogeny_all/*.fasta`
- `phylogeny_all/*.aln.fasta`
- `phylogeny_all/*.fasttree.nwk`
- `phylogeny_all/*.family_colors.pdf`
- `phylogeny_all/*.subfamily_colors.pdf`

## Technical Notes

- Entrez calls include a configurable throttle (`--sleep`) to reduce request pressure.
- Selection logic keeps one representative hit per gene (longest sequence + accession preference rules).
- Unclassified records are retained and can be partially recovered during phylogeny pre-processing.

## Reproducibility & Good Practice

- Keep raw PHMMER files versioned in `phmmer_inputs/`.
- Store cache JSONs to ensure deterministic reruns.
- Record exact MAFFT mode and run prefix for each analysis batch.

## Career Context

This project demonstrates practical skills expected in remote bioinformatics roles:

- Python pipeline engineering for biological data
- API-driven data collection (NCBI/Entrez)
- data wrangling and reporting automation
- phylogenetic workflow orchestration
- reproducible, scriptable computational biology

## Author

**Joao Barros**  
GitHub: [JB5000](https://github.com/JB5000)
