# GPCR B1 Family Pipeline

A comprehensive bioinformatics pipeline for processing, filtering, and annotating GPCR B1 Family protein sequences from PHMMER results. This tool extracts detailed genomic and proteomic information from the NCBI Entrez database and classifies sequences into family categories.

## 🎯 Overview

This pipeline processes HMM search results (PHMMER output) to:
- **Identify and deduplicate** protein sequences
- **Retrieve genomic locations** and exon information from NCBI Gene database
- **Fetch protein sequences** for detailed analysis
- **Classify sequences** into GPCR B1 family categories (CALCR, CRHR, PTHR, GCGR, VIP_SCTR)
- **Generate annotated TSV files** with comprehensive genomic and proteomic data

## 🔧 What This Pipeline Can Do

### 1. **Sequence Processing**
- Reads PHMMER HMM search output files
- Deduplicates protein accessions (retains unique sequences only)
- Groups sequences by Gene ID
- Filters to the longest/best representative sequence per gene
- Prefers NP (RefSeq) accessions over XP (model) accessions

### 2. **NCBI Data Retrieval**
- **Gene ID Mapping**: Links protein accessions to NCBI Gene database entries
- **Exon Count**: Extracts the number of exons for each gene
- **Chromosomal Coordinates**: Retrieves:
  - Chromosome location
  - Chromosome accession version
  - Genomic start position (bp)
  - Genomic stop position (bp)
- **Protein Sequences**: Downloads full-length amino acid sequences in FASTA format

### 3. **Family Classification**
Automatically categorizes sequences into 15 GPCR B1 receptor subfamilies, each in its own folder:
- **SCTR** - Secretin receptor
- **GCGR** - Glucagon receptor
- **GIPR** - Gastric inhibitory peptide receptor
- **GLP1R** - Glucagon-like peptide-1 receptor
- **GLP2R** - Glucagon-like peptide-2 receptor
- **GHRHR** - Growth hormone-releasing hormone receptor
- **PTH1R** - Parathyroid hormone receptor 1
- **PTH2R** - Parathyroid hormone receptor 2
- **CRHR1** - Corticotropin-releasing hormone receptor 1
- **CRHR2** - Corticotropin-releasing hormone receptor 2
- **ADCYAP1R1** - PAC1 receptor (PACAP type I)
- **VIPR1** - Vasoactive intestinal peptide receptor 1
- **VIPR2** - Vasoactive intestinal peptide receptor 2
- **CALCR** - Calcitonin receptor
- **CALCRL** - Calcitonin receptor-like receptor
- **UNCLASSIFIED** - Sequences that don't match any subfamily keywords

### 4. **Output Generation**
Produces tab-separated value (TSV) files with the following columns:

| Column | Description |
|--------|-------------|
| `accession` | NCBI protein accession (e.g., NP_001098689.1) |
| `tlen` | Target sequence length (amino acids) |
| `GeneID` | NCBI Gene ID |
| `exon_count` | Number of exons in the gene |
| `chromosome` | Chromosome location (e.g., 12, X, MT) |
| `chr_accver` | Chromosome sequence accession version (e.g., NC_042296.1) |
| `chr_start` | Genomic start position (base pairs) |
| `chr_stop` | Genomic stop position (base pairs) |
| `protein_sequence` | Full amino acid sequence with custom header (see below) |
| `description` | Original PHMMER description |

## 📋 Input Requirements

### Input File Format
- **Format**: PHMMER output (tab-delimited with comment lines)
- **Structure**: Standard HMMER3 table output with sequence headers as first column
- **Expected columns**:
  - Column 1: Protein accession
  - Column 3: Target length (tlen)
  - Column 20+: Description text

### Example Input (Galaxy14_PHMMER.txt)
```
# Comment line (skipped)
NP_001098689.1       -            794 FAMILYB1    ...  calcitonin receptor precursor [Takifugu rubripes]
NP_001098689.1       -            794 FAMILYB1    ...  calcitonin receptor precursor [Takifugu rubripes]
XP_029700367.1       -            794 FAMILYB1    ...  calcitonin receptor isoform X1 [Takifugu rubripes]
```

### 🧬 Sequence Header Format (Beginner Friendly)

Each output sequence is written in a **single field** with a custom header and a space before the sequence:

```
>Tru_NP_001098689.1_CALCR MKSAGYSCILWLLLMMVTDTESLSEPSLSPGQ...
```

**How the header is built (automatic):**
- `>` symbol first
- **Species abbreviation** derived from the organism in square brackets `[]`
  - Example organism: `[Takifugu rubripes]`
  - Abbreviation rule: **first letter of genus (uppercase)** + **first two letters of species (lowercase)**
  - Result: `Tru`
- Then an underscore `_`
- Then the **accession number**
- Then an underscore `_`
- Then the **subfamily symbol** (e.g., CALCR, CRHR1, GLP1R)
- Then a **single space**
- Then the **protein sequence**

If the organism name is missing or invalid, the abbreviation becomes `NA`. Unknown subfamilies are marked as `UNKNOWN`.

---

## 🌐 How to Obtain the PHMMER Input (Galaxy.eu) — Step by Step

This pipeline expects a **PHMMER/HMMER output file** with per-domain hits. The easiest way to generate it is with **Galaxy.eu**:

### ✅ Step 1 — Open HMMER Search Tool

Use this direct link:
https://usegalaxy.eu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fhmmer_hmmsearch%2Fhmmer_hmmsearch%2F3.4%2Bgalaxy0&version=latest

Or on Galaxy.eu:
1. Go to https://usegalaxy.eu/
2. Search for **hmmsearch** in the tool search bar
3. Open the tool named **HMMER hmmsearch**

### ✅ Step 2 — Upload Your Data

You need two inputs:
1. **Reference proteome** (best: curated RefSeq animal species proteome)
2. **HMM model** (your GPCR B1 family HMM built previously)

Upload both files to Galaxy:
- Click **Upload Data**
- Select the proteome FASTA and your HMM file

### ✅ Step 3 — Configure hmmsearch

Set inputs as follows:
- **HMM file**: your GPCR B1 family HMM
- **Sequence database**: your proteome FASTA

### ✅ Step 4 — Choose Output Type

In **Output options**, select:
- ✅ **Per-domain hits** (this is important!)

This output contains the fields used by the pipeline.

### ✅ Step 5 — Run and Download

1. Click **Run Tool**
2. When finished, download the **per-domain hits** file

Rename it (example):
```
Galaxy14_PHMMER.txt
```

### ✅ Step 6 — Move File to the Pipeline Folder

Place the downloaded file in the same folder as `master_pipeline_gpcr.py`:

```
gpcr-b1family-pipeline/
├── master_pipeline_gpcr.py
├── Galaxy14_PHMMER.txt   ← place the file here
```

Now you can run the pipeline.

## 🚀 Installation & Usage

## ✅ Easy to Install (Minimal Steps)

**Linux/macOS (3 commands):**
1. Clone:
  - `git clone https://github.com/JB5000/gpcr-b1family-pipeline.git`
2. Enter folder:
  - `cd gpcr-b1family-pipeline`
3. Run setup:
  - `./setup.sh`

**Windows (3 commands):**
1. Clone:
  - `git clone https://github.com/JB5000/gpcr-b1family-pipeline.git`
2. Enter folder:
  - `cd gpcr-b1family-pipeline`
3. Run setup:
  - `setup.bat`

**Run the pipeline (after setup):**
1. Edit `master_pipeline_gpcr.py` and set your email:
  - `Entrez.email = "your_email@example.com"`
2. Put your PHMMER file in the project folder (example: `Galaxy14_PHMMER_Takifugu_rubripes.txt`)
3. Run:
  - `python master_pipeline_gpcr.py`

**Generate Excel + FASTA outputs:**
  - `python generate_outputs.py`

That’s it. The results will be in `FINAL_OUTPUT/` plus:
- `GPCR_B1_All_Subfamilies.xlsx`
- `GPCR_B1_Classified_Sequences.fasta`

### ✅ Copy‑Paste Quick Start

**Linux/macOS (all in one):**
```bash
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git && \
cd gpcr-b1family-pipeline && \
./setup.sh && \
python master_pipeline_gpcr.py && \
python generate_outputs.py
```

**Windows (PowerShell):**
```powershell
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline
setup.bat
python master_pipeline_gpcr.py
python generate_outputs.py
```

### Requirements
- **Python 3.7+** (Download: https://www.python.org/downloads/)
- **Internet connection** (needed for NCBI queries)
- **BioPython** (installed automatically - just one command!)

### ⚡ Quick Start (Automatic Setup)

#### Linux / macOS
```bash
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline
chmod +x setup.sh
./setup.sh
```

#### Windows
```bash
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline
setup.bat
```

That's it! The scripts handle everything. ✅

### 📚 Detailed Setup Instructions

**Need more details or having issues?** 

👉 See [SETUP_GUIDE.md](SETUP_GUIDE.md) for **complete step-by-step instructions** including:
- Manual setup for any OS
- Troubleshooting common problems
- Understanding virtual environments
- Verification checklists

### Manual Setup (If you prefer)

```bash
# 1. Clone repository
git clone https://github.com/JB5000/gpcr-b1family-pipeline.git
cd gpcr-b1family-pipeline

# 2. Create virtual environment
python -m venv venv

# 3. Activate virtual environment
source venv/bin/activate              # Linux/macOS
# OR on Windows:
venv\Scripts\activate

# 4. Install dependencies
pip install -r requirements.txt
```

### Configuration

Edit the constants at the top of `master_pipeline_gpcr.py`:

```python
Entrez.email = "your_email@example.com"  # ⚠️ REQUIRED! NCBI needs this
INPUT_FILE = "Galaxy14_PHMMER.txt"       # Your PHMMER input file
OUTPUT_DIR = "FINAL_OUTPUT"              # Output directory (created if needed)
SLEEP = 0.34                             # Delay between NCBI requests (seconds)
PROGRESS_STEP = 10                       # Progress bar update frequency (%)
```

⚠️ **Important**: The `Entrez.email` field is **mandatory**. NCBI uses it to contact users if there are issues with their queries.

### Running the Pipeline

```bash
python master_pipeline_gpcr.py
```

### Expected Output

The pipeline will:
1. **Progress Indicator**: Show real-time progress for protein fetching and gene annotation
2. **Output Structure**: Create organized folder structure in `FINAL_OUTPUT/`:
   - One folder per subfamily (SCTR, GCGR, GIPR, GLP1R, GLP2R, GHRHR, PTH1R, PTH2R, CRHR1, CRHR2, ADCYAP1R1, VIPR1, VIPR2, CALCR, CALCRL)
   - Each folder contains a TSV file with the subfamily name
   - Unclassified sequences in `UNCLASSIFIED/UNCLASSIFIED.tsv`

```
Total unique accessions: 240

[10%] 24/240 | elapsed: 0:00:50 | ETA: 0:07:33
[20%] 48/240 | elapsed: 0:02:19 | ETA: 0:09:16
...
[100%] 240/240 | elapsed: 0:11:32 | ETA: 0:00:00

Total unique GeneIDs: 66

[Gene 10%] 7/66 | elapsed: 0:00:06 | ETA: 0:00:51
...
[Gene 100%] 66/66 | elapsed: 0:00:49 | ETA: 0:00:00

PIPELINE COMPLETE ✅
Results written to: FINAL_OUTPUT/
```

## 📊 Analysis Capabilities

### Genomic Analysis
- **Exon Structure**: Understand gene complexity through exon counts
- **Chromosomal Distribution**: Map receptor genes to chromosomal locations
- **Comparative Genomics**: Identify orthologs across organisms

### Proteomic Analysis
- **Sequence Alignment**: Use downloaded sequences for BLAST, alignment studies
- **Domain Analysis**: Perform transmembrane domain prediction (7-TM characteristic of GPCRs)
- **Evolutionary Studies**: Track sequence divergence and family relationships

### Quality Control
- **Deduplication**: Removes redundant sequences (same protein from different sources)
- **Best Representative**: Selects longest/highest-quality sequence per gene
- **Preference Order**: Prioritizes RefSeq (NP_) over predicted (XP_) sequences

## 🔬 Use Cases

1. **GPCR Classification**: Automatically categorize newly found receptors into subfamilies
2. **Genomic Mapping**: Create comparative maps of receptor genes across species
3. **Phylogenetic Studies**: Use sequences and genomic positions for evolutionary analysis
4. **Drug Target Identification**: Compile comprehensive receptor protein databases
5. **Structural Prediction**: Use sequences for transmembrane topology and structure modeling
6. **Exon-Intron Analysis**: Study splicing patterns and alternative isoforms

## 🛠️ Customization

### Adding New Subfamilies or Refining Keywords

Edit the `FAMILIES` dictionary in `master_pipeline_gpcr.py`:

```python
FAMILIES = {
    "SCTR": ["secretin receptor", "secretin", "sctr"],
    "GCGR": ["glucagon receptor", "gcgr"],
    "GLP1R": [
        "glucagon-like peptide-1",
        "glucagon like peptide 1",
        "glp1r",
    ],
    # ... add more as needed
    "NEW_SUBFAMILY": ["keyword1", "keyword2", "alternative_name"],  # Add here
}
```

The pipeline searches the description field (case-insensitive) for any matching keyword. Add more keyword variants to improve classification accuracy.

### Adjusting NCBI Query Rates

Modify `SLEEP` parameter (default: 0.34 seconds between requests):

```python
SLEEP = 0.34  # Friendly rate: ~3 requests/second
SLEEP = 0.5   # Conservative rate: ~2 requests/second
SLEEP = 0.1   # Aggressive rate: ~10 requests/second (use with caution)
```

⚠️ **Note**: NCBI recommends 3 requests per second maximum for responsible API usage.

## 📈 Performance Metrics

- **Typical Runtime**: ~15 minutes for 240 unique accessions + 66 gene summaries
- **NCBI Rate**: Respects NCBI guidelines with configurable sleep delays
- **Data Retrieval**:
  - Protein sequences: ~0.05 seconds each
  - Gene summaries: ~0.05 seconds each
  - Gene mapping: ~0.05 seconds each

## ⚠️ Important Notes

1. **NCBI Email Required**: The `Entrez.email` field is mandatory per NCBI policy. They use it to contact users if there are issues.

2. **Network Dependency**: All data is fetched in real-time from NCBI. Internet connectivity is required.

3. **Rate Limiting**: The pipeline includes intelligent retry handling and respects NCBI's service policies.

4. **Error Handling**: Failed queries (network issues, invalid IDs) return "NA" in output rather than crashing.

5. **Data Quality**: Not all accessions will have complete information. Some organisms or genes may not have full genomic annotation.

## 📁 File Structure

```
gpcr-b1family-pipeline/
├── master_pipeline_gpcr.py         # Main pipeline script
├── Galaxy14_PHMMER.txt             # Example input file
├── Galaxy23_PHMMER.txt             # Example input file
├── README.md                        # This file
├── requirements.txt                 # Python dependencies
├── setup.sh                         # Linux/macOS setup script
├── setup.bat                        # Windows setup script
├── SETUP_GUIDE.md                   # Detailed setup instructions
└── FINAL_OUTPUT/                    # Output directory (created on first run)
    ├── SCTR/
    │   └── SCTR.tsv
    ├── GCGR/
    │   └── GCGR.tsv
    ├── GIPR/
    │   └── GIPR.tsv
    ├── GLP1R/
    │   └── GLP1R.tsv
    ├── GLP2R/
    │   └── GLP2R.tsv
    ├── GHRHR/
    │   └── GHRHR.tsv
    ├── PTH1R/
    │   └── PTH1R.tsv
    ├── PTH2R/
    │   └── PTH2R.tsv
    ├── CRHR1/
    │   └── CRHR1.tsv
    ├── CRHR2/
    │   └── CRHR2.tsv
    ├── ADCYAP1R1/
    │   └── ADCYAP1R1.tsv
    ├── VIPR1/
    │   └── VIPR1.tsv
    ├── VIPR2/
    │   └── VIPR2.tsv
    ├── CALCR/
    │   └── CALCR.tsv
    ├── CALCRL/
    │   └── CALCRL.tsv
    └── UNCLASSIFIED/
        └── UNCLASSIFIED.tsv
```

## 🔍 Troubleshooting

### Common Issues

**Problem**: `ModuleNotFoundError: No module named 'Bio'`
```bash
pip install biopython
```

**Problem**: NCBI connection timeouts
- Increase `SLEEP` value
- Check internet connection
- Verify Entrez email is correct

**Problem**: Some genes have "NA" values
- Not all organisms have complete genomic annotation in NCBI
- Some accessions may be obsolete or withdrawn
- Check NCBI Gene database directly for manual verification

**Problem**: Pipeline runs very slowly
- Reduce PROGRESS_STEP for less frequent updates
- Increase SLEEP if you have network issues (but respectfully)
- Consider running on a faster internet connection

## 📚 References

- **NCBI Entrez API**: https://www.ncbi.nlm.nih.gov/books/NBK25499/
- **BioPython Documentation**: https://biopython.org/
- **GPCR Database (GPCRdb)**: https://gpcrdb.org/
- **PHMMER**: https://www.ebi.ac.uk/Tools/phmmer/

## 📄 Citation

If you use this pipeline in your research, please cite:

```
GPCR B1 Family Pipeline v1.0
Github: https://github.com/JB5000/gpcr-b1family-pipeline
```

## 📝 License

[Specify your license here - MIT, GPL, etc.]

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push and create a pull request

## 📧 Contact & Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Contact: a71364@ualg.pt

---

**Last Updated**: February 2026  
**Version**: 2.0  
**Status**: Production Ready ✅  
**Major Update**: Now includes 15 distinct GPCR B1 subfamily categories with organized folder structure and improved FASTA header formatting.
