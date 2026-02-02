# GPCR Pipeline Modifications Summary

## Overview
The `master_pipeline_gpcr.py` has been updated to organize GPCR isoforms into 15 distinct receptor subfamily categories with individual folder structures.

## Changes Made

### 1. **Expanded FAMILIES Dictionary (Lines 23-38)**
Updated from 5 general families to 15 specific receptor categories:

- **SCTR** - Secretin receptor
- **GCGR** - Glucagon receptor
- **GIPR** - Gastric inhibitory polypeptide receptor
- **GLP1R** - Glucagon-like peptide-1 receptor
- **GLP2R** - Glucagon-like peptide-2 receptor
- **GHRHR** - Growth hormone-releasing hormone receptor
- **PTH1R** - Parathyroid hormone receptor 1
- **PTH2R** - Parathyroid hormone receptor 2 *(NEW)*
- **CRHR1** - Corticotropin-releasing hormone receptor 1
- **CRHR2** - Corticotropin-releasing hormone receptor 2
- **ADCYAP1R1** - PAC1 receptor (PACAP type I)
- **VIPR1** - Vasoactive intestinal peptide receptor 1
- **VIPR2** - Vasoactive intestinal peptide receptor 2
- **CALCR** - Calcitonin receptor
- **CALCRL** - Calcitonin receptor-like receptor

### 2. **Updated FASTA Header Format (Lines 89-93)**
Changed from: `>species_abbr accession sequence`
Changed to: `>species_abbr_accession_subfamily sequence`

Example output:
```
>Tru_NP_000151.1_GCGR MSTGWTLLTGLLLGTCGSVGVACGGCNLAKPPPPTWCCGSSSFSAPPPGFQRFRQAGP...
```

### 3. **Created Subfamily Folders (Lines 237-246)**
Each subfamily now has its own folder in `FINAL_OUTPUT/`:
```
FINAL_OUTPUT/
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

### 4. **Updated Output Logic (Lines 248-281)**
- Generates subfamily-specific FASTA headers for each classified record
- Each TSV file contains the formatted protein sequences with subfamily identifiers
- Unclassified sequences are marked with "UNKNOWN" subfamily tag

## File Format
Each subfamily TSV file contains:
```
accession	tlen	GeneID	exon_count	chromosome	chr_accver	chr_start	chr_stop	protein_sequence	description
```

Where `protein_sequence` column contains the FASTA header in the new format:
```
>species_abbr_accession_subfamily protein_sequence_data
```

## Usage
Run the modified pipeline with:
```bash
python master_pipeline_gpcr.py
```

Results will be organized in `FINAL_OUTPUT/` with each subfamily in its own directory.
