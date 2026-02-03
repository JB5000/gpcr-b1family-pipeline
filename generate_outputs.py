#!/usr/bin/env python3
"""
Generate Excel workbook with multiple sheets and combined FASTA file
"""
import os
import re
import argparse
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils.dataframe import dataframe_to_rows

DEFAULT_OUTPUT_DIR = "FINAL_OUTPUT"

# Define subfamilies in order
SUBFAMILIES = [
    "SCTR", "GCGR", "GIPR", "GLP1R", "GLP2R", "GHRHR",
    "PTH1R", "PTH2R", "CRHR1", "CRHR2", "ADCYAP1R1",
    "VIPR1", "VIPR2", "CALCR", "CALCRL", "UNCLASSIFIED"
]

def sanitize_name(value):
    safe = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return safe or "output"


def resolve_output_dir(input_file):
    base = os.path.splitext(os.path.basename(input_file))[0]
    safe_base = sanitize_name(base)
    return f"FINAL_OUTPUT_{safe_base}"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate Excel + FASTA outputs for GPCR B1 pipeline"
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        default=None,
        help="PHMMER input file used in the pipeline",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Pipeline output directory (overrides --input)",
    )
    return parser.parse_args()


args = parse_args()
if args.output_dir:
    OUTPUT_DIR = args.output_dir
elif args.input_file:
    OUTPUT_DIR = resolve_output_dir(args.input_file)
else:
    OUTPUT_DIR = DEFAULT_OUTPUT_DIR

EXCEL_OUTPUT = os.path.join(OUTPUT_DIR, "GPCR_B1_All_Subfamilies.xlsx")
FASTA_OUTPUT = os.path.join(OUTPUT_DIR, "GPCR_B1_Classified_Sequences.fasta")

print("🔄 Starting output generation...")
print(f"📊 Creating Excel workbook: {EXCEL_OUTPUT}")
print(f"🧬 Creating FASTA file: {FASTA_OUTPUT}")

# Create Excel workbook
wb = Workbook()
wb.remove(wb.active)  # Remove default sheet

fasta_sequences = []
total_classified = 0

for subfamily in SUBFAMILIES:
    tsv_path = os.path.join(OUTPUT_DIR, subfamily, f"{subfamily}.tsv")
    
    if not os.path.exists(tsv_path):
        print(f"⚠️  Skipping {subfamily} - file not found")
        continue
    
    # Read TSV
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Create sheet in Excel
    ws = wb.create_sheet(title=subfamily)
    
    # Add header with formatting
    header_fill = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
    header_font = Font(bold=True, color="FFFFFF")
    
    for r_idx, row in enumerate(dataframe_to_rows(df, index=False, header=True), 1):
        for c_idx, value in enumerate(row, 1):
            cell = ws.cell(row=r_idx, column=c_idx, value=value)
            
            if r_idx == 1:  # Header row
                cell.fill = header_fill
                cell.font = header_font
                cell.alignment = Alignment(horizontal="center", vertical="center")
    
    # Auto-adjust column widths
    for column in ws.columns:
        max_length = 0
        column_letter = column[0].column_letter
        for cell in column:
            try:
                if cell.value:
                    max_length = max(max_length, len(str(cell.value)))
            except:
                pass
        adjusted_width = min(max_length + 2, 50)  # Cap at 50
        ws.column_dimensions[column_letter].width = adjusted_width
    
    # Freeze header row
    ws.freeze_panes = "A2"
    
    print(f"✅ Added sheet: {subfamily} ({len(df)} sequences)")
    
    # Extract FASTA sequences (skip UNCLASSIFIED)
    if subfamily != "UNCLASSIFIED" and 'protein_sequence' in df.columns:
        for seq in df['protein_sequence'].dropna():
            if seq.startswith('>'):
                fasta_sequences.append(seq)
                total_classified += 1

# Save Excel
wb.save(EXCEL_OUTPUT)
print(f"\n📊 Excel saved: {EXCEL_OUTPUT}")

# Save FASTA
with open(FASTA_OUTPUT, 'w') as f:
    for seq in fasta_sequences:
        f.write(seq + '\n')

print(f"🧬 FASTA saved: {FASTA_OUTPUT}")
print(f"   Total classified sequences: {total_classified}")
print("\n✅ ALL DONE!")
