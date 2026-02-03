#!/usr/bin/env python3
"""
Script para corrigir genes mal classificados e prevenir futuros erros.

PROBLEMA IDENTIFICADO:
- Keyword genérica "vasoactive intestinal polypeptide receptor" em VIPR1
  captura também VIPR2 (que contém "receptor 2")
- Solução: Priorizar keywords específicas (com números) antes das genéricas

GENES PARA CORRIGIR:
1. XP_003968215.1 - VIPR1 → VIPR2 (Takifugu rubripes)
2. NP_571854.1 - VIPR1 → VIPR2 (Danio rerio)  
3. XP_029687315.1 - VIPR1 → VIPR2 (Takifugu rubripes)
"""

import os
import shutil
from pathlib import Path

# Configuração
BASE_DIR = Path(__file__).parent
OUTPUTS = [
    "FINAL_OUTPUT_Galaxy14_PHMMER_Takifugu_rubripes",
    "FINAL_OUTPUT_Galaxy28_PHMMER_Danio_rerio",
    "FINAL_OUTPUT_Galaxy23_PHMMER_Mytilus_galloprovincialis",
]

# Genes a corrigir (AccessionID: nova_subfamília)
CORRECTIONS = {
    "XP_003968215.1": "VIPR2",
    "NP_571854.1": "VIPR2",
    "XP_029687315.1": "VIPR2",
}

def fix_classification():
    """Move genes mal classificados para subfamílias corretas"""
    
    print("=" * 80)
    print("CORREÇÃO DE GENES MAL CLASSIFICADOS")
    print("=" * 80)
    
    for output_dir in OUTPUTS:
        output_path = BASE_DIR / output_dir
        if not output_path.exists():
            continue
            
        print(f"\n📁 Processando: {output_dir}")
        
        # Ler FASTA principal
        fasta_path = output_path / "GPCR_B1_Classified_Sequences.fasta"
        if not fasta_path.exists():
            print(f"  ⚠ FASTA não encontrado")
            continue
            
        with open(fasta_path, 'r') as f:
            fasta_lines = f.readlines()
        
        # Processar cada gene
        for acc_id, correct_family in CORRECTIONS.items():
            # Procurar no FASTA
            seq_start = None
            for i, line in enumerate(fasta_lines):
                if acc_id in line:
                    seq_start = i
                    header = line.strip()
                    # Extrair família atual do header
                    current_family = header.split('_')[-1]
                    
                    if current_family == correct_family:
                        print(f"  ✓ {acc_id}: Já está em {correct_family}")
                        break
                    
                    print(f"  🔄 {acc_id}: {current_family} → {correct_family}")
                    
                    # Atualizar header no FASTA
                    organism = header.split('_')[0].replace('>', '')
                    new_header = f">{organism}_{acc_id}_{correct_family}"
                    fasta_lines[i] = new_header + '\n'
                    
                    # Mover arquivos TSV entre pastas
                    old_tsv = output_path / current_family / f"{current_family}.tsv"
                    new_tsv = output_path / correct_family / f"{correct_family}.tsv"
                    
                    if old_tsv.exists() and new_tsv.exists():
                        # Ler TSV antigo e encontrar linha do gene
                        with open(old_tsv, 'r') as f:
                            old_lines = f.readlines()
                        
                        gene_line = None
                        new_old_lines = []
                        for line in old_lines:
                            if acc_id in line:
                                gene_line = line
                            else:
                                new_old_lines.append(line)
                        
                        if gene_line:
                            # Remover do TSV antigo
                            with open(old_tsv, 'w') as f:
                                f.writelines(new_old_lines)
                            
                            # Adicionar ao TSV novo
                            with open(new_tsv, 'a') as f:
                                f.write(gene_line)
                            
                            print(f"    ✓ TSV atualizado")
                    
                    break
        
        # Reescrever FASTA atualizado
        with open(fasta_path, 'w') as f:
            f.writelines(fasta_lines)
        
        print(f"  ✓ FASTA atualizado")

def print_keyword_rules():
    """Documenta regras para evitar futuros erros"""
    
    print("\n" + "=" * 80)
    print("REGRAS PARA KEYWORDS (para evitar erros futuros)")
    print("=" * 80)
    print("""
⚠️ KEYWORDS PROBLEMÁTICAS:

1. NÃO usar keywords genéricas que capturam subfamílias específicas:
   ❌ "vasoactive intestinal polypeptide receptor" (sem número)
      → Captura VIPR1 E VIPR2!
   
2. SEMPRE especificar número do receptor quando existir subfamílias:
   ✅ "vasoactive intestinal peptide receptor 1"
   ✅ "vasoactive intestinal polypeptide receptor 1"
   ✅ "vipr1"
   
3. ORDEM DE PRIORIDADE nas keywords:
   1º) Keywords com números específicos (receptor 1, receptor 2)
   2º) Keywords genéricas (apenas se não houver subfamílias numeradas)
   
4. SUBFAMÍLIAS COM NÚMEROS (verificar keywords):
   - VIPR1 / VIPR2
   - GLP1R / GLP2R
   - PTH1R / PTH2R
   - CRHR1 / CRHR2
   
5. SOLUÇÃO NO CÓDIGO:
   - Verificar PRIMEIRO se descrição contém "receptor 2" ou "receptor-2"
   - Só depois verificar keyword genérica
   - OU: Remover keywords genéricas de subfamílias numeradas

EXEMPLO DE CORREÇÃO NO master_pipeline_gpcr.py:

    "VIPR1": [
        "vasoactive intestinal peptide receptor 1",
        "vasoactive intestinal polypeptide receptor 1",
        # "vasoactive intestinal polypeptide receptor",  ← REMOVER!
        "vipr1",
    ],
    "VIPR2": [
        "vasoactive intestinal peptide receptor 2",
        "vasoactive intestinal polypeptide receptor 2",
        "vipr2",
    ],
""")

if __name__ == "__main__":
    fix_classification()
    print_keyword_rules()
    print("\n✅ Correções aplicadas! Agora execute o pipeline novamente para verificar.")
