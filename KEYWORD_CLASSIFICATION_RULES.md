# KEYWORDS CLASSIFICATION RULES - GPCR B1 Family Pipeline

## ⚠️ PROBLEMA IDENTIFICADO (03/02/2026)

Três genes foram **mal classificados** como VIPR1 quando deveriam ser VIPR2:
- `XP_003968215.1` (Takifugu rubripes) - "vasoactive intestinal polypeptide receptor 2"
- `NP_571854.1` (Danio rerio) - "vasoactive intestinal polypeptide receptor 2"  
- `XP_029687315.1` (Takifugu rubripes) - "vasoactive intestinal polypeptide receptor 2"

### Causa Raiz
A keyword genérica **"vasoactive intestinal polypeptide receptor"** (sem número) estava em `VIPR1`, capturando também sequências que continham "receptor 2".

### Solução Aplicada
Removida a keyword genérica de `VIPR1` no arquivo `master_pipeline_gpcr.py`:

```python
# ANTES (❌ ERRADO):
"VIPR1": [
    "vasoactive intestinal peptide receptor 1",
    "vasoactive intestinal polypeptide receptor 1",
    "vasoactive intestinal polypeptide receptor",  # ← PROBLEMA!
    "vipr1",
],

# DEPOIS (✅ CORRETO):
"VIPR1": [
    "vasoactive intestinal peptide receptor 1",
    "vasoactive intestinal polypeptide receptor 1",
    # REMOVED: generic keyword that captured VIPR2
    "vipr1",
],
```

---

## 📋 REGRAS PARA KEYWORDS

### 1. **NUNCA usar keywords genéricas em subfamílias numeradas**
   
❌ **ERRADO**:
```python
"VIPR1": ["vasoactive intestinal polypeptide receptor", ...]  # Captura VIPR2!
"GLP1R": ["glucagon-like peptide receptor", ...]               # Capturaria GLP2R!
"PTH1R": ["parathyroid hormone receptor", ...]                 # Capturaria PTH2R!
```

✅ **CORRETO**:
```python
"VIPR1": ["vasoactive intestinal polypeptide receptor 1", "vipr1"]
"VIPR2": ["vasoactive intestinal polypeptide receptor 2", "vipr2"]
```

---

### 2. **Subfamílias COM números (cuidado extra)**

Estas subfamílias têm variantes numeradas e **devem** especificar o número nas keywords:

| Subfamília | Keywords OBRIGATÓRIAS | Não Usar |
|------------|----------------------|----------|
| **VIPR1/VIPR2** | "receptor 1", "receptor 2", "vipr1", "vipr2" | "vasoactive intestinal polypeptide receptor" |
| **GLP1R/GLP2R** | "glucagon-like peptide-1", "glucagon-like peptide-2", "glp1r", "glp2r" | "glucagon-like peptide receptor" |
| **PTH1R/PTH2R** | "parathyroid hormone receptor 1", "parathyroid hormone receptor 2", "pth1r", "pth2r" | "parathyroid hormone receptor" |
| **CRHR1/CRHR2** | "crhr1", "crhr2", "corticotropin releasing hormone receptor 1/2" | "corticotropin-releasing" |

---

### 3. **Ordem de prioridade na classificação**

O código deve verificar na seguinte ordem:

1. **Keywords específicas com números** (e.g., "receptor 2", "vipr2")
2. **Keywords específicas sem números** (apenas se não houver subfamílias numeradas)
3. **Keywords genéricas** (apenas para famílias sem subfamílias)

**Implementação recomendada**:
```python
def classify_sequence(description):
    desc_lower = description.lower()
    
    # 1º: Verificar PRIMEIRO subfamílias específicas (com números)
    if "receptor 2" in desc_lower or "polypeptide receptor 2" in desc_lower:
        if "vasoactive intestinal" in desc_lower:
            return "VIPR2"
        elif "glucagon-like" in desc_lower:
            return "GLP2R"
        elif "parathyroid" in desc_lower:
            return "PTH2R"
    
    # 2º: Depois verificar subfamílias com "1"
    if "receptor 1" in desc_lower or "polypeptide receptor 1" in desc_lower:
        if "vasoactive intestinal" in desc_lower:
            return "VIPR1"
        elif "glucagon-like" in desc_lower:
            return "GLP1R"
        elif "parathyroid" in desc_lower:
            return "PTH1R"
    
    # 3º: Keywords genéricas (apenas se não houver número)
    # ...
```

---

### 4. **Casos especiais - Sinônimos problemáticos**

Algumas descrições do NCBI usam variantes que podem causar confusão:

| Descrição NCBI | Parece | Na verdade é |
|----------------|--------|--------------|
| "vasoactive intestinal **polypeptide** receptor 2" | VIPR1? | **VIPR2** ✅ |
| "vasoactive intestinal **peptide** receptor 2" | VIPR1? | **VIPR2** ✅ |
| "parathyroid hormone/**parathyroid hormone-related** peptide receptor" | PTH1R ou PTH2R? | **PTH1R** ✅ (receptor híbrido) |
| "LOW QUALITY PROTEIN: parathyroid hormone..." | Ignorar? | **PTH1R** ✅ (manter mas marcar) |

**Solução**: Sempre procurar pelo **número do receptor** primeiro, independentemente da variante de "peptide/polypeptide".

---

### 5. **Verificação pós-classificação**

Após classificar todas as sequências, executar:

```bash
# Verificar se há VIPR2 classificados como VIPR1
grep "receptor 2" FINAL_OUTPUT_*/VIPR1/VIPR1.tsv

# Verificar se há GLP2R classificados como GLP1R
grep "peptide-2" FINAL_OUTPUT_*/GLP1R/GLP1R.tsv

# Verificar se há PTH2R classificados como PTH1R
grep "hormone receptor 2" FINAL_OUTPUT_*/PTH1R/PTH1R.tsv
```

---

## 🔧 COMANDOS DE CORREÇÃO

Se encontrar genes mal classificados:

```python
# 1. Identificar genes problemáticos
python3 << 'EOF'
from Bio import Entrez
Entrez.email = "a71364@ualg.pt"

# Buscar descrições do NCBI
for acc_id in ["XP_003968215.1", "NP_571854.1", "XP_029687315.1"]:
    handle = Entrez.efetch(db="protein", id=acc_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    print(f"{acc_id}: {record.description}")
EOF

# 2. Corrigir FASTA headers
sed -i 's/>Tru_XP_003968215.1_VIPR1/>Tru_XP_003968215.1_VIPR2/g' FINAL_OUTPUT_*/GPCR_B1_Classified_Sequences.fasta
sed -i 's/>Dre_NP_571854.1_VIPR1/>Dre_NP_571854.1_VIPR2/g' FINAL_OUTPUT_*/GPCR_B1_Classified_Sequences.fasta
sed -i 's/>Tru_XP_029687315.1_VIPR1/>Tru_XP_029687315.1_VIPR2/g' FINAL_OUTPUT_*/GPCR_B1_Classified_Sequences.fasta
```

---

## ✅ CHECKLIST PRÉ-EXECUÇÃO

Antes de executar o pipeline em **novos organismos**:

- [ ] Verificar se há subfamílias numeradas (VIPR1/2, GLP1R/2R, PTH1R/2R, CRHR1/2)
- [ ] Garantir que keywords **não têm** termos genéricos capturando números
- [ ] Testar com 1-2 sequências de cada subfamília numerada
- [ ] Executar verificação pós-classificação (grep nos TSVs)
- [ ] Confirmar com NCBI se houver dúvidas (Entrez.efetch)

---

## 📚 REFERÊNCIAS

- **Arquivo corrigido**: `master_pipeline_gpcr.py` (linhas 68-77)
- **Genes corrigidos**: Commitados em 03/02/2026
- **Script de correção**: `fix_misclassified_genes.py`

---

**Última atualização**: 03/02/2026  
**Autor**: Correção aplicada durante análise filogenética GPCR B1
