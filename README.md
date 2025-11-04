# MR-Phewas-Finngen-pipeline

Minimal utilities for working with **FinnGen R12** summary statistics:
1) filter the public FinnGen R12 manifest and batch-download sumstats  
2) extract rows by a **list of rsIDs** (“Filter SNPs”)  
3) (coming soon) run **MR-PheWAS** on the filtered variants

No GWAS data are included in this repository — only small helper scripts.

---

## Quickstart

### Step 1 — Download (after filtering the manifest)

**Features**
- keep phenotypes with `num_cases >= 1000` and `num_controls != 0`  
- save a filtered manifest as TSV (same columns as original)  
- download `.gz` from `path_https` (skip existing files)  
- works locally or via Slurm on HPC

**Run**
    
    python scripts/downloadFinngenGwas.py

or（HPC）：
    
    sbatch download_finngen.sh

---

### Step 2 — Filter SNPs (rsID list from your MR table’s first column `SNP`)

**Goal**  
Given a list of rsIDs, extract matching rows from each FinnGen summary-statistics file.

**Run**
    
    python scripts/Filter_SNPs.py # edit paths at the top of the script

or（HPC）：
    
    sbatch Filter_SNPs.sh     


### Step 3 — Run MR-PheWAS 

### Step 4 — Make Figure

---

## Requirements
- Python 3.9+
- Packages:
  - `pandas >= 1.3`
  - `requests >= 2.31`
  - `openpyxl >= 3.1` (only needed if you read Excel `.xlsx`)
- Sufficient disk space for large GWAS files
