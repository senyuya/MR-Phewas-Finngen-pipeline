# MR-Phewas-Finngen-pipeline

Minimal utilities for working with **FinnGen R12** summary statistics:
(1) filter the public FinnGen R12 manifest and batch-download summary statistics from the `path_https` column.  
(2)  extract rows by a **list of rsIDs** (“Filter SNPs”),
(3) (coming soon) run **MR-PheWAS** on the filtered variants.

No GWAS data are included in this repository—only small helper scripts.

## Quickstart
### Step 1 — Download (after filtering the manifest)
### Features
- Keep phenotypes with `num_cases >= 1000` and `num_controls != 0`
- Save a filtered manifest as TSV (same columns as the original)
- Download `.gz` files listed in `path_https` (skips files that already exist)
- Ready for local runs or Slurm submission on HPC

python scripts/downloadFinngenGwas.py
or: sbatch download_finngen.sh

### Step 2 — Filter SNPs (rsID list from your MR table's column `SNP/RSID/...`)
**Goal:** Given a list of rsIDs, extract matching rows from each FinnGen summary-statistics file.
python scripts/Filter_SNPs.py     # edit paths at the top of the script
or: sbatch Filter_SNPs.sh

### Step 3 — Run MR-PheWAS (placeholder; to be added)
### Rscript scripts/run_mr_phewas.R --inputs ...  (TBD)

---

## Requirements
- Python 3.9+
- Packages:
  - `pandas >= 1.3`
  - `requests >= 2.31`
  - `openpyxl >= 3.1` (only needed if you read Excel `.xlsx`)
- Sufficient disk space for large GWAS files
