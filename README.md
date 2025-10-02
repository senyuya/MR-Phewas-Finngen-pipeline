# MR-Phewas-Finngen-pipeline
## FinnGen R12 GWAS Downloader 

Minimal utilities to 
(1) filter the public FinnGen R12 manifest and (2) batch-download summary statistics from the `path_https` column.  
This repository **does not** include any GWAS data—only small helper scripts.

---

## Features
- Keep phenotypes with `num_cases >= 1000` and `num_controls != 0`
- Save a filtered manifest as TSV (same columns as the original)
- Download `.gz` files listed in `path_https` (skips files that already exist)
- Ready for local runs or Slurm submission on HPC

---

## Requirements
- Python 3.9+
- Packages:
  - `pandas >= 1.3`
  - `requests >= 2.31`
  - `openpyxl >= 3.1` (only needed if you read Excel `.xlsx`)
- Sufficient disk space for large GWAS files
