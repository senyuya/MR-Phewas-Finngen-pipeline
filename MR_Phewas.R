# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(TwoSampleMR)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(export)
})

# ============================================================================
# 0. Parameters (exposure name and input files)
# ============================================================================
TARGET_NAME   <- "Target"                            # hide real target name
EXPOSURE_FILE <- "target_snps_phewas.txt"           # your exposure SNP file
MANIFEST_FILE <- "../../RARRES2_FGL1/finngen_R12_manifest.tsv"
OUTCOME_DIR   <- "extracted_snps_panos"              # folder with extracted Finngen R12 files

# ============================================================================
# 1. MR analysis helper
# ============================================================================
MR_function <- function(exposure.data, outcome.data) {
  # Harmonize
  final.data <- harmonise_data(exposure.data, outcome.data, action = 1)
  final.data <- final.data[!duplicated(final.data), ]
  
  # MR methods
  results <- mr(final.data, method_list = c("mr_ivw_fe", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  
  # Heterogeneity & pleiotropy
  heterogeneity <- mr_heterogeneity(final.data)
  pleiotropy    <- mr_pleiotropy_test(final.data)
  
  # Join diagnostics
  ivw_qpval <- heterogeneity$Q_pval[heterogeneity$method == "Inverse variance weighted"]
  results$IVW.Qpval             <- ivw_qpval
  results$pleiotropy.pval       <- pleiotropy$pval
  results$pleiotropy.intercept  <- pleiotropy$egger_intercept
  results$pleiotropy.se         <- pleiotropy$se
  
  results
}

# ============================================================================
# 2. Load exposure data
# ============================================================================
EXPOSURE <- fread(EXPOSURE_FILE, data.table = FALSE)

# ============================================================================
# 3. Load outcome data
# ============================================================================
manifest_r12 <- fread(MANIFEST_FILE, data.table = FALSE)

merge_files_correctly <- function(file_list, mapping_df, source_type) {
  combined_df <- data.frame()
  for (file_path in file_list) {
    dat <- read.csv(file_path, stringsAsFactors = FALSE, sep = "\t")
    file_name <- basename(file_path)
    file_base <- tools::file_path_sans_ext(file_name)
    phenocode <- sub("^(extracted?_finngen_R12_)(.*)$", "\\2", file_base)
    
    mapping_info <- mapping_df %>% dplyr::filter(phenocode == !!phenocode)
    if (nrow(mapping_info) == 0) {
      warning(paste("Phenocode", phenocode, "not found in manifest"))
      next
    }
    
    dat$phenocode <- mapping_info$phenocode
    dat$phenotype <- mapping_info$phenotype
    dat$category  <- mapping_info$category
    dat$source    <- source_type
    
    combined_df <- bind_rows(combined_df, dat)
  }
  combined_df
}

Outcome_files <- list.files(OUTCOME_DIR, pattern = "extracted_finngen_R12_.*\\.(csv|txt)$", full.names = TRUE)
OUTCOME_COMBINED <- merge_files_correctly(Outcome_files, manifest_r12, TARGET_NAME)

# ============================================================================
# 4. Format outcome and run MR
# ============================================================================
outcome_df <- format_data(
  OUTCOME_COMBINED,
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  pval_col = "pval",
  effect_allele_col = "alt",
  eaf_col = "af_alt",
  other_allele_col = "ref"
)

MR_RES <- MR_function(EXPOSURE, outcome_df)

# ============================================================================
# 5. Select appropriate IVW flavor and apply FDR
# ============================================================================
filtered_ivw <- MR_RES %>%
  group_by(exposure, outcome) %>%
  filter(if (any(method == "Inverse variance weighted" & IVW.Qpval < 0.05)) {
    method == "Inverse variance weighted"
  } else {
    method == "Inverse variance weighted (fixed effects)"
  }) %>%
  ungroup()

MR_TARGET_ivwfe <- filtered_ivw[, c("outcome", "exposure", "method", "nsnp", "b", "se", "pval")]
MR_TARGET_ivwfe$pval_adjust <- p.adjust(MR_TARGET_ivwfe$pval, method = "BH")

# Add outcome categories
match_idx <- match(MR_TARGET_ivwfe$outcome, OUTCOME_COMBINED$phenotype)
MR_TARGET_ivwfe$group_narrow     <- OUTCOME_COMBINED$category[match_idx]
MR_TARGET_ivwfe$neg_log10_qpvalue <- -log10(MR_TARGET_ivwfe$pval_adjust)

# ============================================================================
# 6. Mapping to simplified and broad disease groups
# ============================================================================
simplified_names <- c(
  "Alcohol related diseases" = "Alcohol-related",
  "Asthma and related endpoints" = "Asthma-related",
  "Cardiometabolic endpoints" = "Cardiometabolic",
  "Comorbidities of Asthma" = "Asthma comorbidities",
  "Comorbidities of COPD" = "COPD comorbidities",
  "Comorbidities of Diabetes" = "Diabetes comorbidities",
  "Comorbidities of Gastrointestinal endpoints" = "GI comorbidities",
  "Comorbidities of Interstitial lung disease endpoints" = "Interstitial lung comorbidities",
  "Comorbidities of Neurological endpoints" = "Neurological comorbidities",
  "COPD and related endpoints" = "COPD-related",
  "Diabetes endpoints" = "Diabetes",
  "Diseases marked as autimmune origin" = "Autoimmune diseases",
  "Drug purchase endpoints" = "Drug-related",
  "Gastrointestinal endpoints" = "Gastrointestinal",
  "I Certain infectious and parasitic diseases (AB1_)" = "Infectious diseases",
  "II Neoplasms from hospital discharges (CD2_)" = "Neoplasms (hospital)",
  "II Neoplasms, from cancer register (ICD-O-3)" = "Neoplasms (register)",
  "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)" = "Blood/immune",
  "Interstitial lung disease endpoints" = "Interstitial lung",
  "IV Endocrine, nutritional and metabolic diseases (E4_)" = "Endocrine/metabolic",
  "IX Diseases of the circulatory system (I9_)" = "Circulatory system",
  "Miscellaneous, not yet classified endpoints" = "Miscellaneous",
  "Neurological endpoints" = "Neurological",
  "Other, not yet classified endpoints (same as #MISC)" = "Other classified",
  "Psychiatric endpoints from Katri Räikkönen" = "Psychiatric",
  "Quantitative endpoints" = "Quantitative",
  "Rheuma endpoints" = "Rheuma",
  "V Mental and behavioural disorders (F5_)" = "Mental/behavioral",
  "VI Diseases of the nervous system (G6_)" = "Nervous system",
  "VII Diseases of the eye and adnexa (H7_)" = "Eye disorders",
  "VIII Diseases of the ear and mastoid process (H8_)" = "Ear disorders",
  "X Diseases of the respiratory system (J10_)" = "Respiratory system",
  "XI Diseases of the digestive system (K11_)" = "Digestive system",
  "XII Diseases of the skin and subcutaneous tissue (L12_)" = "Skin disorders",
  "XIII Diseases of the musculoskeletal system and connective tissue (M13_)" = "Musculoskeletal system",
  "XIV Diseases of the genitourinary system (N14_)" = "Genitourinary",
  "XIX Injury, poisoning and certain other consequences of external causes (ST19_)" = "Injuries/poisoning",
  "XV Pregnancy, childbirth and the puerperium (O15_)" = "Pregnancy-related",
  "XVI Certain conditions originating in the perinatal period (P16_)" = "Perinatal conditions",
  "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)" = "Congenital anomalies",
  "XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R18_)" = "Symptoms/signs",
  "XXI Factors influencing health status and contact with health services (Z21_)" = "Health factors",
  "XXII Codes for special purposes (U22_)" = "Special codes"
)

broad_groups <- c(
  "Asthma-related" = "Respiratory diseases",
  "COPD-related" = "Respiratory diseases",
  "COPD comorbidities" = "Respiratory diseases",
  "Asthma comorbidities" = "Respiratory diseases",
  "Respiratory system" = "Respiratory diseases",
  "Interstitial lung" = "Respiratory diseases",
  "Interstitial lung comorbidities" = "Respiratory diseases",
  
  "Diabetes" = "Endocrine/Metabolic-related diseases",
  "Diabetes comorbidities" = "Endocrine/Metabolic-related diseases",
  "Endocrine/metabolic" = "Endocrine/Metabolic-related diseases",
  
  "Mental/behavioral" = "Mental and neurological disorders",
  "Nervous system" = "Mental and neurological disorders",
  "Psychiatric" = "Mental and neurological disorders",
  "Neurological comorbidities" = "Mental and neurological disorders",
  "Neurological" = "Mental and neurological disorders",
  
  "Circulatory system" = "Cardiovascular diseases",
  "Cardiometabolic" = "Cardiovascular diseases",
  
  "Autoimmune diseases" = "Autoimmune and inflammatory diseases",
  "Rheuma" = "Autoimmune and inflammatory diseases",
  
  "Neoplasms (hospital)" = "Cancer and neoplasms",
  "Neoplasms (register)" = "Cancer and neoplasms",
  
  "Gastrointestinal" = "Digestive diseases",
  "GI comorbidities" = "Digestive diseases",
  "Digestive system" = "Digestive diseases",
  
  "Genitourinary" = "Genitourinary diseases",
  
  "Skin disorders" = "Dermatologic diseases",
  
  "Injuries/poisoning" = "Injuries, Symptoms & Signs, Musculoskeletal",
  "Symptoms/signs" = "Injuries, Symptoms & Signs, Musculoskeletal",
  "Musculoskeletal system" = "Injuries, Symptoms & Signs, Musculoskeletal",
  
  "Infectious diseases" = "Infectious diseases",
  "Blood/immune" = "Blood/immune diseases",
  
  "Alcohol-related" = "Miscellaneous",
  "Drug-related" = "Miscellaneous",
  "Health factors" = "Miscellaneous",
  "Eye disorders" = "Eye/ear disorders",
  "Ear disorders" = "Eye/ear disorders",
  "Pregnancy-related" = "Pregnancy/Perinatal-related diseases",
  "Perinatal conditions" = "Pregnancy/Perinatal-related diseases",
  "Congenital anomalies" = "Developmental disorders",
  "Special codes" = "Miscellaneous",
  "Quantitative" = "Miscellaneous",
  "Miscellaneous" = "Miscellaneous",
  "Other classified" = "Miscellaneous"
)

MR_TARGET_ivwfe <- MR_TARGET_ivwfe %>%
  mutate(
    group_narrow_simplified = simplified_names[group_narrow],
    broad_group = broad_groups[group_narrow_simplified]
  )

check_mapping <- function(data, data_name) {
  unmapped <- unique(data$group_narrow_simplified[is.na(data$broad_group)])
  if (length(unmapped) > 0) {
    cat(paste("Unmapped groups in", data_name, ":\n"))
    print(unmapped)
  } else {
    cat(paste("All groups mapped successfully in", data_name, "\n"))
  }
  
  counts <- data %>%
    group_by(broad_group) %>%
    summarise(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  cat(paste("Broad group counts for", data_name, ":\n"))
  print(counts)
  cat("\n")
}

check_mapping(MR_TARGET_ivwfe, TARGET_NAME)

# ============================================================================
# 7. Manhattan-style ordering helper
# ============================================================================
desired_order <- c(
  "Infectious diseases",
  "Blood/immune diseases",
  "Digestive diseases",
  "Mental and neurological disorders",
  "Endocrine/Metabolic-related diseases",
  "Eye/ear disorders",
  "Injuries, Symptoms & Signs, Musculoskeletal",
  "Respiratory diseases",
  "Cardiovascular diseases",
  "Cancer and neoplasms",
  "Dermatologic diseases",
  "Genitourinary diseases",
  "Pregnancy/Perinatal-related diseases",
  "Autoimmune and inflammatory diseases",
  "Developmental disorders",
  "Miscellaneous"
)

# Use MR results, merge manifest, and filter for binary endpoints with sufficient cases
manifest_core <- manifest_r12 %>%
  dplyr::select(phenocode, phenotype, category, num_cases, num_controls) %>%
  distinct()

mr_res2 <- MR_RES %>%
  left_join(manifest_core, by = c("outcome" = "phenotype")) %>%
  filter(!grepl("Quantitative", category)) %>%
  filter(!is.na(num_cases) & num_cases >= 1000)

cat("Eligible unique outcomes after filters: ",
    mr_res2 %>% distinct(outcome) %>% nrow(), "\n")

# (Table 1) IVW-FE main analysis + FDR within exposure
MR_TARGET_ivwfe <- mr_res2 %>%
  filter(method == "Inverse variance weighted (fixed effects)") %>%
  group_by(exposure) %>%
  mutate(FDR = p.adjust(pval, method = "fdr"),
         neg_log10_qpvalue = -log10(FDR)) %>%
  ungroup() %>%
  mutate(group_narrow = category,
         group_narrow_simplified = unname(simplified_names[group_narrow]),
         broad_group = unname(broad_groups[group_narrow_simplified]))

# (Table 2) Wide table with diagnostics
method_key <- c(
  "Inverse variance weighted (fixed effects)" = "ivwfe",
  "Inverse variance weighted"                 = "ivw",
  "MR Egger"                                  = "egger",
  "Weighted median"                           = "wm"
)

diag_cols <- mr_res2 %>%
  group_by(exposure, id.exposure, outcome, id.outcome) %>%
  summarise(
    nsnp               = first(nsnp),
    phenocode          = first(phenocode),
    category           = first(category),
    num_cases          = first(num_cases),
    num_controls       = first(num_controls),
    IVW_Qpval          = suppressWarnings(first(IVW.Qpval)),
    egger_intercept    = suppressWarnings(first(pleiotropy.intercept)),
    egger_intercept_se = suppressWarnings(first(pleiotropy.se)),
    egger_intercept_p  = suppressWarnings(first(pleiotropy.pval)),
    .groups = "drop"
  )

effects_wide <- mr_res2 %>%
  mutate(method_clean = recode(method, !!!method_key)) %>%
  select(exposure, id.exposure, outcome, id.outcome, method_clean, b, se, pval) %>%
  distinct() %>%
  tidyr::pivot_wider(
    id_cols   = c(exposure, id.exposure, outcome, id.outcome),
    names_from  = method_clean,
    values_from = c(b, se, pval),
    names_glue  = "{.value}_{method_clean}",
    values_fn   = list(b = ~ first(.), se = ~ first(.), pval = ~ first(.))
  )

MR_TARGET_wide <- diag_cols %>%
  left_join(effects_wide, by = c("exposure","id.exposure","outcome","id.outcome")) %>%
  group_by(exposure) %>%
  mutate(FDR_ivwfe = p.adjust(pval_ivwfe, method = "fdr"),
         neg_log10_qpvalue = -log10(FDR_ivwfe)) %>%
  ungroup() %>%
  mutate(group_narrow = category,
         group_narrow_simplified = unname(simplified_names[group_narrow]),
         broad_group = unname(broad_groups[group_narrow_simplified]))

# Sanity check: each (exposure, outcome) should be unique now
MR_TARGET_wide %>%
  count(exposure, outcome) %>%
  filter(n > 1) %>%
  print(n = 20)

# (Table 3) Strict robust-significance filter
#  a) IVW-FE: FDR < 0.05
#  b) WM same direction and p < 0.05
#  c) Egger intercept p > 0.05
#  d) If heterogeneity (Q p < 0.05): WM & Egger same direction as IVW-FE, and |b_egger - b_ivwfe|/|b_ivwfe| ≤ 0.20
MR_TARGET_signif <- MR_TARGET_wide %>%
  mutate(
    cond_a = !is.na(FDR_ivwfe) & FDR_ivwfe < 0.05,
    cond_b = !is.na(b_wm) & !is.na(b_ivwfe) &
      sign(b_wm) == sign(b_ivwfe) &
      !is.na(pval_wm) & pval_wm < 0.05,
    cond_c = is.na(egger_intercept_p) | egger_intercept_p > 0.05,
    hetero = !is.na(IVW_Qpval) & IVW_Qpval < 0.05,
    dir_consistent_egger = !is.na(b_egger) & !is.na(b_ivwfe) &
      sign(b_egger) == sign(b_ivwfe),
    egger_within20 = !is.na(b_egger) & !is.na(b_ivwfe) &
      abs(b_ivwfe) > 0 &
      abs(b_egger - b_ivwfe) / abs(b_ivwfe) <= 0.20,
    cond_d = (!hetero) | (hetero & dir_consistent_egger & cond_b & egger_within20),
    PASS_robust = cond_a & cond_b & cond_c & cond_d
  ) %>%
  filter(PASS_robust) %>%
  arrange(exposure, FDR_ivwfe) %>%
  mutate(group_narrow = category,
         group_narrow_simplified = unname(simplified_names[group_narrow]),
         broad_group = unname(broad_groups[group_narrow_simplified]))

cat("Robustly significant entries: ", nrow(MR_TARGET_signif), "\n")

# Optional: pre-processed data for plotting order
data_processed_target <- MR_TARGET_ivwfe %>%
  mutate(broad_group = factor(broad_group, levels = desired_order)) %>%
  arrange(broad_group, neg_log10_qpvalue) %>%
  mutate(Chromosome = dplyr::row_number())

# ============================================================================
# 8. Export
# ============================================================================
write.xlsx(
  list(
    "IVWFE_only"      = MR_TARGET_ivwfe,
    "AllMethods_wide" = MR_TARGET_wide,
    "Robust_signif"   = MR_TARGET_signif
  ),
  file = "MR_Target_FinngenR12_results.xlsx",
  overwrite = TRUE
)
