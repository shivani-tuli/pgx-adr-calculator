# PGx-ADR Risk Calculator - Core Algorithm
# ==========================================
# Quantitative pharmacogenomic ADR risk scoring with two modes:
#
# 1. PharmCAT Mode (recommended): VCF → PharmCAT star allele calling →
#    diplotype-based risk scoring using CPIC-aligned phenotype mappings
#
# 2. rsID Mode (fallback): Direct rsID → risk ratio lookup for manual
#    variant entry when PharmCAT is not available
#
# Key innovation: Evidence-weighted geometric mean model that outputs
# numerical risk ratios from either diplotype or variant-level data.

# ---- Evidence Level Weights ----
# Maps PharmGKB/CPIC evidence levels to numerical weights
EVIDENCE_WEIGHTS <- c(
  "1A" = 1.0,   # CPIC guideline, FDA label
  "1B" = 0.9,   # Strong evidence, replicated studies
  "2A" = 0.7,   # VIP gene, PharmGKB annotation
  "2B" = 0.5,   # PharmGKB annotation, moderate evidence
  "3"  = 0.3,   # Single significant study
  "4"  = 0.1    # Case report or in vitro
)

# Alternative weight schemes for sensitivity analysis
WEIGHT_SCHEMES <- list(
  "CPIC-aligned" = c("1A"=1.0, "1B"=0.9, "2A"=0.7, "2B"=0.5, "3"=0.3, "4"=0.1),
  "Binary"       = c("1A"=1.0, "1B"=1.0, "2A"=1.0, "2B"=0.5, "3"=0.0, "4"=0.0),
  "Linear"       = c("1A"=1.0, "1B"=0.83, "2A"=0.67, "2B"=0.5, "3"=0.33, "4"=0.17)
)

# ---- Gene-Specific Risk Models ----
# For genes where the multiplicative (risk^allele_count) model is wrong,
# use published heterozygote and homozygote risk ratios directly.
# Sources: CPIC guidelines + meta-analyses (PMIDs in comments)
GENE_RISK_MODELS <- list(
  # DPYD: het ~2.54x, hom ~25.6x (not 2.54^2=6.45x)
  # PMID:27296574 (het), PMID:Henricks 2019 JCO (hom)
  "DPYD"   = list(het = 2.54, hom = 25.6),
  # TPMT: het ~3.9x, hom ~6.67x (PMID:31342537, 25201288)
  "TPMT"   = list(het = 3.9,  hom = 6.67),
  # NUDT15: het ~3.0x, hom ~25.0x (PMID:26878728)
  "NUDT15" = list(het = 3.0,  hom = 25.0),
  # UGT1A1: het ~1.9x, hom ~5.34x (PMID:23529007, 28529590)
  "UGT1A1" = list(het = 1.9,  hom = 5.34),
  # G6PD: hemizygous males are fully affected; risk stored is for deficient state
  "G6PD"   = list(het = 1.5,  hom = 8.0)
)

# ---- Gene-Specific Variant Risk ----
# Replaces the naive risk^allele_count model
get_variant_risk <- function(gene, base_rr, allele_count, genetic_model = "multiplicative") {
  if (allele_count == 0) return(1.0)

  if (genetic_model == "activity_score" && gene %in% names(GENE_RISK_MODELS)) {
    model <- GENE_RISK_MODELS[[gene]]
    if (allele_count == 1) return(model$het)
    if (allele_count >= 2) return(model$hom)
    return(1.0)
  }

  if (genetic_model == "x_linked" && gene %in% names(GENE_RISK_MODELS)) {
    model <- GENE_RISK_MODELS[[gene]]
    # For X-linked: even het can be fully affected (hemizygous males)
    # Use the base risk ratio from data (which represents the deficient state)
    if (allele_count >= 1) return(base_rr)
    return(1.0)
  }

  # Default: multiplicative (codominant) model
  return(base_rr ^ allele_count)
}

# ---- Combined CI via Delta Method ----
# Propagates uncertainty from individual variant CIs to the combined risk score
compute_combined_ci <- function(log_rrs, ci_lowers, ci_uppers, weights) {
  # Estimate SE on log scale from CIs
  log_ses <- mapply(function(rr, lo, hi) {
    if (is.na(lo) || is.na(hi) || lo <= 0 || hi <= 0) {
      return(0.5)  # conservative default SE when CI unavailable
    }
    (log(hi) - log(lo)) / (2 * 1.96)
  }, log_rrs, ci_lowers, ci_uppers)

  w <- weights / sum(weights)
  log_combined <- sum(log_rrs * w)

  # Variance of weighted mean (assumes independent sources)
  log_var <- sum((w^2) * (log_ses^2))
  log_se  <- sqrt(log_var)

  has_any_ci <- !all(is.na(ci_lowers) | is.na(ci_uppers))

  list(
    combined_rr = exp(log_combined),
    ci_lower    = exp(log_combined - 1.96 * log_se),
    ci_upper    = exp(log_combined + 1.96 * log_se),
    has_ci      = has_any_ci
  )
}

# ---- Sensitivity Analysis ----
# Runs risk scoring under 3 different weight schemes to check robustness
run_sensitivity_analysis <- function(matched_data) {
  if (is.null(matched_data) || nrow(matched_data) == 0) return(NULL)

  results <- lapply(names(WEIGHT_SCHEMES), function(scheme_name) {
    w <- WEIGHT_SCHEMES[[scheme_name]]
    matched_data %>%
      dplyr::mutate(sw = w[evidence_level]) %>%
      dplyr::group_by(drug_name) %>%
      dplyr::summarise(
        rr = exp(sum(log(variant_risk) * sw) / sum(sw)),
        scheme = scheme_name,
        .groups = "drop"
      )
  })

  combined <- dplyr::bind_rows(results)

  # For each drug, check if risk category is the same across all schemes
  stability <- combined %>%
    dplyr::mutate(category = categorize_risk(rr)) %>%
    dplyr::group_by(drug_name) %>%
    dplyr::summarise(
      n_categories = dplyr::n_distinct(category),
      min_rr = min(rr),
      max_rr = max(rr),
      sensitivity_badge = ifelse(dplyr::n_distinct(category) == 1,
                                  "Robust", "Weight-Sensitive"),
      .groups = "drop"
    )

  return(stability)
}

# ---- Risk Categories ----
categorize_risk <- function(risk_ratio) {
  dplyr::case_when(
    risk_ratio >= 5.0 ~ "Very High",
    risk_ratio >= 2.5 ~ "High",
    risk_ratio >= 1.5 ~ "Moderate",
    risk_ratio >= 1.2 ~ "Slightly Elevated",
    TRUE              ~ "Normal"
  )
}

risk_color <- function(category) {
  dplyr::case_when(
    category == "Very High"         ~ "#dc2626",
    category == "High"              ~ "#ea580c",
    category == "Moderate"          ~ "#ca8a04",
    category == "Slightly Elevated" ~ "#65a30d",
    category == "Normal"            ~ "#16a34a",
    TRUE                            ~ "#6b7280"
  )
}

# ---- VCF Parser ----
parse_vcf <- function(vcf_path) {
  lines <- readLines(vcf_path)
  # Skip header lines
  data_lines <- lines[!grepl("^##", lines)]
  if (length(data_lines) < 2) return(NULL)

  # Parse header
  header <- strsplit(data_lines[1], "\t")[[1]]
  header[1] <- sub("^#", "", header[1])

  # Parse data
  data_rows <- data_lines[-1]
  if (length(data_rows) == 0) return(NULL)

  df <- do.call(rbind, lapply(strsplit(data_rows, "\t"), function(x) {
    length(x) <- length(header)
    x
  }))
  colnames(df) <- header
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # Extract rsIDs and genotypes
  if (!"ID" %in% colnames(df)) return(NULL)

  # Find sample column (typically column 10+)
  sample_cols <- setdiff(colnames(df),
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))

  if (length(sample_cols) == 0) return(NULL)

  sample_col <- sample_cols[1]
  format_col <- df$FORMAT

  # Extract GT field
  genotypes <- sapply(seq_len(nrow(df)), function(i) {
    fmt <- strsplit(format_col[i], ":")[[1]]
    vals <- strsplit(df[[sample_col]][i], ":")[[1]]
    gt_idx <- which(fmt == "GT")
    if (length(gt_idx) > 0 && gt_idx <= length(vals)) vals[gt_idx] else NA
  })

  result <- data.frame(
    rsid     = df$ID,
    genotype = genotypes,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(genotype), rsid != ".")

  # Count alt alleles
result$allele_count <- sapply(result$genotype, function(gt) {
    alleles <- unlist(strsplit(gt, "[/|]"))
    sum(alleles != "0" & alleles != ".", na.rm = TRUE)
  })

  return(result)
}

# ---- Manual Entry Parser ----
parse_manual_entry <- function(entry_text) {
  # Expected format: "rs1799853:0/1, rs4149056:1/1"
  entries <- trimws(unlist(strsplit(entry_text, "[,;\n]")))
  entries <- entries[nchar(entries) > 0]

  if (length(entries) == 0) return(NULL)

  result <- do.call(rbind, lapply(entries, function(e) {
    parts <- trimws(unlist(strsplit(e, ":")))
    if (length(parts) != 2) return(NULL)
    rsid <- parts[1]
    gt   <- parts[2]
    alleles <- unlist(strsplit(gt, "[/|]"))
    allele_count <- sum(alleles != "0" & alleles != ".", na.rm = TRUE)
    data.frame(rsid = rsid, genotype = gt, allele_count = allele_count,
               stringsAsFactors = FALSE)
  }))

  return(result)
}

# ---- Core Risk Scoring Algorithm ----
#
# For each drug, the combined risk ratio is:
#   Combined_Risk = exp( Σ [log(variant_risk_i) × w_i] / Σ w_i )
#
# Key improvements over naive model:
#   1. Gene-specific risk: uses activity_score/x_linked models where appropriate
#   2. LD-aware: deduplicates variants in same haplotype group per drug
#   3. CI propagation: delta method on log-scale geometric mean
#   4. Sensitivity analysis: runs 3 weight schemes to check robustness

calculate_risk_scores <- function(patient_variants, pgx_associations) {
  if (is.null(patient_variants) || nrow(patient_variants) == 0) {
    return(NULL)
  }

  # Match patient variants to known PGx associations
  matched <- dplyr::inner_join(
    patient_variants,
    pgx_associations,
    by = c("rsid" = "variant_rsid")
  )

  if (nrow(matched) == 0) return(NULL)

  # Handle missing new columns gracefully (backward compat)
  if (!"genetic_model" %in% names(matched)) {
    matched$genetic_model <- "multiplicative"
  }
  if (!"haplotype_group" %in% names(matched)) {
    matched$haplotype_group <- paste0(matched$gene, "_", matched$allele_name)
  }

  # ---- LD-Aware Deduplication ----
  # Within each drug + haplotype_group, keep only the highest-risk variant
  # to prevent double-counting rsIDs that tag the same star allele
  matched <- matched %>%
    dplyr::group_by(drug_name, haplotype_group) %>%
    dplyr::slice_max(order_by = risk_ratio * allele_count, n = 1,
                      with_ties = FALSE) %>%
    dplyr::ungroup()

  # ---- Gene-Specific Risk Calculation ----
  matched <- matched %>%
    dplyr::mutate(
      evidence_weight = EVIDENCE_WEIGHTS[evidence_level],
      variant_risk    = mapply(get_variant_risk, gene, risk_ratio,
                                allele_count, genetic_model),
      weighted_log_risk = log(pmax(variant_risk, 1e-6)) * evidence_weight
    )

  # ---- Sensitivity Analysis ----
  sensitivity <- run_sensitivity_analysis(matched)

  # ---- Aggregate by Drug with CI Propagation ----
  drug_scores <- matched %>%
    dplyr::group_by(drug_name, drug_class) %>%
    dplyr::summarise(
      combined_risk_ratio = exp(sum(weighted_log_risk) / sum(evidence_weight)),
      ci_lower = compute_combined_ci(
        log(pmax(variant_risk, 1e-6)),
        suppressWarnings(as.numeric(ci_lower)),
        suppressWarnings(as.numeric(ci_upper)),
        evidence_weight
      )$ci_lower,
      ci_upper = compute_combined_ci(
        log(pmax(variant_risk, 1e-6)),
        suppressWarnings(as.numeric(ci_lower)),
        suppressWarnings(as.numeric(ci_upper)),
        evidence_weight
      )$ci_upper,
      has_ci = compute_combined_ci(
        log(pmax(variant_risk, 1e-6)),
        suppressWarnings(as.numeric(ci_lower)),
        suppressWarnings(as.numeric(ci_upper)),
        evidence_weight
      )$has_ci,
      n_variants          = dplyr::n(),
      genes_involved      = paste(unique(gene), collapse = ", "),
      models_used         = paste(unique(genetic_model), collapse = ", "),
      primary_adr         = dplyr::first(primary_adr),
      adr_severity        = dplyr::first(adr_severity),
      max_evidence        = dplyr::first(sort(unique(evidence_level))),
      variant_details     = paste(
        sprintf("%s (%s %s, RR=%.1fx)",
                rsid, gene, genotype, variant_risk),
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      risk_category = categorize_risk(combined_risk_ratio),
      risk_color    = risk_color(risk_category),
      ci_display = ifelse(has_ci,
        sprintf("%.2f (%.2f–%.2f)", combined_risk_ratio, ci_lower, ci_upper),
        sprintf("%.2f (CI unavailable)", combined_risk_ratio)
      )
    ) %>%
    dplyr::arrange(dplyr::desc(combined_risk_ratio))

  # Attach sensitivity badges
  if (!is.null(sensitivity)) {
    drug_scores <- dplyr::left_join(drug_scores, sensitivity,
                                     by = "drug_name")
  }

  return(drug_scores)
}

# ---- Gene-Level Summary ----
summarize_by_gene <- function(patient_variants, pgx_associations) {
  if (is.null(patient_variants) || nrow(patient_variants) == 0) return(NULL)

  matched <- dplyr::inner_join(
    patient_variants,
    pgx_associations,
    by = c("rsid" = "variant_rsid")
  )

  if (nrow(matched) == 0) return(NULL)

  # Handle missing columns
  if (!"genetic_model" %in% names(matched)) {
    matched$genetic_model <- "multiplicative"
  }

  gene_summary <- matched %>%
    dplyr::mutate(
      computed_risk = mapply(get_variant_risk, gene, risk_ratio,
                             allele_count, genetic_model)
    ) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_variants    = dplyr::n_distinct(rsid),
      drugs_affected = paste(unique(drug_name), collapse = ", "),
      n_drugs       = dplyr::n_distinct(drug_name),
      variants      = paste(unique(rsid), collapse = ", "),
      max_risk      = max(computed_risk, na.rm = TRUE),
      .groups       = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(max_risk))

  return(gene_summary)
}

# ---- Population Frequency Analysis ----
# Compares patient variant frequencies to population baselines
population_context <- function(patient_variants, pgx_associations) {
  if (is.null(patient_variants) || nrow(patient_variants) == 0) return(NULL)

  matched <- dplyr::inner_join(
    patient_variants,
    pgx_associations,
    by = c("rsid" = "variant_rsid")
  )

  if (nrow(matched) == 0) return(NULL)

  pop_summary <- matched %>%
    dplyr::select(rsid, gene, allele_count, genotype,
                  global_maf, eur_af, afr_af, eas_af, sas_af, amr_af) %>%
    dplyr::distinct()

  return(pop_summary)
}

# ============================================================================
# Phenoconversion Module
# ============================================================================
# Models how co-administered CYP inhibitors/inducers shift genetic metabolizer
# phenotypes, adjusting ADR risk scores accordingly.
#
# Example: CYP2D6 Normal Metabolizer + fluoxetine (strong CYP2D6 inhibitor)
#          → effectively Poor Metabolizer → codeine RR jumps from 1.0x to 2.5x

# ---- Phenotype Ladder ----
# Ordered from highest to lowest activity; index position used for shifting
PHENOTYPE_LADDER <- c(
  "Ultrarapid Metabolizer",
  "Rapid Metabolizer",
  "Normal Metabolizer",
  "Intermediate Metabolizer",
  "Poor Metabolizer"
)

# ---- Apply Phenoconversion ----
# Takes genetic phenotypes (from PharmCAT or manual) and co-medications,
# returns adjusted phenotypes per enzyme.
apply_phenoconversion <- function(genetic_phenotypes, concomitant_meds,
                                  inhibitor_db) {
  # genetic_phenotypes: data.frame with columns: gene, phenotype
  # concomitant_meds: character vector of drug names
  # inhibitor_db: data.frame (cyp_inhibitors_inducers.csv)

  if (length(concomitant_meds) == 0 || is.null(concomitant_meds)) {
    return(genetic_phenotypes %>%
      dplyr::mutate(
        adjusted_phenotype = phenotype,
        phenoconversion_active = FALSE,
        interacting_drugs = "",
        shift_reason = ""
      ))
  }

  # Find all interactions for the patient's concomitant meds
  active_interactions <- inhibitor_db %>%
    dplyr::filter(tolower(drug_name) %in% tolower(concomitant_meds))

  if (nrow(active_interactions) == 0) {
    return(genetic_phenotypes %>%
      dplyr::mutate(
        adjusted_phenotype = phenotype,
        phenoconversion_active = FALSE,
        interacting_drugs = "",
        shift_reason = ""
      ))
  }

  # For each gene, compute net phenotype shift
  results <- genetic_phenotypes

  results$adjusted_phenotype <- results$phenotype
  results$phenoconversion_active <- FALSE
  results$interacting_drugs <- ""
  results$shift_reason <- ""

  for (i in seq_len(nrow(results))) {
    gene <- results$gene[i]
    genetic_pheno <- results$phenotype[i]

    # Find interactions targeting this enzyme
    gene_interactions <- active_interactions %>%
      dplyr::filter(enzyme == gene)

    if (nrow(gene_interactions) == 0) next

    # Use the strongest effect (most negative shift for inhibitors)
    # If both inhibitor and inducer present, inhibitor wins (conservative)
    inhibitors <- gene_interactions %>% dplyr::filter(effect == "inhibitor")
    inducers <- gene_interactions %>% dplyr::filter(effect == "inducer")

    net_shift <- 0
    reason_parts <- c()

    if (nrow(inhibitors) > 0) {
      # Take strongest inhibitor (most negative shift)
      strongest <- inhibitors %>%
        dplyr::arrange(phenotype_shift) %>%
        dplyr::slice(1)
      net_shift <- as.numeric(strongest$phenotype_shift)
      reason_parts <- c(reason_parts,
        sprintf("%s (%s %s inhibitor)",
                strongest$drug_name, strongest$potency, gene))
    } else if (nrow(inducers) > 0) {
      # Take strongest inducer (most positive shift)
      strongest <- inducers %>%
        dplyr::arrange(dplyr::desc(phenotype_shift)) %>%
        dplyr::slice(1)
      net_shift <- as.numeric(strongest$phenotype_shift)
      reason_parts <- c(reason_parts,
        sprintf("%s (%s %s inducer)",
                strongest$drug_name, strongest$potency, gene))
    }

    if (net_shift != 0) {
      # Find current position on ladder
      current_pos <- match(genetic_pheno, PHENOTYPE_LADDER)
      if (is.na(current_pos)) {
        # Try partial matching for non-standard phenotype names
        current_pos <- 3  # Default to NM if unknown
      }

      # Shift position (negative shift = toward PM = higher index)
      new_pos <- current_pos - net_shift
      new_pos <- max(1, min(length(PHENOTYPE_LADDER), new_pos))

      new_pheno <- PHENOTYPE_LADDER[new_pos]

      if (new_pheno != genetic_pheno) {
        results$adjusted_phenotype[i] <- new_pheno
        results$phenoconversion_active[i] <- TRUE
        results$interacting_drugs[i] <- paste(
          gene_interactions$drug_name, collapse = ", ")
        results$shift_reason[i] <- paste(reason_parts, collapse = "; ")
      }
    }
  }

  results
}

# ---- Calculate Adjusted Risk ----
# Re-scores risk using phenoconverted phenotypes instead of genetic phenotypes
calculate_adjusted_risk <- function(phenoconverted, diplotype_data) {
  # phenoconverted: output of apply_phenoconversion()
  # diplotype_data: diplotype risk associations database

  if (is.null(phenoconverted) || nrow(phenoconverted) == 0) return(NULL)
  if (!any(phenoconverted$phenoconversion_active)) return(NULL)

  # Create lookup using adjusted phenotypes
  adjusted_results <- phenoconverted %>%
    dplyr::mutate(
      lookup_key = paste0(gene, ":", adjusted_phenotype)
    )

  # Re-run diplotype risk calculation with adjusted phenotypes
  calculate_diplotype_risk(adjusted_results, diplotype_data)
}

# ---- Summarize Interactions ----
# Returns a summary table for the Drug Interactions tab
summarize_interactions <- function(concomitant_meds, inhibitor_db,
                                    genetic_phenotypes = NULL) {
  if (length(concomitant_meds) == 0) return(NULL)

  interactions <- inhibitor_db %>%
    dplyr::filter(tolower(drug_name) %in% tolower(concomitant_meds))

  if (nrow(interactions) == 0) return(NULL)

  summary_df <- interactions %>%
    dplyr::group_by(enzyme) %>%
    dplyr::summarize(
      interacting_drugs = paste(unique(drug_name), collapse = ", "),
      strongest_effect = dplyr::first(effect[abs(phenotype_shift) ==
                                               max(abs(phenotype_shift))]),
      strongest_potency = dplyr::first(potency[abs(phenotype_shift) ==
                                                 max(abs(phenotype_shift))]),
      max_shift = max(abs(phenotype_shift)),
      .groups = "drop"
    )

  # Add genetic and adjusted phenotypes if available
  if (!is.null(genetic_phenotypes)) {
    pheno_map <- stats::setNames(genetic_phenotypes$phenotype,
                                  genetic_phenotypes$gene)
    adj_map <- if ("adjusted_phenotype" %in% names(genetic_phenotypes)) {
      stats::setNames(genetic_phenotypes$adjusted_phenotype,
                       genetic_phenotypes$gene)
    } else { pheno_map }

    summary_df$genetic_phenotype <- pheno_map[summary_df$enzyme]
    summary_df$adjusted_phenotype <- adj_map[summary_df$enzyme]
  }

  summary_df
}


# ============================================================================
# PharmCAT Integration Layer
# ============================================================================
# These functions bridge PharmCAT's star allele calling with the risk scoring
# engine. When a VCF is uploaded and PharmCAT is available, the pipeline is:
#   VCF → PharmCAT wrapper (Python) → star allele TSV → diplotype risk scoring

# ---- PharmCAT Paths ----
PHARMCAT_DIR <- file.path(getwd(), "pharmcat")
PHARMCAT_WRAPPER <- file.path(getwd(), "scripts", "pharmcat_wrapper.py")
PHARMCAT_PYTHON <- file.path(PHARMCAT_DIR, ".venv", "bin", "python3")
# Detect Java 17 dynamically (env var → platform-specific paths)
JAVA_HOME_17 <- Sys.getenv("JAVA_HOME", unset = "")
if (JAVA_HOME_17 == "" || !dir.exists(JAVA_HOME_17)) {
  java_candidates <- c(
    "/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home",  # macOS ARM
    "/usr/local/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home",     # macOS Intel
    "/usr/lib/jvm/java-17-openjdk-amd64",                              # Debian/Ubuntu
    "/usr/lib/jvm/java-17-openjdk",                                    # RHEL/Fedora
    Sys.getenv("JAVA_HOME", unset = "")
  )
  for (candidate in java_candidates) {
    if (candidate != "" && dir.exists(candidate)) {
      JAVA_HOME_17 <- candidate
      break
    }
  }
}
PHARMGKB_UPDATER <- file.path(getwd(), "scripts", "update_pharmgkb_data.py")
PHARMGKB_REPORT <- file.path(getwd(), "data", "pharmgkb_update_report.json")

# ---- PharmGKB Data Freshness ----
check_data_freshness <- function(csv_path = NULL) {
  # Check the update report file for last update timestamp
  if (file.exists(PHARMGKB_REPORT)) {
    report <- jsonlite::fromJSON(PHARMGKB_REPORT)
    last_updated <- report$timestamp
    if (!is.null(last_updated)) {
      update_time <- as.POSIXct(last_updated, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
      age_days <- as.numeric(difftime(Sys.time(), update_time, units = "days"))
      return(list(
        fresh = age_days < 30,
        last_updated = last_updated,
        age_days = round(age_days, 1),
        source = "PharmGKB_API"
      ))
    }
  }

  # Fallback: check CSV last_updated column
  if (!is.null(csv_path) && file.exists(csv_path)) {
    df <- read.csv(csv_path, nrows = 1, stringsAsFactors = FALSE)
    if ("last_updated" %in% names(df) && nchar(df$last_updated[1]) > 0) {
      return(list(
        fresh = TRUE,
        last_updated = df$last_updated[1],
        age_days = 0,
        source = df$source[1]
      ))
    }
  }

  return(list(fresh = FALSE, last_updated = "Never", age_days = Inf, source = "Manual"))
}

# ---- Trigger PharmGKB Update ----
trigger_data_update <- function() {
  if (!file.exists(PHARMGKB_UPDATER)) {
    return(list(success = FALSE, message = "Updater script not found"))
  }

  result <- system2("python3", args = shQuote(PHARMGKB_UPDATER),
                     stdout = TRUE, stderr = TRUE)
  exit_code <- attr(result, "status")

  if (is.null(exit_code) || exit_code == 0) {
    return(list(success = TRUE, message = paste(result, collapse = "\n")))
  } else {
    return(list(success = FALSE, message = paste(result, collapse = "\n")))
  }
}

# ---- Check PharmCAT Availability ----
check_pharmcat_available <- function() {
  jar_exists <- file.exists(file.path(PHARMCAT_DIR, "pharmcat.jar"))
  wrapper_exists <- file.exists(PHARMCAT_WRAPPER)
  java_exists <- file.exists(file.path(JAVA_HOME_17, "bin", "java"))

  list(
    available = jar_exists && wrapper_exists && java_exists,
    jar = jar_exists,
    wrapper = wrapper_exists,
    java = java_exists,
    version = if (jar_exists && java_exists) "3.2.0" else "not installed"
  )
}

# ---- Run PharmCAT Pipeline ----
run_pharmcat <- function(vcf_path) {
  if (!file.exists(vcf_path)) {
    stop("VCF file not found: ", vcf_path)
  }

  status <- check_pharmcat_available()
  if (!status$available) {
    stop("PharmCAT is not available. Missing: ",
         paste(c(
           if (!status$jar) "pharmcat.jar",
           if (!status$wrapper) "wrapper script",
           if (!status$java) "Java 17"
         ), collapse = ", "))
  }

  # Create temporary output file
  output_tsv <- tempfile(fileext = ".tsv")

  # Determine which Python to use
  python_cmd <- if (file.exists(PHARMCAT_PYTHON)) {
    PHARMCAT_PYTHON
  } else {
    "python3"
  }

  # Set up environment with Java 17
  env_vars <- paste0(
    "JAVA_HOME=", JAVA_HOME_17,
    " PATH=", JAVA_HOME_17, "/bin:",
    file.path(PHARMCAT_DIR, ".venv", "bin"), ":$PATH"
  )

  # Run the wrapper
  cmd <- paste(
    env_vars,
    python_cmd,
    shQuote(PHARMCAT_WRAPPER),
    shQuote(vcf_path),
    shQuote(output_tsv),
    "--pharmcat-dir", shQuote(PHARMCAT_DIR)
  )

  result <- system(cmd, intern = FALSE, ignore.stderr = FALSE)

  if (result != 0) {
    stop("PharmCAT pipeline failed with exit code: ", result)
  }

  if (!file.exists(output_tsv)) {
    stop("PharmCAT output file not created")
  }

  return(output_tsv)
}

# ---- Parse PharmCAT Results ----
parse_pharmcat_results <- function(tsv_path) {
  if (!file.exists(tsv_path)) return(NULL)

  results <- read.delim(tsv_path, stringsAsFactors = FALSE)

  if (nrow(results) == 0) return(NULL)

  # Clean up
  results$activity_score <- ifelse(
    results$activity_score %in% c("N/A", "", "None"),
    NA_real_,
    as.numeric(results$activity_score)
  )

  return(results)
}

# ---- Diplotype-Based Risk Scoring ----
# Maps PharmCAT diplotype calls to the diplotype_risk_associations database
# using gene + phenotype matching (more robust than exact diplotype string match)
calculate_diplotype_risk <- function(pharmcat_results, diplotype_db) {
  if (is.null(pharmcat_results) || nrow(pharmcat_results) == 0) return(NULL)
  if (is.null(diplotype_db) || nrow(diplotype_db) == 0) return(NULL)

  # Strategy: match by gene + phenotype (primary) or gene + diplotype (secondary)
  matched <- data.frame()

  for (i in seq_len(nrow(pharmcat_results))) {
    gene <- pharmcat_results$gene[i]
    phenotype <- pharmcat_results$phenotype[i]
    diplotype <- pharmcat_results$diplotype[i]

    # Try phenotype-based match first (most clinically relevant)
    pheno_match <- diplotype_db[
      diplotype_db$gene == gene &
      tolower(diplotype_db$phenotype) == tolower(phenotype),
    ]

    # Fall back to exact diplotype match
    if (nrow(pheno_match) == 0) {
      pheno_match <- diplotype_db[
        diplotype_db$gene == gene &
        diplotype_db$diplotype == diplotype,
      ]
    }

    if (nrow(pheno_match) > 0) {
      pheno_match$called_diplotype <- diplotype
      pheno_match$called_phenotype <- phenotype
      pheno_match$activity_score_called <- pharmcat_results$activity_score[i]
      pheno_match$source <- "PharmCAT"
      matched <- rbind(matched, pheno_match)
    }
  }

  if (nrow(matched) == 0) return(NULL)

  # Apply evidence weighting
  matched$evidence_weight <- EVIDENCE_WEIGHTS[matched$evidence_level]

  # Score per drug (using evidence-weighted model)
  drug_scores <- matched %>%
    dplyr::group_by(drug_name, drug_class) %>%
    dplyr::summarise(
      combined_risk_ratio = exp(
        sum(log(risk_ratio) * evidence_weight) / sum(evidence_weight)
      ),
      n_genes = dplyr::n_distinct(gene),
      genes_involved = paste(unique(gene), collapse = ", "),
      diplotypes = paste(
        sprintf("%s %s (%s)", gene, called_diplotype, called_phenotype),
        collapse = "; "
      ),
      primary_adr = dplyr::first(primary_adr),
      adr_severity = dplyr::first(adr_severity),
      max_evidence = dplyr::first(sort(unique(evidence_level))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      risk_category = categorize_risk(combined_risk_ratio),
      risk_color = risk_color(risk_category),
      scoring_mode = "PharmCAT (diplotype-based)"
    ) %>%
    dplyr::arrange(dplyr::desc(combined_risk_ratio))

  return(drug_scores)
}

# ---- Diplotype Summary ----
# Produces a gene-level summary of PharmCAT calls for the Diplotype Report tab
summarize_diplotypes <- function(pharmcat_results, diplotype_db) {
  if (is.null(pharmcat_results) || nrow(pharmcat_results) == 0) return(NULL)

  summary_df <- pharmcat_results %>%
    dplyr::select(gene, diplotype, phenotype, activity_score) %>%
    dplyr::distinct()

  # Add drug counts from the database
  summary_df$n_drugs <- sapply(seq_len(nrow(summary_df)), function(i) {
    g <- summary_df$gene[i]
    p <- summary_df$phenotype[i]
    nrow(diplotype_db[
      diplotype_db$gene == g &
      tolower(diplotype_db$phenotype) == tolower(p),
    ])
  })

  # Add max risk
  summary_df$max_risk <- sapply(seq_len(nrow(summary_df)), function(i) {
    g <- summary_df$gene[i]
    p <- summary_df$phenotype[i]
    matches <- diplotype_db[
      diplotype_db$gene == g &
      tolower(diplotype_db$phenotype) == tolower(p),
    ]
    if (nrow(matches) > 0) max(matches$risk_ratio) else NA_real_
  })

  return(summary_df)
}


# Genes on the X chromosome requiring sex-aware interpretation
X_LINKED_GENES <- c("G6PD")

# Drugs where risk is strongly dose-dependent
DOSE_DEPENDENT_DRUGS <- list(
  "Simvastatin"  = "SLCO1B1 myopathy risk is strongly dose-dependent (highest at 80mg). CPIC recommends avoiding simvastatin 80mg in *5 carriers.",
  "Atorvastatin"  = "Myopathy risk is dose-dependent. Consider lower doses in SLCO1B1 *5 carriers.",
  "Rosuvastatin"  = "Myopathy risk is dose-dependent. Consider lower doses in SLCO1B1 *5 carriers.",
  "Warfarin"      = "Risk depends on dose adjustment. Genotype-guided dosing significantly reduces bleeding events.",
  "Acenocoumarol" = "Risk depends on dose adjustment. Genotype-guided dosing significantly reduces bleeding events."
)

# Genes with known structural variation not captured by rsIDs
STRUCTURAL_VARIATION_GENES <- c("CYP2D6")

# Known CYP inhibitors that cause phenoconversion
PHENOCONVERSION_DRUGS <- list(
  "CYP2D6" = c("fluoxetine", "paroxetine", "bupropion", "quinidine", "terbinafine"),
  "CYP2C19" = c("omeprazole", "fluvoxamine", "fluconazole"),
  "CYP3A5" = c("ketoconazole", "itraconazole", "clarithromycin")
)

generate_clinical_warnings <- function(patient_variants, pgx_associations, risk_scores) {
  warnings <- list()
  
  if (is.null(patient_variants) || nrow(patient_variants) == 0) return(warnings)
  
  matched <- dplyr::inner_join(
    patient_variants, pgx_associations,
    by = c("rsid" = "variant_rsid")
  )
  if (nrow(matched) == 0) return(warnings)
  
  # 1. G6PD X-linked warning
  x_linked_hits <- matched %>%
    dplyr::filter(gene %in% X_LINKED_GENES, allele_count > 0)
  if (nrow(x_linked_hits) > 0) {
    drugs_affected <- paste(unique(x_linked_hits$drug_name), collapse = ", ")
    warnings <- c(warnings, list(list(
      type = "critical",
      icon = "exclamation-circle",
      title = "X-Linked Gene Detected: G6PD",
      message = paste0(
        "G6PD is X-linked. Hemizygous males (one variant copy) are FULLY affected, ",
        "equivalent to homozygous females. The displayed risk may UNDERESTIMATE ",
        "the true risk for male patients. Affected drugs: ", drugs_affected, ". ",
        "Verify patient sex before interpreting these results."
      )
    )))
  }
  
  # 2. CYP2D6 structural variation warning
  cyp2d6_hits <- matched %>%
    dplyr::filter(gene %in% STRUCTURAL_VARIATION_GENES)
  if (nrow(cyp2d6_hits) > 0) {
    warnings <- c(warnings, list(list(
      type = "warning",
      icon = "dna",
      title = "CYP2D6 Structural Variation Limitation",
      message = paste0(
        "CYP2D6 has complex structural variation (gene deletions, duplications, ",
        "hybrid alleles) that CANNOT be detected from rsID-based genotyping. ",
        "This tool may miss CYP2D6 gene duplications (ultra-rapid metabolizers) or ",
        "whole-gene deletions. Consider clinical CYP2D6 star-allele testing for ",
        "comprehensive assessment."
      )
    )))
  }
  
  # 3. Compound heterozygote warnings (multiple variants in same gene for same drug)
  multi_variant <- matched %>%
    dplyr::filter(allele_count > 0) %>%
    dplyr::group_by(gene, drug_name) %>%
    dplyr::summarise(
      n_variants = dplyr::n_distinct(rsid),
      variants = paste(unique(rsid), collapse = " + "),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_variants > 1)
  
  if (nrow(multi_variant) > 0) {
    for (i in seq_len(nrow(multi_variant))) {
      warnings <- c(warnings, list(list(
        type = "info",
        icon = "layer-group",
        title = paste0("Compound Variants in ", multi_variant$gene[i],
                       " for ", multi_variant$drug_name[i]),
        message = paste0(
          "Multiple variants detected in ", multi_variant$gene[i], " (",
          multi_variant$variants[i], ") affecting ", multi_variant$drug_name[i],
          ". The combined risk uses a multiplicative model, but compound ",
          "heterozygotes may have synergistic (greater-than-multiplicative) effects. ",
          "Clinical star-allele testing is recommended."
        )
      )))
    }
  }
  
  # 4. Dose-dependent drug warnings
  if (!is.null(risk_scores)) {
    dose_hits <- risk_scores %>%
      dplyr::filter(drug_name %in% names(DOSE_DEPENDENT_DRUGS))
    for (i in seq_len(nrow(dose_hits))) {
      drug <- dose_hits$drug_name[i]
      warnings <- c(warnings, list(list(
        type = "info",
        icon = "pills",
        title = paste0("Dose-Dependent Risk: ", drug),
        message = DOSE_DEPENDENT_DRUGS[[drug]]
      )))
    }
  }
  
  # 5. Phenoconversion: direct users to Drug Interactions module
  affected_genes <- unique(matched$gene[matched$allele_count > 0])
  pheno_genes <- intersect(affected_genes, names(PHENOCONVERSION_DRUGS))
  if (length(pheno_genes) > 0) {
    inhibitor_list <- sapply(pheno_genes, function(g) {
      paste0(g, " (common inhibitors: ", paste(PHENOCONVERSION_DRUGS[[g]], collapse = ", "), ")")
    })
    warnings <- c(warnings, list(list(
      type = "info",
      icon = "pills",
      title = "Phenoconversion Available for Detected Enzymes",
      message = paste0(
        "Enzymes in this patient's results can be affected by co-administered drugs: ",
        paste(inhibitor_list, collapse = "; "), ". ",
        "To model phenoconversion (adjusted risk from drug-drug interactions), ",
        "select the patient's current medications under 'Concomitant Medications' ",
        "in the sidebar, then view the Drug Interactions tab."
      )
    )))
  }
  
  return(warnings)
}
