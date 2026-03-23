# PGx-ADR: Genotype-Based Adverse Drug Reaction Risk Calculator

**A novel pharmacogenomics tool for quantitative ADR risk prediction**

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://cran.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-1.7+-green.svg)](https://shiny.rstudio.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

PGx-ADR is an interactive web application that calculates **quantitative adverse drug reaction (ADR) risk scores** from patient genotype data. Unlike existing tools that provide qualitative guideline recommendations, PGx-ADR generates **numerical risk ratios** that quantify a patient's ADR susceptibility compared to the general population.

### What Makes This Novel?

| Existing Tools | PGx-ADR |
|----------------|---------|
| **PharmCAT**: Outputs CPIC guideline text ("consider alternative therapy") | Outputs numerical risk ratios (e.g., "4.82x population average") |
| **PGxRAG**: Retrieves guideline recommendations via LLM | Predicts quantitative risk probabilities |
| **shinyDeepDR**: Predicts drug efficacy (IC50) | Predicts ADR risk from germline variants |
| **DGANet**: Population-level drug-ADR associations | Individual patient-level risk scoring |
| Binary output (actionable/not) | Continuous risk ratio with evidence weighting |

## Key Features

- **VCF Upload**: Parse standard VCF files to extract pharmacogenomic variants
- **Manual Entry**: Enter rsIDs with genotypes (e.g., `rs1799853:0/1`)
- **Evidence-Weighted Risk Scoring**: Novel multiplicative model using CPIC/PharmGKB evidence levels
- **Multi-Gene Aggregation**: Combines risk across all pharmacogenes affecting each drug
- **Population Context**: Compares variant frequencies across 6 ancestry groups (gnomAD-derived)
- **Interactive Dashboard**: Risk heatmaps, bar charts, gene-level summaries
- **Database Explorer**: Browse all 63 curated variant-drug-ADR associations
- **CSV Export**: Download results for downstream analysis

## Algorithm

### Evidence-Weighted Geometric Mean Risk Score

For each drug, the combined risk ratio is calculated as:

```
Combined_Risk = exp( Σ [log(variant_risk_i) × w_i] / Σ w_i )
```

Where:
- `variant_risk_i = base_risk_ratio ^ allele_count` (per-allele multiplicative model)
- `w_i` = evidence weight from PharmGKB/CPIC evidence levels:

| Evidence Level | Weight | Description |
|---------------|--------|-------------|
| 1A | 1.0 | CPIC guideline, FDA label |
| 1B | 0.9 | Strong evidence, replicated |
| 2A | 0.7 | VIP gene, PharmGKB annotation |
| 2B | 0.5 | Moderate evidence |
| 3  | 0.3 | Single significant study |
| 4  | 0.1 | Case report |

### Risk Categories

| Risk Ratio | Category |
|------------|----------|
| ≥ 5.0x | Very High |
| ≥ 2.5x | High |
| ≥ 1.5x | Moderate |
| ≥ 1.2x | Slightly Elevated |
| < 1.2x | Normal |

## Data Sources

The curated dataset includes **63 variant-drug-ADR associations** across **19 pharmacogenes** from:

- **CPIC Guidelines**: Clinical Pharmacogenetics Implementation Consortium
- **PharmGKB**: Level 1A/1B clinical annotations
- **FDA Drug Labels**: Pharmacogenomic biomarker information
- **gnomAD**: Population allele frequencies (Global, EUR, AFR, EAS, SAS, AMR)

### Covered Genes (19 pharmacogenes)

CYP2D6, CYP2C19, CYP2C9, CYP2B6, CYP3A5, CYP4F2, VKORC1, SLCO1B1, ABCB1, TPMT, NUDT15, DPYD, UGT1A1, G6PD, IFNL3, MTHFR, SLC19A1, OPRM1, PTPN22

### ADR Types Covered

- **Life-threatening**: SJS/TEN, hemolytic anemia, severe myelosuppression, fatal toxicity
- **Severe**: Major bleeding, CNS toxicity, rhabdomyolysis, neutropenia
- **Moderate**: Myopathy, nephrotoxicity, cardiotoxicity, QT prolongation
- **Clinical Impact**: Treatment failure, reduced efficacy, dose adjustment needed

## Installation

### Prerequisites

- R ≥ 4.0
- RStudio (recommended)

### Install Dependencies

```r
install.packages(c(
  "shiny",
  "shinydashboard",
  "DT",
  "ggplot2",
  "dplyr",
  "tidyr",
  "plotly"
))
```

### Run Locally

```r
# Clone the repository
# git clone https://github.com/[username]/pgx-adr-calculator.git

# Set working directory
setwd("pgx-adr-calculator")

# Run the app
shiny::runApp()
```

### Deploy to shinyapps.io

```r
library(rsconnect)
rsconnect::deployApp()
```

## Usage

### Option 1: Upload VCF File
1. Select "Upload VCF" in the sidebar
2. Upload your VCF file (must include rsIDs in the ID column)
3. Click "Analyze"

### Option 2: Manual Entry
1. Select "Manual Entry" in the sidebar
2. Enter variants in format: `rs1799853:0/1` (one per line)
3. Click "Analyze"

### Option 3: Example Data
1. Select "Load Example" (default)
2. Click "Analyze" to see results with 5 sample variants

### Example Output

```
Drug          | Risk Ratio | Category | Primary ADR
--------------|------------|----------|---------------------------
Warfarin      | 4.82x      | High     | Major Bleeding
Simvastatin   | 2.12x      | Moderate | Myopathy/Rhabdomyolysis
Clopidogrel   | 1.76x      | Moderate | Major Cardiovascular Events
Codeine       | 1.87x      | Moderate | Respiratory Depression
```

## Project Structure

```
pgx-adr-calculator/
├── app.R                           # Main Shiny application
├── R/
│   └── risk_calculator.R           # Core risk scoring algorithm
├── data/
│   ├── pgx_adr_associations.csv    # 63 variant-drug-ADR associations
│   └── example_sample.vcf          # Demo VCF file
├── scripts/
│   └── download_pharmgkb_data.py   # Data curation documentation
├── www/                            # Static assets
├── LICENSE
└── README.md
```

## Limitations & Disclaimers

⚠️ **For Research Use Only** — Not intended for clinical decision-making without professional interpretation.

- Risk scores are based on published associations and may not capture all genetic factors
- Population frequencies are derived from gnomAD and may not represent all ancestries equally
- Drug-drug interactions are not currently modeled
- Environmental and non-genetic factors are not included
- Novel variants not in the database will not be scored
- Star allele calling for complex genes (e.g., CYP2D6 CNVs) is simplified

## Future Directions

- [ ] Integration with gnomAD v4 for expanded multi-ancestry population frequencies
- [ ] Drug-drug-gene interaction module
- [ ] FAERS signal integration for expanded ADR coverage
- [ ] API endpoint for programmatic access
- [ ] FHIR/HL7 compatibility for EHR integration
- [ ] Validation against UK Biobank ADR outcomes

## Citation

If you use PGx-ADR in your research, please cite:

```
PGx-ADR: A Quantitative Pharmacogenomic Risk Calculator for Adverse Drug Reactions.
GitHub: https://github.com/[username]/pgx-adr-calculator
```

## References

1. Whirl-Carrillo M, et al. (2021). PharmGKB: An Integrated Resource of Pharmacogenomic Knowledge. *Clin Pharmacol Ther*.
2. Relling MV, Klein TE. (2011). CPIC: Clinical Pharmacogenetics Implementation Consortium. *Clin Pharmacol Ther*.
3. Pirmohamed M. (2023). Pharmacogenomics: Current status and future perspectives. *Nature Reviews Genetics*.
4. Magavern EF, et al. (2025). Pharmacogenetics and adverse drug reports: Insights from a UK national pharmacovigilance database. *PLoS Medicine*.
5. Swen JJ, et al. (2023). A 12-gene pharmacogenetic panel to prevent adverse drug reactions (PREPARE trial). *Lancet*.

## License

MIT License — See [LICENSE](LICENSE) for details.
