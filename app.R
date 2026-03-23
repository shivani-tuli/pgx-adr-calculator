# ============================================================================
# PGx-ADR: Genotype-Based Adverse Drug Reaction Risk Calculator
# ============================================================================
# A novel pharmacogenomics tool that calculates QUANTITATIVE ADR risk scores
# from patient genotype data. Unlike PharmCAT (qualitative guidelines) or
# PGxRAG (guideline retrieval), this tool outputs numerical risk ratios.
#
# Author: Shivani
# License: MIT
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)

# Source the risk calculator
source("R/risk_calculator.R")

# Load PGx-ADR associations database
pgx_data <- read.csv("data/pgx_adr_associations.csv", stringsAsFactors = FALSE)

# Load diplotype risk database (for PharmCAT mode)
diplotype_data <- read.csv("data/diplotype_risk_associations.csv", stringsAsFactors = FALSE)

# Load CYP inhibitor/inducer database (for phenoconversion)
inhibitor_data <- read.csv("data/cyp_inhibitors_inducers.csv", stringsAsFactors = FALSE)
all_concomitant_drugs <- sort(unique(inhibitor_data$drug_name))

# Check PharmCAT availability at startup
pharmcat_status <- check_pharmcat_available()

# ---- Example Data ----
example_variants <- data.frame(
  rsid = c("rs1799853", "rs4149056", "rs9923231", "rs4244285", "rs3892097"),
  genotype = c("0/1", "0/1", "1/1", "0/1", "0/1"),
  allele_count = c(1, 1, 2, 1, 1),
  stringsAsFactors = FALSE
)

# ============================================================================
# UI
# ============================================================================
ui <- dashboardPage(
  skin = "black",

  dashboardHeader(
    title = span(
      icon("dna"), " PGx-ADR Risk Calculator"
    ),
    titleWidth = 350
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "sidebar",
      menuItem("Risk Analysis", tabName = "analysis", icon = icon("chart-bar")),
      menuItem("Diplotype Report", tabName = "diplotypes", icon = icon("dna")),
      menuItem("Drug Interactions", tabName = "interactions", icon = icon("pills")),
      menuItem("Gene Summary", tabName = "genes", icon = icon("microscope")),
      menuItem("Population Context", tabName = "population", icon = icon("globe")),
      menuItem("Database Explorer", tabName = "database", icon = icon("database")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),

    # PharmCAT status badge
    tags$div(
      style = "padding: 5px 15px; font-size: 11px;",
      if (pharmcat_status$available) {
        tags$span(style = "color: #16a34a;", icon("check-circle"),
                  " PharmCAT Available")
      } else {
        tags$span(style = "color: #ca8a04;", icon("exclamation-circle"),
                  " PharmCAT Not Installed (rsID mode)")
      }
    ),

    tags$hr(),

    # Input section
    h4("  Input Genotype Data", style = "padding-left: 15px; color: #ecf0f1;"),

    radioButtons("input_method", "Method:",
      choices = c("Upload VCF" = "vcf",
                  "Manual Entry" = "manual",
                  "Load Example" = "example"),
      selected = "example"
    ),

    conditionalPanel(
      condition = "input.input_method == 'vcf'",
      fileInput("vcf_file", "Upload VCF File:",
                accept = c(".vcf", ".vcf.gz")),
      if (pharmcat_status$available) {
        checkboxInput("use_pharmcat", "Use PharmCAT Star Allele Calling",
                      value = TRUE)
      }
    ),

    conditionalPanel(
      condition = "input.input_method == 'manual'",
      textAreaInput("manual_variants", "Enter Variants:",
        placeholder = "rs1799853:0/1\nrs4149056:1/1\nrs9923231:0/1",
        rows = 5
      ),
      helpText("Format: rsID:genotype (one per line). Uses rsID mode.")
    ),

    actionButton("analyze_btn", "Analyze",
      icon = icon("play"),
      class = "btn-primary",
      style = "margin: 10px 15px; width: 85%;"
    ),

    tags$hr(),
    h4("  Concomitant Medications", style = "padding-left: 15px; color: #ecf0f1;"),
    selectizeInput("concomitant_meds", "Current Medications:",
      choices = all_concomitant_drugs,
      multiple = TRUE,
      options = list(
        placeholder = "Search and select drugs...",
        maxItems = 10
      )
    ),
    helpText("Select co-administered drugs to model phenoconversion ",
             "(CYP inhibitors/inducers).",
             style = "padding: 0 15px; font-size: 11px;"),

    tags$hr(),
    tags$div(
      style = "padding: 0 15px; font-size: 11px; color: #95a5a6;",
      icon("exclamation-triangle"),
      " For Research Use Only. Not for clinical decision-making."
    )
  ),

  dashboardBody(
    # Custom CSS
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #1a1a2e; }
      .box { background-color: #16213e; border-top: 3px solid #0f3460; color: #ecf0f1; }
      .box-header { color: #ecf0f1; }
      .box-body { color: #ecf0f1; }
      .info-box { background-color: #16213e; color: #ecf0f1; }
      .info-box-icon { background-color: #0f3460; }
      .small-box { border-radius: 8px; }
      .skin-black .main-header .logo { background-color: #0f3460; }
      .skin-black .main-header .navbar { background-color: #16213e; }
      .skin-black .main-sidebar { background-color: #0a0a1a; }
      .dataTables_wrapper { color: #ecf0f1; }
      table.dataTable { color: #ecf0f1 !important; }
      table.dataTable thead th { color: #ecf0f1 !important; background-color: #0f3460 !important; }
      table.dataTable tbody td { color: #ecf0f1 !important; }
      .dataTables_info, .dataTables_length, .dataTables_filter { color: #ecf0f1 !important; }
      .risk-very-high { color: #dc2626; font-weight: bold; }
      .risk-high { color: #ea580c; font-weight: bold; }
      .risk-moderate { color: #ca8a04; font-weight: bold; }
      .risk-elevated { color: #65a30d; }
      .risk-normal { color: #16a34a; }
      .nav-tabs-custom > .tab-content { background: #16213e; }
      h3, h4 { color: #ecf0f1; }
      .tab-pane { color: #ecf0f1; }
      .clinical-warning { border-radius: 8px; padding: 12px 16px; margin-bottom: 10px; }
      .warning-critical { background-color: rgba(220,38,38,0.15); border-left: 4px solid #dc2626; color: #fca5a5; }
      .warning-warning { background-color: rgba(234,88,12,0.15); border-left: 4px solid #ea580c; color: #fdba74; }
      .warning-info { background-color: rgba(59,130,246,0.15); border-left: 4px solid #3b82f6; color: #93c5fd; }
      .warning-title { font-weight: bold; font-size: 14px; margin-bottom: 4px; }
      .warning-message { font-size: 12px; line-height: 1.5; }
      .hla-disclaimer { background-color: rgba(139,92,246,0.15); border-left: 4px solid #8b5cf6; border-radius: 8px; padding: 12px 16px; color: #c4b5fd; margin-bottom: 10px; }
    "))),

    tabItems(
      # ---- Risk Analysis Tab ----
      tabItem(tabName = "analysis",
        fluidRow(
          valueBoxOutput("total_variants_box", width = 3),
          valueBoxOutput("high_risk_drugs_box", width = 3),
          valueBoxOutput("genes_affected_box", width = 3),
          valueBoxOutput("max_risk_box", width = 3)
        ),

        fluidRow(
          box(
            title = "ADR Risk Scores by Drug",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("risk_bar_chart", height = "400px")
          )
        ),

        fluidRow(
          box(
            title = "Clinical Warnings & Caveats",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            uiOutput("clinical_warnings_panel")
          )
        ),

        fluidRow(
          box(
            title = "Risk Heatmap by Drug Class",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            plotlyOutput("risk_heatmap", height = "350px")
          ),
          box(
            title = "ADR Severity Distribution",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            plotlyOutput("severity_chart", height = "350px")
          )
        ),

        fluidRow(
          box(
            title = "Detailed Risk Results",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("risk_table"),
            downloadButton("download_results", "Export as CSV",
              style = "margin-top: 10px;")
          )
        )
      ),

      # ---- Diplotype Report Tab ----
      tabItem(tabName = "diplotypes",
        fluidRow(
          box(
            title = "PharmCAT Star Allele Calling",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            uiOutput("pharmcat_mode_indicator"),
            tags$hr(),
            h4("Diplotype Calls by Gene"),
            helpText("Star allele diplotypes called by PharmCAT's Named Allele Matcher,
                      with CPIC-assigned phenotypes and activity scores. These calls
                      replace the simplified rsID-based genotyping with validated
                      haplotype-level resolution."),
            DT::dataTableOutput("diplotype_table")
          )
        ),
        fluidRow(
          box(
            title = "Diplotype-Based Risk Scores",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            helpText("Risk scores calculated from diplotype-level phenotype assignments,
                      mapped to the curated CPIC drug-phenotype-risk database.
                      This provides more accurate risk stratification than per-SNP scoring."),
            DT::dataTableOutput("diplotype_risk_table")
          )
        )
      ),

      # ---- Drug Interactions Tab (Phenoconversion) ----
      tabItem(tabName = "interactions",
        fluidRow(
          box(
            title = "Phenoconversion: Drug-Gene Interactions",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            helpText("Phenoconversion models how co-administered CYP inhibitors and inducers ",
                     "shift a patient's genetic metabolizer phenotype. For example, a CYP2D6 ",
                     "Normal Metabolizer taking fluoxetine (strong CYP2D6 inhibitor) effectively ",
                     "becomes a Poor Metabolizer, dramatically changing ADR risk."),
            uiOutput("phenoconversion_summary"),
            tags$hr(),
            h4("Enzyme-Level Interaction Details"),
            DT::dataTableOutput("interaction_table")
          )
        ),
        fluidRow(
          box(
            title = "Adjusted vs Genetic Risk Scores",
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            helpText("Risk scores recalculated using phenoconverted (adjusted) phenotypes. ",
                     "Compare with genetic-only scores in the Risk Analysis tab."),
            DT::dataTableOutput("adjusted_risk_table")
          )
        )
      ),

      # ---- Gene Summary Tab ----
      tabItem(tabName = "genes",
        fluidRow(
          box(
            title = "Pharmacogene Impact Summary",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("gene_impact_chart", height = "400px")
          )
        ),
        fluidRow(
          box(
            title = "Gene-Level Details",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("gene_table")
          )
        )
      ),

      # ---- Population Context Tab ----
      tabItem(tabName = "population",
        fluidRow(
          box(
            title = "Variant Allele Frequencies Across Populations",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            helpText("Comparison of patient variant allele frequencies against
                      global population databases (gnomAD-derived estimates).
                      Allele frequencies vary significantly across ancestries,
                      which affects ADR risk calibration."),
            tags$div(
              class = "hla-disclaimer",
              tags$div(class = "warning-title", icon("info-circle"), " Interpreting Population Frequencies"),
              tags$div(class = "warning-message",
                "Higher allele frequency means the variant is more common in that population. ",
                "A heterozygous carrier (0/1) of a common variant (e.g., MAF=0.15) has different ",
                "clinical implications than a homozygous carrier (1/1) of a rare variant (MAF=0.01). ",
                "Population context helps calibrate whether a patient's genotype is typical or unusual ",
                "for their ancestry, which may affect the clinical urgency of pharmacogenomic testing."
              )
            ),
            plotlyOutput("population_chart", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Population Frequency Table",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("population_table")
          )
        )
      ),

      # ---- Database Explorer Tab ----
      tabItem(tabName = "database",
        fluidRow(
          box(
            title = "PGx-ADR Association Database",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            helpText(paste0("Curated database of ", nrow(pgx_data),
              " variant-drug-ADR associations from PharmGKB, CPIC guidelines,",
              " and FDA pharmacogenomic labels.")),
            DT::dataTableOutput("database_table")
          )
        )
      ),

      # ---- About Tab ----
      tabItem(tabName = "about",
        fluidRow(
          box(
            title = "About PGx-ADR Risk Calculator",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            h3("Novel Approach"),
            p("PGx-ADR is the first tool to provide QUANTITATIVE adverse drug
              reaction risk scores from patient genotype data. Existing tools
              provide qualitative outputs:"),
            tags$ul(
              tags$li(strong("PharmCAT:"), " Outputs CPIC guideline text
                      (e.g., 'consider alternative therapy')"),
              tags$li(strong("PGxRAG:"), " Retrieves guideline recommendations
                      via LLM"),
              tags$li(strong("shinyDeepDR:"), " Predicts drug efficacy (IC50),
                      not ADR risk")
            ),
            p(strong("PGx-ADR outputs numerical risk ratios"), " (e.g.,
              'Warfarin bleeding risk: 4.82x population average'), enabling
              quantitative clinical decision support."),

            h3("Algorithm"),
            p("For each drug, the combined risk ratio is calculated as an
              evidence-weighted geometric mean:"),
            tags$pre("Combined_Risk = exp( sum(log(variant_risk_i * w_i)) / sum(w_i) )"),
            p("Where variant_risk is the base risk ratio raised to the power of
              the allele count, and weights are derived from PharmGKB/CPIC
              evidence levels (1A=1.0 through 4=0.1)."),

            h3("Data Sources"),
            tags$ul(
              tags$li("CPIC Guidelines (Clinical Pharmacogenetics Implementation Consortium)"),
              tags$li(strong("PharmGKB REST API"), " — variant annotations and evidence levels ",
                     "programmatically pulled via ", tags$code("api.pharmgkb.org")),
              tags$li("FDA Pharmacogenomic Drug Labels"),
              tags$li("Population frequencies from gnomAD"),
              tags$li(strong("Meta-analysis-derived risk ratios"), " — pooled effect sizes from ",
                     "published systematic meta-analyses (10 unique PMIDs)")
            ),
            uiOutput("data_freshness_indicator"),
            actionButton("refresh_pharmgkb", "Refresh PharmGKB Data",
              icon = icon("sync"),
              style = "margin: 5px 0 15px 0;"),
            tags$div(
              style = "font-size: 11px; color: #6b7280; margin-bottom: 15px;",
              "Data from PharmGKB is subject to the ",
              tags$a(href = "http://creativecommons.org/licenses/by-sa/4.0/",
                     "CC-BY-SA 4.0 International License", target = "_blank"),
              ". See ", tags$a(href = "https://www.pharmgkb.org/page/dataUsagePolicy",
                              "PharmGKB Data Usage Policy", target = "_blank"), "."
            ),

            h3("External Validation: FDA FAERS"),
            p("Risk ratios were validated against real-world adverse event data from the ",
              strong("FDA Adverse Event Reporting System (FAERS)"),
              " via the openFDA API. Proportional Reporting Ratios (PRRs) were computed ",
              "for each drug-ADR pair and correlated with predicted risk ratios."),
            uiOutput("faers_validation_summary"),
            DT::dataTableOutput("faers_validation_table"),
            tags$div(
              style = "font-size: 11px; color: #6b7280; margin-top: 10px;",
              "PRRs are pharmacovigilance signal detection metrics, not causal risk ratios. ",
              "A positive correlation validates alignment between predicted risk and ",
              "real-world ADR reporting patterns."
            ),

            h3("Limitations"),
            tags$ul(
              tags$li("For Research Use Only - not for clinical decisions without professional review"),
              tags$li("Population frequencies primarily from European-ancestry datasets"),
              tags$li(strong("HLA alleles not covered:"), " Clinically critical HLA-based tests ",
                      "(HLA-B*57:01/Abacavir, HLA-B*15:02/Carbamazepine, HLA-A*31:01/Carbamazepine) ",
                      "cannot be captured by rsID-based VCF genotyping. These tests require ",
                      "dedicated HLA typing assays."),
              tags$li(strong("Phenoconversion modeled:"), " The Drug Interactions module accounts for ",
                      "co-administered CYP inhibitors/inducers shifting metabolizer phenotypes. ",
                      "Select concomitant medications in the sidebar to activate."),
              tags$li(strong("Dose not modeled:"), " Risk scores do not account for drug dose. ",
                      "Statin myopathy risk, for example, is strongly dose-dependent (highest at 80mg)."),
              tags$li(strong("X-linked genes:"), " G6PD is X-linked; hemizygous males are fully affected. ",
                      "The tool does not model sex-specific inheritance."),
              tags$li(strong("CYP2D6 structural variation:"), " Gene deletions, duplications, and hybrid ",
                      "alleles cannot be detected from rsIDs alone."),
              tags$li("Environmental and non-genetic factors not included"),
              tags$li("Novel variants not in database will not be scored")
            ),

            h3("Citation"),
            tags$pre("PGx-ADR: A Quantitative Pharmacogenomic Risk Calculator
for Adverse Drug Reactions. 2026.
GitHub: https://github.com/shivani-tuli/pgx-adr-calculator")
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {

  # Reactive: patient variants
  patient_data <- reactiveVal(NULL)

  # On analyze button click
  observeEvent(input$analyze_btn, {
    variants <- NULL

    if (input$input_method == "vcf" && !is.null(input$vcf_file)) {
      variants <- parse_vcf(input$vcf_file$datapath)
    } else if (input$input_method == "manual" && nchar(input$manual_variants) > 0) {
      variants <- parse_manual_entry(input$manual_variants)
    } else if (input$input_method == "example") {
      variants <- example_variants
    }

    patient_data(variants)
  })

  # Auto-load example on start
  observe({
    if (is.null(patient_data())) {
      patient_data(example_variants)
    }
  })

  # Reactive: risk scores
  risk_scores <- reactive({
    req(patient_data())
    calculate_risk_scores(patient_data(), pgx_data)
  })

  # Reactive: gene summary
  gene_summary <- reactive({
    req(patient_data())
    summarize_by_gene(patient_data(), pgx_data)
  })

  # Reactive: population context
  pop_context <- reactive({
    req(patient_data())
    population_context(patient_data(), pgx_data)
  })

  # Reactive: clinical warnings
  clinical_warnings <- reactive({
    req(patient_data())
    generate_clinical_warnings(patient_data(), pgx_data, risk_scores())
  })

  # ---- Clinical Warnings Panel ----
  output$clinical_warnings_panel <- renderUI({
    warns <- clinical_warnings()
    if (length(warns) == 0) {
      return(tags$div(
        style = "color: #16a34a; padding: 8px;",
        icon("check-circle"), " No clinical warnings for this genotype profile."
      ))
    }

    warning_divs <- lapply(warns, function(w) {
      css_class <- paste0("clinical-warning warning-", w$type)
      tags$div(
        class = css_class,
        tags$div(class = "warning-title", icon(w$icon), " ", w$title),
        tags$div(class = "warning-message", w$message)
      )
    })

    # Add HLA disclaimer at the end
    hla_div <- tags$div(
      class = "hla-disclaimer",
      tags$div(class = "warning-title", icon("exclamation-triangle"), " HLA Alleles Not Assessed"),
      tags$div(class = "warning-message",
        "Clinically critical HLA-drug associations (HLA-B*57:01/Abacavir hypersensitivity, ",
        "HLA-B*15:02/Carbamazepine SJS/TEN, HLA-A*31:01/Carbamazepine DRESS) cannot be assessed ",
        "from rsID-based VCF data. Dedicated HLA typing is required for these tests."
      )
    )

    do.call(tagList, c(warning_divs, list(hla_div)))
  })

  # ---- PharmCAT Integration ----
  # Store the uploaded VCF path for PharmCAT
  uploaded_vcf_path <- reactiveVal(NULL)

  # Track when VCF is uploaded — save to temp location for PharmCAT
  observeEvent(input$vcf_file, {
    if (!is.null(input$vcf_file)) {
      uploaded_vcf_path(input$vcf_file$datapath)
    }
  })

  # Reactive: PharmCAT results (only runs when PharmCAT mode is active)
  pharmcat_results <- reactive({
    vcf_path <- uploaded_vcf_path()
    req(vcf_path)
    req(input$input_method == "vcf")

    # Only run PharmCAT if the toggle is checked and PharmCAT is available
    use_pc <- if (!is.null(input$use_pharmcat)) input$use_pharmcat else FALSE
    req(use_pc && pharmcat_status$available)

    tryCatch({
      tsv_path <- run_pharmcat(vcf_path)
      parse_pharmcat_results(tsv_path)
    }, error = function(e) {
      showNotification(
        paste("PharmCAT error:", e$message, "- Falling back to rsID mode"),
        type = "warning", duration = 10
      )
      NULL
    })
  })

  # Reactive: diplotype-based risk scores
  diplotype_scores <- reactive({
    pc_results <- pharmcat_results()
    req(pc_results)
    calculate_diplotype_risk(pc_results, diplotype_data)
  })

  # PharmCAT mode indicator
  output$pharmcat_mode_indicator <- renderUI({
    use_pc <- if (!is.null(input$use_pharmcat)) input$use_pharmcat else FALSE

    if (use_pc && pharmcat_status$available) {
      pc_results <- tryCatch(pharmcat_results(), error = function(e) NULL)

      if (!is.null(pc_results) && nrow(pc_results) > 0) {
        tags$div(
          class = "clinical-warning warning-info",
          tags$div(class = "warning-title",
                   icon("check-circle"),
                   sprintf(" PharmCAT Mode Active — %d gene(s) called",
                           nrow(pc_results))),
          tags$div(class = "warning-message",
                   "Star allele diplotypes were called by PharmCAT 3.2.0. ",
                   "Risk scores use diplotype-level phenotype mapping (CPIC-aligned).")
        )
      } else {
        tags$div(
          class = "clinical-warning warning-warning",
          tags$div(class = "warning-title", icon("exclamation-triangle"),
                   " PharmCAT: No Results"),
          tags$div(class = "warning-message",
                   "PharmCAT did not return gene calls. This may be because the VCF ",
                   "does not contain PGx positions or is not GRCh38-aligned. ",
                   "Falling back to rsID mode.")
        )
      }
    } else {
      tags$div(
        class = "clinical-warning warning-warning",
        tags$div(class = "warning-title", icon("info-circle"),
                 " rsID Mode (Basic)"),
        tags$div(class = "warning-message",
                 if (!pharmcat_status$available) {
                   "PharmCAT is not installed. Using simplified rsID-based genotyping. "
                 } else {
                   "PharmCAT is disabled. Using simplified rsID-based genotyping. "
                 },
                 "Enable PharmCAT for validated star allele calling with diplotype-level risk.")
      )
    }
  })

  # Diplotype calls table
  output$diplotype_table <- DT::renderDataTable({
    pc_results <- tryCatch(pharmcat_results(), error = function(e) NULL)

    if (is.null(pc_results) || nrow(pc_results) == 0) {
      return(DT::datatable(
        data.frame(Message = "No PharmCAT results available. Upload a GRCh38-aligned VCF with PharmCAT enabled."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    diplo_summary <- summarize_diplotypes(pc_results, diplotype_data)

    display_df <- data.frame(
      Gene = diplo_summary$gene,
      Diplotype = diplo_summary$diplotype,
      Phenotype = diplo_summary$phenotype,
      `Activity Score` = ifelse(is.na(diplo_summary$activity_score),
                                "N/A",
                                as.character(diplo_summary$activity_score)),
      `Drugs Affected` = diplo_summary$n_drugs,
      `Max Risk Ratio` = ifelse(is.na(diplo_summary$max_risk),
                                "—",
                                sprintf("%.1fx", diplo_summary$max_risk)),
      check.names = FALSE
    )

    DT::datatable(
      display_df,
      options = list(pageLength = 25, dom = "ft"),
      rownames = FALSE
    ) %>%
      DT::formatStyle("Phenotype",
        backgroundColor = DT::styleEqual(
          c("Poor Metabolizer", "Intermediate Metabolizer",
            "Normal Metabolizer", "Rapid Metabolizer",
            "Ultrarapid Metabolizer"),
          c("rgba(220,38,38,0.2)", "rgba(234,88,12,0.2)",
            "rgba(22,163,74,0.2)", "rgba(59,130,246,0.2)",
            "rgba(139,92,246,0.2)")
        )
      )
  })

  # Diplotype-based risk results table
  output$diplotype_risk_table <- DT::renderDataTable({
    ds <- tryCatch(diplotype_scores(), error = function(e) NULL)

    if (is.null(ds) || nrow(ds) == 0) {
      return(DT::datatable(
        data.frame(Message = "No diplotype-based risk scores available."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    display_df <- data.frame(
      Drug = ds$drug_name,
      Class = ds$drug_class,
      `Risk Ratio` = sprintf("%.2fx", ds$combined_risk_ratio),
      Category = ds$risk_category,
      ADR = ds$primary_adr,
      Diplotypes = ds$diplotypes,
      Evidence = ds$max_evidence,
      check.names = FALSE
    )

    DT::datatable(
      display_df,
      options = list(pageLength = 25, dom = "ft"),
      rownames = FALSE
    )
  })

  # ---- Value Boxes ----
  output$total_variants_box <- renderValueBox({
    n <- if (!is.null(patient_data())) nrow(patient_data()) else 0
    valueBox(n, "Variants Analyzed", icon = icon("dna"), color = "blue")
  })

  output$high_risk_drugs_box <- renderValueBox({
    rs <- risk_scores()
    n <- if (!is.null(rs)) sum(rs$risk_category %in% c("High", "Very High")) else 0
    valueBox(n, "High-Risk Drugs", icon = icon("exclamation-triangle"),
             color = if (n > 0) "red" else "green")
  })

  output$genes_affected_box <- renderValueBox({
    gs <- gene_summary()
    n <- if (!is.null(gs)) nrow(gs) else 0
    valueBox(n, "Pharmacogenes Hit", icon = icon("microscope"), color = "purple")
  })

  output$max_risk_box <- renderValueBox({
    rs <- risk_scores()
    max_r <- if (!is.null(rs) && nrow(rs) > 0) {
      sprintf("%.1fx", max(rs$combined_risk_ratio))
    } else { "N/A" }
    valueBox(max_r, "Highest Risk Ratio", icon = icon("arrow-up"),
             color = "orange")
  })

  # ---- Risk Bar Chart ----
  output$risk_bar_chart <- renderPlotly({
    rs <- risk_scores()
    req(rs, nrow(rs) > 0)

    rs <- rs %>% arrange(combined_risk_ratio)
    rs$drug_name <- factor(rs$drug_name, levels = rs$drug_name)

    p <- ggplot(rs, aes(x = drug_name, y = combined_risk_ratio,
                         fill = risk_category,
                         text = paste0(
                           "Drug: ", drug_name,
                           "\nRisk Ratio: ", sprintf("%.2fx", combined_risk_ratio),
                           "\nCategory: ", risk_category,
                           "\nPrimary ADR: ", primary_adr,
                           "\nGenes: ", genes_involved
                         ))) +
      geom_col() +
      geom_hline(yintercept = 1, linetype = "dashed", color = "#6b7280") +
      coord_flip() +
      scale_fill_manual(values = c(
        "Very High" = "#dc2626",
        "High" = "#ea580c",
        "Moderate" = "#ca8a04",
        "Slightly Elevated" = "#65a30d",
        "Normal" = "#16a34a"
      )) +
      labs(x = NULL, y = "Risk Ratio (vs. population average)", fill = "Risk Level") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "#16213e", color = NA),
        panel.background = element_rect(fill = "#16213e", color = NA),
        text = element_text(color = "#ecf0f1"),
        axis.text = element_text(color = "#ecf0f1"),
        legend.background = element_rect(fill = "#16213e"),
        legend.text = element_text(color = "#ecf0f1"),
        panel.grid.major = element_line(color = "#2a2a4a"),
        panel.grid.minor = element_blank()
      )

    ggplotly(p, tooltip = "text") %>%
      layout(paper_bgcolor = "#16213e", plot_bgcolor = "#16213e",
             font = list(color = "#ecf0f1"))
  })

  # ---- Risk Heatmap ----
  output$risk_heatmap <- renderPlotly({
    rs <- risk_scores()
    req(rs, nrow(rs) > 0)

    p <- ggplot(rs, aes(x = drug_class, y = drug_name,
                         fill = combined_risk_ratio,
                         text = paste0(drug_name, ": ", sprintf("%.2fx", combined_risk_ratio)))) +
      geom_tile(color = "#0a0a1a") +
      scale_fill_gradient2(low = "#16a34a", mid = "#ca8a04", high = "#dc2626",
                           midpoint = 2.5, name = "Risk Ratio") +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "#16213e", color = NA),
        panel.background = element_rect(fill = "#16213e", color = NA),
        text = element_text(color = "#ecf0f1", size = 9),
        axis.text.x = element_text(color = "#ecf0f1", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "#ecf0f1"),
        legend.background = element_rect(fill = "#16213e"),
        legend.text = element_text(color = "#ecf0f1")
      )

    ggplotly(p, tooltip = "text") %>%
      layout(paper_bgcolor = "#16213e", plot_bgcolor = "#16213e",
             font = list(color = "#ecf0f1"))
  })

  # ---- Severity Chart ----
  output$severity_chart <- renderPlotly({
    rs <- risk_scores()
    req(rs, nrow(rs) > 0)

    sev_counts <- rs %>%
      count(adr_severity) %>%
      mutate(adr_severity = factor(adr_severity,
        levels = c("Life-threatening", "Severe", "Moderate", "Clinical Impact")))

    p <- ggplot(sev_counts, aes(x = adr_severity, y = n, fill = adr_severity)) +
      geom_col() +
      scale_fill_manual(values = c(
        "Life-threatening" = "#dc2626",
        "Severe" = "#ea580c",
        "Moderate" = "#ca8a04",
        "Clinical Impact" = "#3b82f6"
      )) +
      labs(x = NULL, y = "Number of Drugs", fill = "Severity") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "#16213e", color = NA),
        panel.background = element_rect(fill = "#16213e", color = NA),
        text = element_text(color = "#ecf0f1"),
        axis.text = element_text(color = "#ecf0f1"),
        legend.position = "none",
        panel.grid.major = element_line(color = "#2a2a4a"),
        panel.grid.minor = element_blank()
      )

    ggplotly(p) %>%
      layout(paper_bgcolor = "#16213e", plot_bgcolor = "#16213e",
             font = list(color = "#ecf0f1"))
  })

  # ---- Risk Results Table ----
  output$risk_table <- DT::renderDataTable({
    rs <- risk_scores()
    req(rs, nrow(rs) > 0)

    display_df <- rs %>%
      select(
        Drug = drug_name,
        Class = drug_class,
        `Risk Ratio` = combined_risk_ratio,
        Category = risk_category,
        `Primary ADR` = primary_adr,
        Severity = adr_severity,
        `# Variants` = n_variants,
        Genes = genes_involved,
        Evidence = max_evidence
      ) %>%
      mutate(`Risk Ratio` = sprintf("%.2fx", `Risk Ratio`))

    DT::datatable(display_df,
      options = list(
        pageLength = 20,
        dom = 'frtip',
        initComplete = DT::JS(
          "function(settings, json) {",
          "  $(this.api().table().container()).css({'background-color': '#16213e'});",
          "}"
        )
      ),
      rownames = FALSE
    )
  })

  # ---- Gene Impact Chart ----
  output$gene_impact_chart <- renderPlotly({
    gs <- gene_summary()
    req(gs, nrow(gs) > 0)

    gs <- gs %>% arrange(n_drugs)
    gs$gene <- factor(gs$gene, levels = gs$gene)

    p <- ggplot(gs, aes(x = gene, y = n_drugs, fill = max_risk,
                         text = paste0(gene, "\nDrugs affected: ", n_drugs,
                                       "\nMax risk: ", sprintf("%.1fx", max_risk),
                                       "\nVariants: ", variants))) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(low = "#3b82f6", high = "#dc2626", name = "Max Risk") +
      labs(x = NULL, y = "Number of Drugs Affected") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "#16213e", color = NA),
        panel.background = element_rect(fill = "#16213e", color = NA),
        text = element_text(color = "#ecf0f1"),
        axis.text = element_text(color = "#ecf0f1"),
        legend.background = element_rect(fill = "#16213e"),
        legend.text = element_text(color = "#ecf0f1"),
        panel.grid.major = element_line(color = "#2a2a4a"),
        panel.grid.minor = element_blank()
      )

    ggplotly(p, tooltip = "text") %>%
      layout(paper_bgcolor = "#16213e", plot_bgcolor = "#16213e",
             font = list(color = "#ecf0f1"))
  })

  # ---- Gene Table ----
  output$gene_table <- DT::renderDataTable({
    gs <- gene_summary()
    req(gs, nrow(gs) > 0)

    display_df <- gs %>%
      select(
        Gene = gene,
        `# Variants` = n_variants,
        `# Drugs Affected` = n_drugs,
        `Max Risk Ratio` = max_risk,
        `Drugs` = drugs_affected,
        Variants = variants
      ) %>%
      mutate(`Max Risk Ratio` = sprintf("%.2fx", `Max Risk Ratio`))

    DT::datatable(display_df,
      options = list(pageLength = 20, dom = 'frtip'),
      rownames = FALSE
    )
  })

  # ---- Population Chart ----
  output$population_chart <- renderPlotly({
    pc <- pop_context()
    req(pc, nrow(pc) > 0)

    pop_long <- pc %>%
      select(rsid, gene, Global = global_maf, European = eur_af,
             African = afr_af, `East Asian` = eas_af,
             `South Asian` = sas_af, `Americas` = amr_af) %>%
      pivot_longer(cols = Global:`Americas`,
                   names_to = "population", values_to = "frequency") %>%
      mutate(label = paste0(rsid, " (", gene, ")"))

    p <- ggplot(pop_long, aes(x = label, y = frequency, fill = population)) +
      geom_col(position = "dodge") +
      scale_fill_brewer(palette = "Set2") +
      labs(x = NULL, y = "Allele Frequency", fill = "Population") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "#16213e", color = NA),
        panel.background = element_rect(fill = "#16213e", color = NA),
        text = element_text(color = "#ecf0f1"),
        axis.text.x = element_text(color = "#ecf0f1", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "#ecf0f1"),
        legend.background = element_rect(fill = "#16213e"),
        legend.text = element_text(color = "#ecf0f1"),
        panel.grid.major = element_line(color = "#2a2a4a"),
        panel.grid.minor = element_blank()
      )

    ggplotly(p) %>%
      layout(paper_bgcolor = "#16213e", plot_bgcolor = "#16213e",
             font = list(color = "#ecf0f1"))
  })

  # ---- Population Table ----
  output$population_table <- DT::renderDataTable({
    pc <- pop_context()
    req(pc, nrow(pc) > 0)

    display_df <- pc %>%
      select(
        Variant = rsid, Gene = gene, Genotype = genotype,
        `Allele Count` = allele_count,
        Global = global_maf, European = eur_af, African = afr_af,
        `East Asian` = eas_af, `South Asian` = sas_af, Americas = amr_af
      )

    DT::datatable(display_df,
      options = list(pageLength = 20, dom = 'frtip'),
      rownames = FALSE
    )
  })

  # ---- Database Explorer Table ----
  output$database_table <- DT::renderDataTable({
    display_df <- pgx_data %>%
      select(
        Variant = variant_rsid, Gene = gene, Allele = allele_name,
        Drug = drug_name, Class = drug_class,
        `Primary ADR` = primary_adr, Severity = adr_severity,
        `Risk Ratio` = risk_ratio, Evidence = evidence_level,
        `Global MAF` = global_maf
      )

    DT::datatable(display_df,
      options = list(pageLength = 25, dom = 'frtip'),
      filter = 'top',
      rownames = FALSE
    )
  })

  # ---- Download Handler ----
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("pgx_adr_risk_scores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      rs <- risk_scores()
      if (!is.null(rs)) {
        write.csv(rs, file, row.names = FALSE)
      }
    }
  )

  # ---- FAERS Validation Display ----
  faers_report <- tryCatch({
    jsonlite::fromJSON("data/faers_validation_report.json")
  }, error = function(e) NULL)

  output$faers_validation_summary <- renderUI({
    if (is.null(faers_report)) {
      return(tags$div(
        style = "padding: 10px; background: #fffbeb; border-left: 3px solid #f59e0b;
                 border-radius: 4px; margin: 10px 0;",
        icon("exclamation-triangle"),
        " FAERS validation report not found. Run ",
        tags$code("python3 scripts/validate_faers.py"), " to generate."
      ))
    }

    rho <- faers_report$spearman_rho
    n <- faers_report$n_valid_pairs
    interp <- faers_report$interpretation

    bg_color <- if (!is.null(rho) && rho > 0.3) "#f0fdf4" else "#fffbeb"
    border_color <- if (!is.null(rho) && rho > 0.3) "#16a34a" else "#f59e0b"
    icon_name <- if (!is.null(rho) && rho > 0.3) "check-circle" else "info-circle"

    tags$div(
      style = sprintf("padding: 12px; background: %s; border-left: 3px solid %s;
                       border-radius: 4px; margin: 10px 0;", bg_color, border_color),
      icon(icon_name),
      sprintf(" Spearman \u03c1 = %.3f (N=%d drug-ADR pairs). ", rho, n),
      interp,
      tags$br(),
      tags$span(style = "font-size: 11px; color: #6b7280;",
        sprintf("Validated: %s", faers_report$validation_date))
    )
  })

  output$faers_validation_table <- DT::renderDataTable({
    if (is.null(faers_report) || is.null(faers_report$results)) {
      return(DT::datatable(
        data.frame(Message = "No FAERS validation data available."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    res <- faers_report$results
    display_df <- data.frame(
      Drug = res$drug_name,
      ADR = substr(res$primary_adr, 1, 35),
      Gene = res$gene,
      `Predicted RR` = sprintf("%.2f", res$predicted_rr),
      `FAERS PRR` = sprintf("%.2f", res$faers_prr),
      `FAERS Reports` = format(res$faers_total_reports, big.mark = ","),
      `ADR Reports` = format(res$faers_adr_reports, big.mark = ","),
      check.names = FALSE
    )

    DT::datatable(display_df,
      options = list(dom = "ftp", pageLength = 10, order = list(list(4, "desc"))),
      rownames = FALSE
    )
  })

  # ---- Phenoconversion Module ----
  # Reactive: phenoconverted phenotypes
  phenoconverted_data <- reactive({
    meds <- input$concomitant_meds
    pc_results <- tryCatch(pharmcat_results(), error = function(e) NULL)

    if (is.null(pc_results) || nrow(pc_results) == 0) return(NULL)

    # Extract gene-phenotype pairs from PharmCAT results
    genetic_phenos <- data.frame(
      gene = pc_results$gene,
      phenotype = pc_results$phenotype,
      diplotype = pc_results$diplotype,
      stringsAsFactors = FALSE
    )

    apply_phenoconversion(genetic_phenos, meds, inhibitor_data)
  })

  # Phenoconversion summary indicator
  output$phenoconversion_summary <- renderUI({
    meds <- input$concomitant_meds

    if (length(meds) == 0) {
      return(tags$div(
        style = "padding: 10px; background: #f0f9ff; border-left: 3px solid #3b82f6;
                 border-radius: 4px; margin: 10px 0;",
        icon("info-circle"), " No concomitant medications entered. ",
        "Select medications in the sidebar to model phenoconversion."
      ))
    }

    pc_data <- tryCatch(phenoconverted_data(), error = function(e) NULL)

    if (is.null(pc_data)) {
      # No PharmCAT results — show interaction summary without phenotype shifts
      int_summary <- summarize_interactions(meds, inhibitor_data)
      if (is.null(int_summary) || nrow(int_summary) == 0) {
        return(tags$div(
          style = "padding: 10px; background: #f0fdf4; border-left: 3px solid #16a34a;
                   border-radius: 4px; margin: 10px 0;",
          icon("check-circle"),
          sprintf(" %d medication(s) entered — no known CYP interactions found.",
                  length(meds))
        ))
      }
      return(tags$div(
        style = "padding: 10px; background: #fffbeb; border-left: 3px solid #f59e0b;
                 border-radius: 4px; margin: 10px 0;",
        icon("exclamation-triangle"),
        sprintf(" %d medication(s) affect %d CYP enzyme(s). ",
                length(meds), nrow(int_summary)),
        "Upload a VCF with PharmCAT to see phenoconversion-adjusted risk."
      ))
    }

    n_shifted <- sum(pc_data$phenoconversion_active)
    if (n_shifted > 0) {
      shifted_genes <- pc_data$gene[pc_data$phenoconversion_active]
      tags$div(
        style = "padding: 10px; background: #fef2f2; border-left: 3px solid #dc2626;
                 border-radius: 4px; margin: 10px 0;",
        icon("exclamation-circle"),
        sprintf(" PHENOCONVERSION ACTIVE: %d gene(s) affected (%s). ",
                n_shifted, paste(shifted_genes, collapse = ", ")),
        "Risk scores have been adjusted. See adjusted scores below."
      )
    } else {
      tags$div(
        style = "padding: 10px; background: #f0fdf4; border-left: 3px solid #16a34a;
                 border-radius: 4px; margin: 10px 0;",
        icon("check-circle"),
        sprintf(" %d medication(s) entered — no phenotype shifts for this patient's genotype.",
                length(meds))
      )
    }
  })

  # Interaction details table
  output$interaction_table <- DT::renderDataTable({
    meds <- input$concomitant_meds

    if (length(meds) == 0) {
      return(DT::datatable(
        data.frame(Message = "No concomitant medications selected."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    pc_data <- tryCatch(phenoconverted_data(), error = function(e) NULL)
    int_summary <- summarize_interactions(meds, inhibitor_data, pc_data)

    if (is.null(int_summary) || nrow(int_summary) == 0) {
      return(DT::datatable(
        data.frame(Message = "No CYP interactions found for the selected medications."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    display_df <- data.frame(
      Enzyme = int_summary$enzyme,
      `Interacting Drugs` = int_summary$interacting_drugs,
      Effect = int_summary$strongest_effect,
      Potency = int_summary$strongest_potency,
      `Genetic Phenotype` = ifelse(is.na(int_summary$genetic_phenotype),
                                    "N/A (upload VCF)",
                                    int_summary$genetic_phenotype),
      `Adjusted Phenotype` = ifelse(is.na(int_summary$adjusted_phenotype),
                                     "N/A",
                                     int_summary$adjusted_phenotype),
      check.names = FALSE
    )

    DT::datatable(display_df, options = list(dom = "ft", pageLength = 20),
                  rownames = FALSE) %>%
      DT::formatStyle("Adjusted Phenotype",
        backgroundColor = DT::styleEqual(
          c("Poor Metabolizer", "Intermediate Metabolizer",
            "Normal Metabolizer"),
          c("rgba(220,38,38,0.2)", "rgba(234,88,12,0.2)",
            "rgba(22,163,74,0.2)")
        )
      )
  })

  # Adjusted risk scores table
  output$adjusted_risk_table <- DT::renderDataTable({
    pc_data <- tryCatch(phenoconverted_data(), error = function(e) NULL)

    if (is.null(pc_data) || !any(pc_data$phenoconversion_active)) {
      return(DT::datatable(
        data.frame(Message = "No phenoconversion adjustments active. Genetic risk scores remain unchanged."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    adj_scores <- calculate_adjusted_risk(pc_data, diplotype_data)

    if (is.null(adj_scores) || nrow(adj_scores) == 0) {
      return(DT::datatable(
        data.frame(Message = "No adjusted risk scores available."),
        options = list(dom = "t"), rownames = FALSE
      ))
    }

    display_df <- data.frame(
      Drug = adj_scores$drug_name,
      Class = adj_scores$drug_class,
      `Adjusted Risk` = sprintf("%.2fx", adj_scores$combined_risk_ratio),
      Category = adj_scores$risk_category,
      ADR = adj_scores$primary_adr,
      check.names = FALSE
    )

    DT::datatable(display_df, options = list(dom = "ft", pageLength = 20),
                  rownames = FALSE)
  })

  # ---- PharmGKB Data Freshness ----
  output$data_freshness_indicator <- renderUI({
    freshness <- check_data_freshness()

    color <- if (freshness$fresh) "#16a34a" else "#ca8a04"
    icon_name <- if (freshness$fresh) "check-circle" else "exclamation-triangle"
    status_text <- if (freshness$fresh) {
      sprintf("Data last updated: %s (%s days ago) via %s",
              substr(freshness$last_updated, 1, 10),
              freshness$age_days, freshness$source)
    } else {
      sprintf("Data may be stale (last updated: %s, %s days ago)",
              freshness$last_updated, freshness$age_days)
    }

    tags$div(
      style = sprintf("padding: 8px; margin: 5px 0; border-radius: 4px;
                        background: %s15; border-left: 3px solid %s;
                        font-size: 12px;", color, color),
      icon(icon_name), " ", status_text
    )
  })

  # Refresh PharmGKB data
  observeEvent(input$refresh_pharmgkb, {
    showNotification("Refreshing PharmGKB data... This may take ~60 seconds.",
                     type = "message", duration = 5)

    result <- trigger_data_update()

    if (result$success) {
      showNotification("PharmGKB data refreshed successfully! Reload the app to use updated data.",
                       type = "message", duration = 10)
    } else {
      showNotification(paste("PharmGKB refresh failed:", result$message),
                       type = "error", duration = 10)
    }
  })
}

# ============================================================================
# Run
# ============================================================================
shinyApp(ui = ui, server = server)
