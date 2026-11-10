library(shiny)
library(shinyWidgets)
library(shinyBS)
library(later)

ui <- fluidPage(
  uiOutput("status_overlay"),

  titlePanel(
    title = div(
      h1("Run GWAS Analysis with Regenie"),
      h4("A graphical interface to configure and launch nf-gwas workflows using Regenie."),
      tags$p(
        HTML(
          '📖 New to nf-gwas? Read the <a href="https://genepi.github.io/nf-gwas/" target="_blank">official tutorial and documentation</a>.'
        ),
        style = "font-size:12px; color:#6c757d; margin-top:0px;"
      )
    ),
    windowTitle = "GWASphere"
  ),
  sidebarLayout(
    sidebarPanel(
      width = 9,

      tabsetPanel(
        id = "config_tabs",

        tabPanel(
          "General",
          br(),

          fluidRow(
            column(
              6,
              radioButtons(
                "test_type",
                HTML('Test type <span title="Choose whether to perform association testing at the variant level (single-variant) or at the gene level using variant aggregation methods.">ⓘ</span>'),
                choices = c("Single-variant" = "sv", "Gene-based" = "gb"),
                selected = "sv",
                inline = TRUE
              )
            ),
            column(
              6,
              textInput("project", "Project name", value = "my_project")
            )
          ),

          fluidRow(
            column(
              6,
              textInput("outdir", HTML("Output directory <span title='Directory where the configuration file, logs and results will be saved.'>ⓘ</span>"), placeholder = "/path/to/regenie_output")
            ),
            column(
              6,
              textInput("partition", "Partition", value = "computational")
            )
          ),

          fluidRow(
            column(
              6,
              selectInput(
                "association_build",
                "Association build",
                choices = c("hg38", "hg19"),
                selected = "hg38"
              )
            )
          )
        ),

        tabPanel(
          "Genotypes",
          br(),

          textInput(
            "genotypes_prediction",
            "Array genotypes (single merged PLINK prefix)",
            value = "",
            placeholder = "/path/to/array_genotypes.{bim,bed,fam}"
          ),

          textInput(
            "genotypes_association",
            "Imputed genotypes (VCF/BGEN path or pattern)",
            value = "",
            placeholder = "/path/to/imputed_chr*.vcf.gz"
          ),

          fluidRow(
            column(
              6,
              selectInput(
                "genotypes_association_format",
                "Imputed genotype format",
                choices = c("vcf", "bgen"),
                selected = "vcf"
              )
            ),
            column(
              6,
            textInput("rsids_filename", "rsID file for annotation",
              placeholder = "/path/to/rsids-v154-hg38.index.gz"),
            div(
              actionLink("rsid_info", label = NULL,
                icon = icon("circle-info"),
                style = "font-size:12px; color:#3498db;"),
              span(style="font-size:12px; color:#888;", " How to download the rsID file")
            ),
            bsModal("rsid_modal", "rsID file — Download instructions",
              "rsid_info", size = "large",
              p("For annotation, a rsid file is required. If not set, this will be downloaded each time. You can also prepare the file once and set it in your pipeline. Substitute rsbuild with hg19 or hg38."),
tags$pre(
"wget https://resources.pheweb.org/rsids-v154-${rsbuild}.tsv.gz -O rsids-v154-${rsbuild}.tsv.gz\necho -e \"CHROM\tPOS\tRSID\tREF\tALT\" | bgzip -c > rsids-v154-${rsbuild}.index.gz\nzcat rsids-v154-${rsbuild}.tsv.gz | bgzip -c >> rsids-v154-${rsbuild}.index.gz\ntabix -s1 -b2 -e2 -S1 rsids-v154-${rsbuild}.index.gz"
)
            )
            )
          ),

          h4("Genotype QC"),
          fluidRow(
            column(
              6,
              numericInput("qc_maf", "Minor allele frequency (MAF)", value = 0.01,min = 0, max = 0.5, step = 0.001)
            ),
            column(
              6,
              numericInput("qc_mac", "Minor allele count (MAC)", value = 100, min = 1, step = 1)
            )
          )
        ),

        tabPanel(
          "Phenotypes & Covariates",
          br(),
          div(style="margin-bottom:12px;",
            actionLink("show_format_example", "See format examples",
                       icon = icon("circle-question"),
                       style = "font-size:13px;"),
            span(style="font-size:12px; color:#888; margin-left:6px;",
                 "— view expected format for phenotype and covariate files")
          ),
          h4("Phenotypes"),
          radioButtons(
            "phenotypes_binary_trait",
            "Binary trait?",
            choices = c("False" = "false", "True" = "true"),
            selected = "false",
            inline = TRUE
          ),

          radioButtons(
            "pheno_input_type",
            HTML("Phenotype input <span title= 'Tab-delimited text file containing sample identifiers (FID, IID) and one or more phenotype columns.'>ⓘ</span>"),
            choices = c("Upload file" = "upload", "Use existing path" = "path"),
            selected = "upload",
            inline = TRUE
          ),

          conditionalPanel(
            condition = "input.pheno_input_type == 'upload'",
            fileInput("phenotypes_file", "Upload phenotypes file")
          ),

          conditionalPanel(
            condition = "input.pheno_input_type == 'path'",
            textInput("phenotypes_filename", "Path to phenotypes file", value = "/path/to/phenotype.txt")
          ),

          radioButtons(
            "phenotypes_delete_missings",
            HTML("Phenotype delete missing <span title= 'Removing samples with missing data at any of the phenotypes.'>ⓘ</span>"),
            choices = c("False" = "false", "True" = "true"),
            selected = "false",
            inline = TRUE
          ),
          uiOutput("phenotypes_columns_ui"),

          tags$hr(),

          h4("Covariates"),
          radioButtons(
            "cova_input_type",
            HTML("Covariate input <span title='Tab-delimited text file containing sample identifiers (FID, IID) and one or more covariate columns. Choose &quot;None&quot; if you do not want to use any covariate.'>ⓘ</span>"),
            choices = c("Upload file" = "upload", "Use existing path" = "path", "None" = "none"),
            selected = "none",
            inline = TRUE
          ),

          conditionalPanel(
            condition = "input.cova_input_type == 'upload'",
            fileInput("covariates_file", "Upload covariates file")
          ),

          conditionalPanel(
            condition = "input.cova_input_type == 'path'",
            textInput("covariates_filename", "Path to covariates file", value = "/path/to/covariates.txt")
          ),

          conditionalPanel(
            condition = "input.cova_input_type == 'none'",
            helpText("No covariates will be included in the analysis.")
          ),

          uiOutput("covariates_columns_ui"),
          uiOutput("covariates_cat_columns_ui"),

          tags$hr(),

          radioButtons(
            "phenotypes_apply_rint",
             HTML("Apply RINT to quantitative phenotypes <span title='Apply rank-based inverse normal transformation to quantitative traits before association testing. Use only for continuous phenotypes, not binary traits.'>ⓘ</span>"),
            choices = c("False" = "false", "True" = "true"),
            selected = "false",
            inline = TRUE
          )
        ),

        tabPanel(
          "Single-variant",
          br(),

          conditionalPanel(
            condition = "input.test_type == 'sv'",

            fluidRow(
              column(
                6,
                selectInput(
                  "regenie_test",
                  "Test model",
                  choices = c("additive", "recessive", "dominant"),
                  selected = "additive"
                )
              ),
              column(
                6,
                numericInput(
                  "annotation_min_log10p",
                  HTML("Annotation threshold (-log10 p) <span title='Annotation limit for interactive and static report.'>ⓘ</span>"),
                  value = 5,min = 0
                )
              )
            ),

            radioButtons(
              "prune_enabled",
              "Enable pruning step",
              choices = c("False" = "false", "True" = "true"),
              selected = "true",
              inline = TRUE
            ),

          conditionalPanel(
            condition = "input.prune_enabled == 'true'",
            fluidRow(
              column(
                3,
                numericInput("prune_maf", "MAF filter", value = 0.01, min = 0, max = 0.5, step = 0.001, width = "100%"),
                tags$p("Minimum allele frequency for variants used in pruning.",
                      style = "font-size:11px; color:#6c757d;")
              ),
              column(
                3,
                numericInput("prune_window_kbsize", "Window size", value = 1000, min = 1, width = "100%"),
                tags$p("LD pruning window size in kilobases.",
                      style = "font-size:11px; color:#6c757d;")
              ),
              column(
                3,
                numericInput("prune_step_size", "Step size", value = 100, min = 1, step = 1, width = "100%"),
                tags$p("Number of variants to shift the pruning window at each step.",
                      style = "font-size:11px; color:#6c757d;")
              ),
              column(
                3,
                numericInput("prune_r2_threshold", "R2 threshold", value = 0.9, min = 0, max = 1, step = 0.01, width = "100%"),
                tags$p("Maximum pairwise LD r² allowed between retained variants.",
                      style = "font-size:11px; color:#6c757d;")
              )
            )
          )
          ),

          conditionalPanel(
            condition = "input.test_type != 'sv'",
            helpText("This section is only used for single-variant analyses.")
          )
        ),

        tabPanel(
          "Gene-based",
          br(),

          conditionalPanel(
            condition = "input.test_type == 'gb'",

            radioButtons(
              "regenie_run_gene_based_tests",
              "Activate gene-based testing",
              choices = c("False" = "false", "True" = "true"),
              selected = "true",
              inline = TRUE
            ),

            textInput(
              "regenie_gene_anno",
              HTML("Variant annotation file <span title='File with variant annotations for each set.'>ⓘ</span>"),
              value = "",
              placeholder = "/path/to/gene_annotation.annot"
            ),

            textInput(
              "regenie_gene_setlist",
              HTML("Set list file <span title='File listing variant sets.'>ⓘ</span>"),
              value = "",
              placeholder = "/path/to/gene_setlist.setlist"
            ),

            textInput(
              "regenie_gene_masks",
              HTML("Mask definition file <span title='File with mask definitions using the annotations defined in regenie_gene_anno.'>ⓘ</span>"),
              value = "",
              placeholder = "/path/to/gene_masks.masks"
            ),

            textInput(
              "regenie_gene_test",
              "Gene-based tests",
              value = "skato,acato"
            )
          ),

          conditionalPanel(
            condition = "input.test_type != 'gb'",
            helpText("This section is only used for gene-based analyses.")
          )
        )
      ),

      br(),
      div(style="margin-bottom:8px;",
        actionButton("validate_inputs", "Validate inputs", class = "btn-default",
                     icon = icon("circle-check")),
        span(style="font-size:12px; color:#888; margin-left:8px;",
             "Check that sample IDs match across genotype, phenotype and covariate files before launching.")
      ),
      uiOutput("validation_summary_ui"),
      br(),
      uiOutput("run_pipeline_ui")
    ),

    mainPanel(width = 3)
  )
)

server <- function(input, output, session) {

    # CHANGE DEFAULT INPUT BASED ON INPUT DATASET


    # read columns of the phenotypes file
  output$phenotypes_columns_ui <- renderUI({
    path <- NULL

    if (input$pheno_input_type == "upload") {
      req(input$phenotypes_file)
      path <- input$phenotypes_file$datapath
    } else {
      req(input$phenotypes_filename)
      path <- input$phenotypes_filename
    }

    if (!file.exists(path)) {
      return(helpText("⚠️ File not found."))
    }

    cols <- tryCatch({
      all_cols <- colnames(read.table(path, header = TRUE, nrows = 1, sep = "", check.names = FALSE))
      exclude <- c("fid", "iid", "aid", "aid_num")
      all_cols[!tolower(all_cols) %in% exclude]
    }, error = function(e) NULL)

    if (is.null(cols) || length(cols) == 0) {
      return(helpText("⚠️ No columns found."))
    }

    pickerInput(
      inputId = "phenotypes_columns",
      label = "List of phenotypes",
      choices = cols,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE
      )
    )
  })

  # read columns of the covariates file (skipped entirely when no covariate file is used)
  output$covariates_columns_ui <- renderUI({

    if (is.null(input$cova_input_type) || input$cova_input_type == "none") {
      return(NULL)
    }

    path <- NULL

    if (input$cova_input_type == "upload") {
      req(input$covariates_file)
      path <- input$covariates_file$datapath
    } else {
      req(input$covariates_filename)
      path <- input$covariates_filename
    }

    if (!file.exists(path)) {
      return(helpText("⚠️ File not found."))
    }

    cols <- tryCatch({
      all_cols <- colnames(read.table(path, header = TRUE, nrows = 1, sep = "", check.names = FALSE))
      exclude <- c("fid", "iid", "aid", "aid_num")
      all_cols[!tolower(all_cols) %in% exclude]
    }, error = function(e) NULL)

    if (is.null(cols) || length(cols) == 0) {
      return(helpText("⚠️ No columns found."))
    }

    selectInput(
      "covariates_columns",
      "List of covariates",
      choices = cols,
      multiple = TRUE
    )
  })

    output$covariates_cat_columns_ui <- renderUI({

      if (is.null(input$cova_input_type) || input$cova_input_type == "none") {
        return(NULL)
      }

      path <- NULL

      if (input$cova_input_type == "upload") {
        req(input$covariates_file)
        path <- input$covariates_file$datapath
      } else {
        req(input$covariates_filename)
        path <- input$covariates_filename
      }

      # if (!file.exists(path)) {
      #   return(helpText("⚠️ File not found."))
      # }

      all_cols <- tryCatch({
        colnames(read.table(path, header = TRUE, nrows = 1, sep = "", check.names = FALSE))
      }, error = function(e) NULL)

      # if (is.null(all_cols)) {
      #   return(helpText("⚠️ No columns found."))
      # }

      exclude <- c("fid", "iid", "aid", "aid_num")
      all_cols <- all_cols[!tolower(all_cols) %in% exclude]

      # differenza: categoriche = tutto meno le covariate selezionate
      cat_cols <- setdiff(all_cols, input$covariates_columns)

      selectInput(
        "covariates_cat_columns",
        "Categorical covariates",
        choices = cat_cols,
        multiple = TRUE
      )
    })

  #status_text <- reactiveVal("No job submitted yet.")
  status_text <- reactiveVal("")

  # ── Pre-flight validation ──────────────────────────────────────────────────
  format_preview <- function(df, label) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    rows_txt <- apply(df, 1, function(r) paste(r, collapse = "  "))
    div(
      strong(label),
      tags$pre(style = "font-size:11px; background:#fff3cd; padding:6px; border-radius:4px; margin:4px 0 0 0;",
               paste(rows_txt, collapse = "
"))
    )
  }
  run_validation <- function(pheno_path, cova_path, geno_prefix,
                             selected_pheno_cols = NULL, selected_cova_cols = NULL) {
    results <- list()

    # Read .fam file — handle {bim,bed,fam} glob pattern
    fam_path <- sub("[{].*$", "", geno_prefix)
    fam_path <- sub("[.]$", "", fam_path)
    fam_path <- paste0(fam_path, ".fam")

    fam_ok <- file.exists(fam_path)
    if (fam_ok) {
      fam <- tryCatch(read.table(fam_path, header = FALSE, stringsAsFactors = FALSE, colClasses = "character"),
                      error = function(e) NULL)
      fam_ids <- if (!is.null(fam)) paste(fam[[1]], fam[[2]], sep = "_") else character(0)
      results$n_geno <- length(fam_ids)
    } else {
      fam_ids <- character(0)
      results$n_geno <- NA
      results$fam_warning <- paste0("⚠️ .fam file not found: ", fam_path)
    }

    # Read phenotype file
    pheno_ok <- file.exists(pheno_path)
    if (pheno_ok) {
      pheno <- tryCatch(read.table(pheno_path, header = TRUE, stringsAsFactors = FALSE, colClasses = "character"),
                        error = function(e) NULL)
      if (!is.null(pheno) && all(c("FID","IID") %in% names(pheno))) {
        pheno_ids <- paste(pheno$FID, pheno$IID, sep = "_")
        results$n_pheno_total <- length(pheno_ids)
        if (length(fam_ids) > 0) {
          matched_pheno <- sum(pheno_ids %in% fam_ids)
          results$n_pheno_matched <- matched_pheno
          results$n_pheno_excluded <- length(pheno_ids) - matched_pheno
          # If no match or very low match, store preview rows
          if (matched_pheno == 0 || matched_pheno / length(pheno_ids) < 0.01) {
            results$low_match_pheno <- TRUE
            results$pheno_preview  <- head(pheno[, c("FID","IID"), drop=FALSE], 5)
            results$fam_preview    <- if (!is.null(fam)) head(fam[, 1:2], 5) else NULL
          }
        }
      } else {
        results$pheno_warning <- "⚠️ Phenotype file missing FID/IID columns"
        pheno_ids <- character(0)
      }
    } else {
      results$pheno_warning <- paste0("❌ Phenotype file not found: ", pheno_path)
      pheno_ids <- character(0)
    }

    # Read covariate file (optional)
    if (!is.null(cova_path) && nzchar(cova_path) && file.exists(cova_path)) {
      cova <- tryCatch(read.table(cova_path, header = TRUE, stringsAsFactors = FALSE, colClasses = "character"),
                       error = function(e) NULL)
      if (!is.null(cova) && all(c("FID","IID") %in% names(cova))) {
        cova_ids <- paste(cova$FID, cova$IID, sep = "_")
        results$n_cova_total <- length(cova_ids)

        # Step 1: covariate vs phenotype
        if (length(pheno_ids) > 0) {
          matched_cova_pheno <- sum(cova_ids %in% pheno_ids)
          results$n_cova_matched_pheno  <- matched_cova_pheno
          results$n_cova_excluded_pheno <- length(cova_ids) - matched_cova_pheno
          cova_pheno_ids <- cova_ids[cova_ids %in% pheno_ids]
          if (matched_cova_pheno == 0 || matched_cova_pheno / length(cova_ids) < 0.01) {
            results$low_match_cova_pheno <- TRUE
            results$cova_preview  <- head(cova[, c("FID","IID"), drop=FALSE], 5)
            results$pheno_preview_for_cova <- if (!is.null(pheno)) head(pheno[, c("FID","IID"), drop=FALSE], 5) else NULL
          }
        } else {
          cova_pheno_ids <- cova_ids
        }

        # Step 2: of those, vs genotypes
        if (length(fam_ids) > 0) {
          matched_cova_geno <- sum(cova_pheno_ids %in% fam_ids)
          results$n_cova_matched_geno  <- matched_cova_geno
          results$n_cova_excluded_geno <- length(cova_pheno_ids) - matched_cova_geno
          if (matched_cova_geno == 0 || (length(cova_pheno_ids) > 0 && matched_cova_geno / length(cova_pheno_ids) < 0.01)) {
            results$low_match_cova <- TRUE
            results$cova_preview   <- head(cova[, c("FID","IID"), drop=FALSE], 5)
            results$fam_preview    <- if (!is.null(fam)) head(fam[, 1:2], 5) else NULL
          }
        }

        # Final N = intersection of all three
        if (length(pheno_ids) > 0 && length(fam_ids) > 0 && exists("cova_ids") && length(cova_ids) > 0) {
          # Filter covariates to selected columns only, exclude rows with NA in selected cols
          cova_cols_to_check <- if (!is.null(selected_cova_cols) && length(selected_cova_cols) > 0)
            intersect(selected_cova_cols, names(cova))
          else
            setdiff(names(cova), c("FID", "IID"))
          if (length(cova_cols_to_check) > 0) {
            cova_complete <- apply(cova[, cova_cols_to_check, drop = FALSE], 1, function(r) {
              !any(is.na(r) | as.character(r) == "NA")
            })
            cova_ids_complete <- cova_ids[cova_complete]
            # NA per covariate among matched individuals (pheno ∩ geno)
            cova_matched_mask <- paste(cova$FID, cova$IID, sep = "_") %in%
              pheno_ids[pheno_ids %in% fam_ids]
            cova_matched <- cova[cova_matched_mask, cova_cols_to_check, drop = FALSE]
            results$na_per_cova <- sapply(cova_cols_to_check, function(col) {
              v <- cova_matched[[col]]
              sum(is.na(v) | as.character(v) == "NA")
            })
          } else {
            cova_ids_complete <- cova_ids
          }
          final_mask <- pheno_ids %in% fam_ids & pheno_ids %in% cova_ids_complete
          results$n_final <- sum(final_mask)
          # NA counts per trait within Final N (only selected traits)
          trait_cols <- intersect(selected_pheno_cols, names(pheno))
          if (length(trait_cols) == 0) {
            results$pheno_col_warning <- "⚠️ No phenotype columns selected."
          }
          if (length(trait_cols) > 0) {
            pheno_final <- pheno[final_mask, , drop = FALSE]
            results$na_per_trait <- sapply(trait_cols, function(col) {
              v <- pheno_final[[col]]
              sum(is.na(v) | as.character(v) == "NA")
            })
          }
        }
      } else {
        results$cova_warning <- "⚠️ Covariate file missing FID/IID columns"
      }
    } else if (!is.null(cova_path) && nzchar(cova_path)) {
      results$cova_warning <- paste0("⚠️ Covariate file not found: ", cova_path)
    }

    # Final N without covariates
    if (is.null(results$n_final) && length(pheno_ids) > 0 && length(fam_ids) > 0) {
      final_mask2 <- pheno_ids %in% fam_ids
      results$n_final <- sum(final_mask2)
      trait_cols2 <- intersect(selected_pheno_cols, names(pheno))
      if (length(trait_cols2) == 0) {
        results$pheno_col_warning <- "⚠️ No phenotype columns selected."
      }
      if (length(trait_cols2) > 0) {
        pheno_final2 <- pheno[final_mask2, , drop = FALSE]
        results$na_per_trait <- sapply(trait_cols2, function(col) {
          v <- pheno_final2[[col]]
          sum(is.na(v) | as.character(v) == "NA")
        })
      }
    }

    results
  }

  # ── Format example modal ───────────────────────────────────────────────────
  observeEvent(input$show_format_example, {
    showModal(modalDialog(
      title = "Input file format examples",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      tagList(
        p(style="font-size:13px; color:#555;",
          "Files must be tab- or space-separated with a header row. ",
          "The first two columns must be ", strong("FID"), " and ", strong("IID"),
          " (can be the same value). Missing values must be coded as ", strong("NA"), "."),
        hr(),
        h5("Phenotype file — quantitative traits"),
        p(style="font-size:12px; color:#888;",
          "One column per quantitative trait. Any numeric value. Missing values as NA."),
        tags$table(
          class = "table table-bordered table-condensed",
          style = "font-size:12px; width:auto;",
          tags$thead(tags$tr(
            tags$th("FID"), tags$th("IID"),
            tags$th("bmi"), tags$th("height"), tags$th("hdl")
          )),
          tags$tbody(
            tags$tr(tags$td("ID001"), tags$td("ID001"), tags$td("24.5"), tags$td("172.3"), tags$td("1.2")),
            tags$tr(tags$td("ID002"), tags$td("ID002"), tags$td("NA"),   tags$td("165.1"), tags$td("0.9")),
            tags$tr(tags$td("ID003"), tags$td("ID003"), tags$td("27.1"), tags$td("NA"),    tags$td("1.5")),
            tags$tr(tags$td("ID004"), tags$td("ID004"), tags$td("31.2"), tags$td("180.0"), tags$td("NA")),
            tags$tr(tags$td("ID005"), tags$td("ID005"), tags$td("22.8"), tags$td("158.7"), tags$td("1.1"))
          )
        ),
        br(),
        h5("Phenotype file — binary traits"),
        p(style="font-size:12px; color:#888;",
          "One column per binary trait. 0 = control, 1 = case, NA = missing. Keep binary traits in a separate file from quantitative traits."),
        tags$table(
          class = "table table-bordered table-condensed",
          style = "font-size:12px; width:auto;",
          tags$thead(tags$tr(
            tags$th("FID"), tags$th("IID"),
            tags$th("t2d"), tags$th("hypertension")
          )),
          tags$tbody(
            tags$tr(tags$td("ID001"), tags$td("ID001"), tags$td("0"), tags$td("1")),
            tags$tr(tags$td("ID002"), tags$td("ID002"), tags$td("1"), tags$td("0")),
            tags$tr(tags$td("ID003"), tags$td("ID003"), tags$td("NA"), tags$td("0")),
            tags$tr(tags$td("ID004"), tags$td("ID004"), tags$td("0"), tags$td("NA")),
            tags$tr(tags$td("ID005"), tags$td("ID005"), tags$td("1"), tags$td("1"))
          )
        ),
        hr(),
        h5("Covariate file"),
        p(style="font-size:12px; color:#888;",
          "One column per covariate. Sex: typically 1 = male, 2 = female. Missing values as NA."),
        tags$table(
          class = "table table-bordered table-condensed",
          style = "font-size:12px; width:auto;",
          tags$thead(tags$tr(
            tags$th("FID"), tags$th("IID"),
            tags$th("age"), tags$th("sex"), tags$th("PC1"), tags$th("PC2")
          )),
          tags$tbody(
            tags$tr(tags$td("ID001"), tags$td("ID001"), tags$td("45"), tags$td("1"), tags$td("0.0023"),  tags$td("-0.0011")),
            tags$tr(tags$td("ID002"), tags$td("ID002"), tags$td("52"), tags$td("2"), tags$td("-0.0014"), tags$td("0.0032")),
            tags$tr(tags$td("ID003"), tags$td("ID003"), tags$td("38"), tags$td("1"), tags$td("0.0041"),  tags$td("0.0017")),
            tags$tr(tags$td("ID004"), tags$td("ID004"), tags$td("NA"), tags$td("2"), tags$td("0.0019"),  tags$td("-0.0008")),
            tags$tr(tags$td("ID005"), tags$td("ID005"), tags$td("61"), tags$td("1"), tags$td("-0.0031"), tags$td("0.0025"))
          )
        ),
        hr(),
        div(style="background:#fff3cd; border:1px solid #ffc107; border-radius:6px; padding:10px 14px; font-size:12px; color:#856404;",
          icon("triangle-exclamation"), strong(" Important:"),
          tags$ul(style="margin:4px 0 0 0;",
            tags$li("Missing values must be coded as ", strong("NA"), " — not 0, -9, or empty cells."),
            tags$li("For binary traits, use ", strong("0"), " = control and ", strong("1"), " = case."),
            tags$li("FID and IID must match exactly with the genotype .fam file.")
          )
        )
      )
    ))
  })

  observeEvent(input$validate_inputs, {
    pheno_path <- if (input$pheno_input_type == "upload" && !is.null(input$phenotypes_file))
      input$phenotypes_file$datapath else input$phenotypes_filename
    cova_path <- if (input$cova_input_type == "upload" && !is.null(input$covariates_file))
      input$covariates_file$datapath
    else if (input$cova_input_type == "path") input$covariates_filename
    else NULL

    geno_prefix <- input$genotypes_prediction
    selected_pheno <- input$phenotypes_columns
    selected_cova  <- unique(c(input$covariates_columns, input$covariates_cat_columns))
    res <- run_validation(pheno_path, cova_path, geno_prefix, selected_pheno, selected_cova)
    validation_results(res)
  })

  validation_results <- reactiveVal(NULL)

  output$validation_summary_ui <- renderUI({
    res <- validation_results()
    if (is.null(res)) return(NULL)

    rows <- list()

    if (!is.null(res$fam_warning)) {
      rows <- c(rows, list(div(style="color:#856404;", res$fam_warning)))
    } else {
      rows <- c(rows, list(div(
        icon("dna"), strong(" Genotypes: "),
        paste0(format(res$n_geno, big.mark=","), " individuals in .fam")
      )))
    }

    if (!is.null(res$pheno_col_warning)) {
      rows <- c(rows, list(div(style="color:#856404;", res$pheno_col_warning)))
    }
    if (!is.null(res$pheno_warning)) {
      rows <- c(rows, list(div(style="color:#721c24;", res$pheno_warning)))
    } else if (!is.null(res$n_pheno_matched)) {
      col <- if (isTRUE(res$n_pheno_excluded == 0)) "#155724" else "#856404"
      rows <- c(rows, list(div(style=paste0("color:",col,";"),
        icon("table"), strong(" Phenotype: "),
        paste0(format(res$n_pheno_matched, big.mark=","), " / ",
               format(res$n_pheno_total, big.mark=","),
               " individuals match genotypes"),
        if (isTRUE(res$n_pheno_excluded > 0))
          span(style="color:#856404;",
               paste0(" (", format(res$n_pheno_excluded, big.mark=","), " will be excluded)"))
      )))
      # Show preview if match is 0 or very low
      if (isTRUE(res$low_match_pheno)) {
        preview_block <- div(
          style = "margin-top:8px; font-size:12px; color:#721c24;",
          strong("⚠️ Very few or no ID matches detected. Please check ID format:"),
          br(),
          div(style="display:flex; gap:20px; margin-top:6px;",
            format_preview(res$pheno_preview, "Phenotype file (FID, IID):"),
            format_preview(res$fam_preview,   ".fam file (col1=FID, col2=IID):")
          )
        )
        rows <- c(rows, list(preview_block))
      }
    }

    if (!is.null(res$cova_warning)) {
      rows <- c(rows, list(div(style="color:#856404;", res$cova_warning)))
    } else if (!is.null(res$n_cova_matched_pheno) || !is.null(res$n_cova_matched_geno)) {
      col <- if (isTRUE(identical(res$n_cova_excluded_geno, 0L) || identical(res$n_cova_excluded_geno, 0)) &&
                 isTRUE(identical(res$n_cova_excluded_pheno, 0L) || identical(res$n_cova_excluded_pheno, 0))) "#155724" else "#856404"
      rows <- c(rows, list(div(style=paste0("color:",col,";"),
        icon("sliders"), strong(" Covariates: "),
        if (!is.null(res$n_cova_matched_pheno))
          tagList(
            paste0(format(res$n_cova_matched_pheno, big.mark=","), " / ",
                   format(res$n_cova_total, big.mark=","),
                   " match phenotype"),
            if (isTRUE(res$n_cova_excluded_pheno > 0))
              span(style="color:#856404;",
                   paste0(" (", format(res$n_cova_excluded_pheno, big.mark=","), " excluded)")),
            if (!is.null(res$n_cova_matched_geno))
              tagList(
                br(),
                span(style="padding-left:16px;",
                  paste0("→ of those: ", format(res$n_cova_matched_geno, big.mark=","), " / ",
                         format(res$n_cova_matched_pheno, big.mark=","),
                         " match genotypes"),
                  if (isTRUE(res$n_cova_excluded_geno > 0))
                    span(style="color:#856404;",
                         paste0(" (", format(res$n_cova_excluded_geno, big.mark=","), " excluded)"))
                )
              )
          )
      )))
      if (isTRUE(res$low_match_cova_pheno)) {
        preview_block <- div(
          style = "margin-top:8px; font-size:12px; color:#721c24;",
          strong("⚠️ Very few or no ID matches detected between covariates and phenotype. Please check ID format:"),
          br(),
          div(style="display:flex; gap:20px; margin-top:6px;",
            format_preview(res$cova_preview, "Covariate file (FID, IID):"),
            format_preview(res$pheno_preview_for_cova, "Phenotype file (FID, IID):")
          )
        )
        rows <- c(rows, list(preview_block))
      } else if (isTRUE(res$low_match_cova)) {
        preview_block <- div(
          style = "margin-top:8px; font-size:12px; color:#721c24;",
          strong("⚠️ Very few or no ID matches detected between covariates and genotypes. Please check ID format:"),
          br(),
          div(style="display:flex; gap:20px; margin-top:6px;",
            format_preview(res$cova_preview, "Covariate file (FID, IID):"),
            format_preview(res$fam_preview,  ".fam file (col1=FID, col2=IID):")
          )
        )
        rows <- c(rows, list(preview_block))
      }
    }

    if (!is.null(res$na_per_cova) && length(res$na_per_cova) > 0) {
      na_cova_lines <- sapply(names(res$na_per_cova), function(cv) {
        paste0(cv, ": ", res$na_per_cova[[cv]], " NA")
      })
      rows <- c(rows, list(div(style="font-size:12px; color:#555; padding-left:20px;",
        strong("Missing values in selected covariates:"), br(),
        HTML(paste(na_cova_lines, collapse = "<br>"))
      )))
    }

    if (!is.null(res$n_final)) {
      final_block <- list(
        tags$hr(style="margin:8px 0;"),
        div(style="font-weight:600; color:#0c5460;",
          icon("users"), " Final N for analysis: ",
          format(res$n_final, big.mark=","), " individuals"
        )
      )
      if (!is.null(res$na_per_trait) && length(res$na_per_trait) > 0) {
        na_lines <- sapply(names(res$na_per_trait), function(tr) {
          n_na <- res$na_per_trait[[tr]]
          n_analysed <- res$n_final - n_na
          paste0(tr, ": ", n_na, " NA (", format(n_analysed, big.mark=","), " analysed)")
        })
        final_block <- c(final_block, list(
          div(style="font-size:12px; color:#555; margin-top:4px; padding-left:20px;",
            strong("Missing values per trait:"), br(),
            HTML(paste(na_lines, collapse = "<br>"))
          )
        ))
      }
      rows <- c(rows, list(tagList(final_block)))
    }

    div(
      style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:6px;
               padding:12px 16px; margin-top:12px; font-size:13px; line-height:2;",
      tagList(rows)
    )
  })

  output$run_pipeline_ui <- renderUI({
    res <- validation_results()
    all_ok <- !is.null(res) &&
      is.null(res$fam_warning) &&
      is.null(res$pheno_warning) &&
      isTRUE(res$n_pheno_matched > 0) &&
      (is.null(res$n_cova_matched_pheno) || isTRUE(res$n_cova_matched_pheno > 0)) &&
      (is.null(res$n_cova_matched_geno)  || isTRUE(res$n_cova_matched_geno > 0))
    if (all_ok) {
      tagList(
        div(style="background:#d4edda; border:1px solid #c3e6cb; border-radius:6px;
                   padding:10px 16px; margin-bottom:10px; font-size:13px; color:#155724;",
          icon("circle-check"), strong(" Input files validated — ready to run."),
          " 🚀 Ready to go? Hit the button below to launch the pipeline!"
        ),
        actionButton("run_job", "Run pipeline", class = "btn-success", icon = icon("play"))
      )
    } else {
      tagList(
        div(style="font-size:12px; color:#888; margin-bottom:6px;",
          "Validate inputs above before launching"
        ),
        actionButton("run_job", "Run pipeline", class = "btn-primary", icon = icon("play"))
      )
    }
  })

  observeEvent(input$run_job, {

    req(input$outdir)

    # ── Check Java and Nextflow before submitting ──────────────────────────
    java_bin <- if (nzchar(Sys.getenv("NXF_JAVA_HOME")))
      file.path(Sys.getenv("NXF_JAVA_HOME"), "bin", "java") else "java"
    nf_bin <- if (nzchar(Sys.getenv("NEXTFLOW_BIN")))
      Sys.getenv("NEXTFLOW_BIN") else "nextflow"

    java_check <- tryCatch(
      system2(java_bin, args = "-version", stdout = TRUE, stderr = TRUE),
      error = function(e) NULL
    )
    nf_check <- tryCatch(
      system2(nf_bin, args = "-version", stdout = TRUE, stderr = TRUE),
      error = function(e) NULL
    )

    errors <- c()
    if (is.null(java_check) || length(java_check) == 0)
      errors <- c(errors, paste0("❌ Java not found at: ", java_bin, ". Set NXF_JAVA_HOME in run.sh."))
    if (is.null(nf_check) || length(nf_check) == 0)
      errors <- c(errors, paste0("❌ Nextflow not found at: ", nf_bin, ". Set NEXTFLOW_BIN in run.sh."))

    if (length(errors) > 0) {
      showModal(modalDialog(
        title = "⚠️ Pre-flight check failed",
        size = "m",
        easyClose = TRUE,
        footer = modalButton("Close"),
        div(style="font-size:13px;",
          p("The following issues were detected before submitting the job:"),
          tags$ul(lapply(errors, tags$li)),
          p(style="color:#888; margin-top:12px;",
            "Edit run.sh to set the correct paths, then restart the app.")
        )
      ))
      return()
    }

    output_dir <- normalizePath(input$outdir, winslash = "/", mustWork = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # phenotype file
    if (input$pheno_input_type == "upload") {
      req(input$phenotypes_file)

      uploaded_path <- input$phenotypes_file$datapath
      uploaded_name <- input$phenotypes_file$name

      pheno_final_path <- file.path(output_dir, uploaded_name)
      file.copy(uploaded_path, pheno_final_path, overwrite = TRUE)

    } else {
      req(input$phenotypes_filename)
      pheno_final_path <- input$phenotypes_filename
    }

    # covariate file (OPTIONAL: if cova_input_type == "none", no covariates are used
    # and cova_final_path stays NULL, so no covariate params are written to the config)
    cova_final_path <- NULL

    if (input$cova_input_type == "upload") {
      req(input$covariates_file)

      uploaded_path <- input$covariates_file$datapath
      uploaded_name <- input$covariates_file$name

      cova_final_path <- file.path(output_dir, uploaded_name)
      file.copy(uploaded_path, cova_final_path, overwrite = TRUE)

    } else if (input$cova_input_type == "path") {
      req(input$covariates_filename)
      cova_final_path <- input$covariates_filename
    }
    # if input$cova_input_type == "none" -> cova_final_path remains NULL

    config_file <- file.path(output_dir, paste0(input$project, ".config"))

    config_lines <- c(
      "/* General genotype configuration */",
      paste0("params.project=\"", input$project, "\""),
      paste0("params.outdir=\"", input$outdir, "\""),
      paste0("params.association_build=\"", input$association_build, "\""),
      paste0("params.genotypes_prediction=\"", input$genotypes_prediction, "\""),
      paste0("params.genotypes_association=\"", input$genotypes_association, "\""),
      paste0("params.genotypes_association_format=\"", input$genotypes_association_format, "\""),
      paste0("params.rsids_filename=\"", input$rsids_filename, "\"")
    )

    if (input$test_type == "sv") {
      config_lines <- c(
        config_lines,
        paste0("params.regenie_test=\"", input$regenie_test, "\""),
        paste0("params.annotation_min_log10p=", input$annotation_min_log10p),
        paste0("params.prune_enabled=", input$prune_enabled),
        paste0("params.prune_maf=", input$prune_maf),
        paste0("params.prune_window_kbsize=", input$prune_window_kbsize),
        paste0("params.prune_step_size=", input$prune_step_size),
        paste0("params.prune_r2_threshold=", input$prune_r2_threshold)
      )
    } else {
      config_lines <- c(
        config_lines,
        paste0("params.regenie_run_gene_based_tests=", input$regenie_run_gene_based_tests),
        paste0("params.regenie_gene_anno=\"", input$regenie_gene_anno, "\""),
        paste0("params.regenie_gene_setlist=\"", input$regenie_gene_setlist, "\""),
        paste0("params.regenie_gene_masks=\"", input$regenie_gene_masks, "\"")
      )
    }

    config_lines <- c(
      config_lines,
      "/* Genotype QC parameters */",
      paste0("params.qc_maf=", input$qc_maf),
      paste0("params.qc_mac=", input$qc_mac),
      "/* Phenotype specific parameters */",
      paste0("params.phenotypes_binary_trait=", input$phenotypes_binary_trait),
      paste0("params.phenotypes_delete_missings=", input$phenotypes_delete_missings),
      paste0("params.phenotypes_apply_rint=", input$phenotypes_apply_rint),
      paste0("params.phenotypes_filename=\"", pheno_final_path, "\""),
      paste0("params.phenotypes_columns=\"", paste(input$phenotypes_columns, collapse = ","), "\"")
    )

    # Covariates are optional in nf-gwas (default: empty). Only write them to the
    # config file if the user actually provided a covariate file/path.
    if (!is.null(cova_final_path) && nzchar(cova_final_path)) {
      config_lines <- c(
        config_lines,
        "/* Covariates specific parameters */",
        paste0("params.covariates_filename=\"", cova_final_path, "\""),
        paste0("params.covariates_columns=\"", paste(input$covariates_columns, collapse = ","), "\""),
        paste0("params.covariates_cat_columns=\"", paste(input$covariates_cat_columns, collapse = ","), "\"")
      )
    }

    process_block <- c(
      "process{",
      "    executor='slurm'",
      "    clusterOptions = '--requeue'",
      paste0("    queue='", input$partition, "'"),
      "    // clusterOptions = '--time=48:00:00'",
      "",
      "    withName: 'IMPUTED_TO_PLINK2' {",
      "        memory = { 32.GB * task.attempt }",
      "    }",
      "",
      "    withLabel: 'required_memory_report' {",
      "        memory = 32.GB",
      "    }",
      "",
      "    withName: 'MERGE_RESULTS' {",
      "        memory = 16.GB",
      "    }",
      "",
      "    withName: 'REGENIE_STEP1_RUN|REGENIE_STEP2_RUN' {",
      "        cpus   = 8",
      "        memory = { 16.GB * task.attempt }",
      "    }",
      "",
      "    withName: 'ANNOTATE_RESULTS' {",
      "        memory = { 16.GB * task.attempt }",
      "    }",
      "",
      "    withName: 'REGENIE_STEP2_RUN_GENE_TESTS' {",
      "        memory = 32.GB",
      "    }",
      "}"
    )

    all_config_lines <- c(config_lines, "", process_block)
    writeLines(all_config_lines, config_file)

    job_name <- gsub("[^A-Za-z0-9_-]", "_", input$project)

    res <- tryCatch(
      system2(
        "sbatch",
        args = c(
          "--job-name", job_name,
          "--partition", input$partition,
          "--output", file.path(output_dir, paste0(job_name, "-%j.out")),
          "--error", file.path(output_dir, paste0(job_name, "-%j.err")),
          file.path(REGENIE_APP_DIR, "regenie_pipeline/bin/run_pipeline.sh"),
          config_file
        ),
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) paste("Submission failed:", e$message)
    )

    status_text(paste(
      c(
        paste("Config file saved:", basename(config_file)),
        paste("Job submitted with output directory:", output_dir),
        res
      ),
      collapse = "\n"
    ))

    later::later(function() status_text(""), delay = 10)
  })

  output$status_overlay <- renderUI({
    req(nzchar(status_text()))

    tagList(
      div(
        style = "
          position: fixed;
          top: 0;
          left: 0;
          width: 100%;
          height: 100%;
          background-color: rgba(0,0,0,0.3);
          z-index: 9998;
        "
      ),
      div(
        style = "
          position: fixed;
          top: 50%;
          left: 50%;
          transform: translate(-50%, -50%);
          background-color: #d4edda;
          color: #155724;
          padding: 25px 30px;
          border-radius: 10px;
          border-left: 6px solid #28a745;
          box-shadow: 0px 4px 20px rgba(0,0,0,0.2);
          z-index: 9999;
          text-align: center;
          min-width: 320px;
          white-space: pre-line;
          font-weight: 500;
        ",
        "✅ ",
        status_text()
      )
    )
  })
}

shinyApp(ui, server)