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
      actionButton("run_job", "Run pipeline", class = "btn-primary")
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

  observeEvent(input$run_job, {

    req(input$outdir)

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