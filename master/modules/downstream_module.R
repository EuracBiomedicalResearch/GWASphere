# =============================================================================
# modules/downstream_module.R
# =============================================================================

downstreamUI <- function(id) {
  ns <- NS(id)

  fluidPage(
    div(
      style = "background:#eafaf1; border-left:4px solid #27ae60;
               padding:8px 14px; margin-bottom:16px; font-size:13px;",
      icon("circle-info"),
      " Configure paths and launch the post-GWAS analysis."
    ),

    fluidRow(

      # --- Left column: configuration ------------------------------------
      column(6,

        # Required inputs
        wellPanel(
          h4(icon("folder-open"), " Input & output directories"),

          textInput(ns("input_dir"),
            HTML('REGENIE results directory <span title="Folder containing .regenie.gz files">ⓘ</span>'),
            value = "", placeholder = "/path/to/regenie_output/results/"),

          textInput(ns("out_dir"),
            HTML('Post-GWAS output directory <span title="Directory where all downstream analysis results will be saved. Will be created if it does not exist">ⓘ</span>'),
            value = "", placeholder = "/path/to/downstream_output/"),

          h4(icon("dna"), " Reference files"),

          textInput(ns("ref_prefix"),
            HTML('PLINK reference file prefix <span title="For per-chromosome files: prefix without number (e.g. /path/to/imputed_chr → finds imputed_chr1.bed, imputed_chr2.bed...). For a single merged file: full prefix (e.g.  /path/to/imputed_all → finds imputed_all.bed). Auto-detected.">ⓘ</span>'),
            placeholder = "/path/to/imputed_chr"),

          h4(icon("filter"), " Clumping & locus definition"),

          fluidRow(
            column(6,
              textInput(ns("p1"),
                HTML('Genome-wide significance threshold (P1) <span title="P-value threshold for independent significant SNPs">ⓘ</span>'),
                value = "5e-8")
            ),
            column(6,
              textInput(ns("ptag"),
                HTML('Candidate SNP P threshold (PTAG) <span title="P-value threshold for candidate/tagged SNPs">ⓘ</span>'),
                value = "0.05")
            )
          ),
          fluidRow(
            column(4,
              numericInput(ns("step1_r2"),
                HTML('Step 1 r² <span title="LD r² threshold for independent significant SNPs">ⓘ</span>'),
                value = 0.6, min = 0.1, max = 1, step = 0.05)
            ),
            column(4,
              numericInput(ns("step2_r2"),
                HTML('Step 2 r² <span title="LD r² threshold for lead SNPs">ⓘ</span>'),
                value = 0.1, min = 0.05, max = 0.5, step = 0.05)
            ),
            column(4,
              numericInput(ns("merge_kb"),
                HTML('Merge distance (kb) <span title="Distance to merge overlapping loci">ⓘ</span>'),
                value = 250, min = 0, max = 1000, step = 50)
            )
          ),
          fluidRow(
            column(6,
              numericInput(ns("maf"),
                HTML('MAF filter <span title="Minimum minor allele frequency">ⓘ</span>'),
                value = 0.001, min = 0, max = 0.1, step = 0.001)
            )
          ),

          h4(icon("sliders"), " SLURM"),

          fluidRow(
            column(6, textInput(ns("partition"), "Partition", value = "computational")),
            column(6, numericInput(ns("mem_gb"), "Memory (GB)", value = 128, min = 8, max = 512, step = 8))
          ),
          fluidRow(
            column(6, numericInput(ns("cpus"), "CPUs", value = 8, min = 1, max = 64, step = 1))
          )
        ),

        # --- Optional: MAGMA ------------------------------------------
        wellPanel(
          style = "border-left: 4px solid #3498db;",
          checkboxInput(ns("run_magma"),
            span(icon("brain"), strong(" Run MAGMA"), " — gene, gene-set and tissue analysis"),
            value = TRUE),

          conditionalPanel(
            condition = paste0("input['", ns("run_magma"), "'] == true"),
            div(style = "margin-top:10px;",
              radioButtons(ns("genome_build"),
                HTML('Genome build <span title="Select the genome build used in your GWAS. This sets the MAGMA gene location file and SnpEff genome automatically.">ⓘ</span>'),
                choices = c("hg38 (GRCh38)" = "hg38", "hg19 (GRCh37)" = "hg19"),
                selected = "hg38", inline = TRUE),
              textInput(ns("gene_loc"),
                "Gene location file",
                value = file.path(Sys.getenv("GWASPHERE_TOOLS_DIR", file.path(Sys.getenv("HOME"), "gwasphere_tools")), "MAGMA/NCBI38/NCBI38.gene.loc")),
              textInput(ns("geneset"),
                "Gene set file (MSigDB .gmt)",
                value = file.path(Sys.getenv("GWASPHERE_TOOLS_DIR", file.path(Sys.getenv("HOME"), "gwasphere_tools")), "MAGMA/MSigDB/c5.all.v2025.1.Hs.entrez.gmt")),
              textInput(ns("gtex"),
                "GTEx expression file",
                value = file.path(Sys.getenv("GWASPHERE_TOOLS_DIR", file.path(Sys.getenv("HOME"), "gwasphere_tools")), "MAGMA/GTEx/gtex_v8_ts_general_avg_log2TPM_entrez.txt")),
              numericInput(ns("magma_parallel"), "MAGMA parallel array jobs", value = 8, min = 1, max = 32)
            )
          )
        ),

        # --- Optional: Novelty + PubMed --------------------------------
        wellPanel(
          style = "border-left: 4px solid #e67e22;",
          checkboxInput(ns("run_novelty"),
            span(icon("magnifying-glass"), strong(" Run Novelty / GWAS Catalog + PubMed"), " — literature evidence"),
            value = TRUE),

          conditionalPanel(
            condition = paste0("input['", ns("run_novelty"), "'] == true"),
            div(style = "margin-top:10px;",
              textInput(ns("catalog_tsv"),
                "GWAS Catalog TSV file",
                value = file.path(Sys.getenv("GWASPHERE_TOOLS_DIR", file.path(Sys.getenv("HOME"), "gwasphere_tools")), "GWAS_catalog/gwas-catalog-download-associations-alt-full.tsv")),
          textInput(ns("default_trait_string"),
            "Trait string (single trait only)",
            placeholder = "body height"),
          textInput(ns("trait_map_xlsx"),
            HTML('Trait to EFO ID mapping file <span title="Columns: trait, EFO_ID — accepts .xlsx, .csv, .tsv">ⓘ</span>'),
            value = "", placeholder = "/path/to/trait_map.xlsx / .csv / .tsv"),
          helpText(style = "font-size:12px; color:#888;",
            icon("circle-info"),
            " For a single trait, enter the trait name in the field above.",
            " For multiple traits, provide a mapping file with three columns:",
            " trait (trait code as in .regenie.gz files),",
            " EFO_ID (trait identifier in the GWAS Catalog),",
            " GWAS_catalog_trait (trait name as written in the GWAS Catalog, used for PubMed searches).")
            )
          )
        ),

      ),

      # --- Right column: launch & live log ----------------------------
      column(6,

        wellPanel(
          h4(icon("play"), " Launch pipeline"),

          uiOutput(ns("config_summary")),

          br(),
          actionButton(ns("launch_pipeline"), "Launch downstream",
                       icon = icon("rocket"),
                       class = "btn-success btn-lg",
                       style = "width:100%;"),

          br(),
          uiOutput(ns("launch_status"))
        ),

        wellPanel(
          h4(icon("terminal"), " Live output"),
          p(style = "font-size:12px; color:#666;",
            "Refreshes every 5 seconds. Shows the last 50 lines of the log."),
          verbatimTextOutput(ns("live_log"), placeholder = TRUE),
          uiOutput(ns("step_progress"))
        )
      )
    )
  )
}


downstreamServer <- function(id, pipeline_script, global) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_job <- reactiveVal(NULL)

    # --- Configuration summary ------------------------------------------
    output$config_summary <- renderUI({
      missing <- c()
      if (!nzchar(input$input_dir))  missing <- c(missing, "REGENIE input directory")
      if (!nzchar(input$out_dir))    missing <- c(missing, "Output directory")
      if (!nzchar(input$ref_prefix)) missing <- c(missing, "PLINK reference prefix")
      if (isTRUE(input$run_magma)) {
        if (!nzchar(input$gene_loc)) missing <- c(missing, "MAGMA gene location file")
        if (!nzchar(input$geneset))  missing <- c(missing, "MAGMA gene set file")
        if (!nzchar(input$gtex))     missing <- c(missing, "MAGMA GTEx file")
      }
      if (isTRUE(input$run_novelty)) {
        if (!nzchar(input$catalog_tsv)) missing <- c(missing, "GWAS Catalog TSV")
      }

      if (length(missing) > 0) {
        div(
          style = "background:#fef9e7; border:1px solid #f39c12;
                   border-radius:4px; padding:8px; font-size:12px;",
          icon("triangle-exclamation", style = "color:#e67e22;"),
          strong(" Required fields missing:"),
          tags$ul(lapply(missing, tags$li))
        )
      } else {
        div(
          style = "background:#eafaf1; border:1px solid #27ae60;
                   border-radius:4px; padding:8px; font-size:12px;",
          icon("check", style = "color:#27ae60;"),
          strong(" Configuration complete."), br(),
          span(style = "color:#555;", paste("Input:", input$input_dir)), br(),
          span(style = "color:#555;", paste("Output:", input$out_dir)), br(),
          span(style = "color:#555;",
               paste("Steps:",
                     paste(c(
                       "Core",
                       if (isTRUE(input$run_magma))   "MAGMA",
                       if (isTRUE(input$run_novelty)) "Novelty+PubMed"
                     ), collapse = " · ")))
        )
      }
    })

    # --- Launch pipeline ------------------------------------------------
    # Auto-update gene_loc when genome build changes
    observeEvent(input$genome_build, {
      if (input$genome_build == "hg19") {
        updateTextInput(session, "gene_loc",
          value = file.path(Sys.getenv("GWASPHERE_TOOLS_DIR", file.path(Sys.getenv("HOME"), "gwasphere_tools")), "MAGMA/NCBI37.3/NCBI37.3.gene.loc"))
      } else {
        updateTextInput(session, "gene_loc",
          value = file.path(Sys.getenv("GWASPHERE_TOOLS_DIR", file.path(Sys.getenv("HOME"), "gwasphere_tools")), "MAGMA/NCBI38/NCBI38.gene.loc"))
      }
    })

    observeEvent(input$launch_pipeline, {
      req(nzchar(input$input_dir), nzchar(input$out_dir), nzchar(input$ref_prefix))

      if (!file.exists(pipeline_script)) {
        showNotification(paste("Script not found:", pipeline_script),
                         type = "error", duration = 8)
        return()
      }

      log_path <- file.path(input$out_dir, "logs")
      dir.create(log_path, recursive = TRUE, showWarnings = FALSE)

      sbatch_args <- c(
        "--parsable",
        "--job-name", "post_gwas_pip",
        "--partition", input$partition,
        "--mem",       paste0(input$mem_gb, "G"),
        "--cpus-per-task", as.character(input$cpus),
        "--nodes",     "1",
        "--output",    file.path(log_path, "pipeline.%j.out"),
        "--error",     file.path(log_path, "pipeline.%j.err"),
        "--export", paste0(
          "ALL",
          ",OUT_DIR=",              input$out_dir,
          ",INPUT_DIR=",            input$input_dir,
          ",PIPELINE_DIR=",         dirname(pipeline_script),
          ",REF_PREFIX=",           input$ref_prefix,
          ",P1=",                   input$p1,
          ",PTAG=",                 input$ptag,
          ",STEP1_R2=",             input$step1_r2,
          ",STEP2_R2=",             input$step2_r2,
          ",MERGE_KB=",             input$merge_kb,
          ",MAF=",                  input$maf,
          ",RUN_MAGMA=",            tolower(as.character(isTRUE(input$run_magma))),
          ",RUN_NOVELTY=",          tolower(as.character(isTRUE(input$run_novelty))),

          ",GENE_LOC=",             if (isTRUE(input$run_magma)) input$gene_loc else "",
          ",SNPEFF_GENOME=",       if (input$genome_build == "hg19") "GRCh37.75" else "GRCh38.99",
          ",GENESET=",              if (isTRUE(input$run_magma)) input$geneset  else "",
          ",GTEX=",                 if (isTRUE(input$run_magma)) input$gtex     else "",
          ",MAGMA_PARALLEL=",       if (isTRUE(input$run_magma)) input$magma_parallel else "1",
          ",TRAIT_MAP_XLSX=",       if (isTRUE(input$run_novelty)) input$trait_map_xlsx else "",
          ",CATALOG_TSV=",          if (isTRUE(input$run_novelty)) input$catalog_tsv   else "",
          ",DEFAULT_TRAIT_STRING=", if (isTRUE(input$run_novelty)) gsub(" ", "_", input$default_trait_string) else "",
          ",SLURM_PARTITION=",       input$partition
        ),
        pipeline_script
      )

      result <- tryCatch(
        system2("sbatch", args = sbatch_args, stdout = TRUE, stderr = TRUE),
        error = function(e) paste("Error:", e$message)
      )

      job_id <- gsub("[^0-9]", "", paste(result, collapse = ""))
      current_job(list(id = job_id, log_path = log_path, started = Sys.time()))

      global$last_downstream_outdir <- input$out_dir
      new_jobs <- global$running_jobs
      new_jobs[[job_id]] <- list(id = job_id, name = "downstream", status = "running")
      global$running_jobs <- new_jobs

      showNotification(paste("Job submitted:", paste(result, collapse = " ")),
                       type = "message", duration = 6)
    })

    # --- Job status panel -----------------------------------------------
    output$launch_status <- renderUI({
      job <- current_job()
      if (is.null(job)) return(NULL)
      div(
        style = "background:#eaf4fb; border-radius:4px; padding:8px; font-size:12px;",
        icon("circle-check", style = "color:#3498db;"),
        strong(paste(" Job ID:", job$id)), br(),
        span(style = "color:#555;", paste("Started:", format(job$started, "%H:%M:%S"))), br(),
        span(style = "color:#555;", paste("Log dir:", job$log_path))
      )
    })

    # --- Live log -------------------------------------------------------
    log_timer <- reactiveTimer(5000)

    output$live_log <- renderText({
      log_timer()
      job <- current_job()
      if (is.null(job)) return("No job running yet.")

      log_files <- list.files(job$log_path, pattern = "pipeline\\..*\\.out$",
                               full.names = TRUE)
      if (length(log_files) == 0) return("Waiting for log file...")

      latest <- log_files[which.max(file.mtime(log_files))]
      lines  <- tryCatch(readLines(latest, warn = FALSE), error = function(e) character(0))
      if (length(lines) == 0) return("Log is empty or still being written...")
      paste(tail(lines, 50), collapse = "\n")
    })

    # --- Step progress --------------------------------------------------
    output$step_progress <- renderUI({
      log_timer()
      job <- current_job()
      if (is.null(job)) return(NULL)

      log_files <- list.files(job$log_path, pattern = "pipeline\\..*\\.out$",
                               full.names = TRUE)
      if (length(log_files) == 0) return(NULL)

      latest <- log_files[which.max(file.mtime(log_files))]
      lines  <- tryCatch(readLines(latest, warn = FALSE), error = function(e) character(0))

      step_lines <- grep("^STEP [0-9]", lines, value = TRUE)
      if (length(step_lines) == 0) return(NULL)

      all_steps  <- c("STEP 1", "STEP 1b", "STEP 1c", "STEP 2", "STEP 3", "STEP 4", "STEP 5", "STEP 6", "STEP 6b", "STEP 6c", "STEP 7", "STEP 8", "STEP 8b")
      done_steps <- unique(gsub(":.*", "", step_lines))

      div(
        style = "margin-top:8px;",
        h6("Step progress:"),
        div(
          style = "display:flex; flex-wrap:wrap; gap:4px;",
          lapply(all_steps, function(s) {
            is_done <- any(startsWith(done_steps, s))
            span(
              style = paste0(
                "padding:2px 8px; border-radius:12px; font-size:11px; ",
                if (is_done) "background:#27ae60; color:#fff;"
                else         "background:#ecf0f1; color:#999;"
              ),
              s
            )
          })
        )
      )
    })

  })
}
