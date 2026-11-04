# =============================================================================
# GWASphere — Master App
# =============================================================================

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinyBS)

source("config.R")
source("modules/regenie_module.R")
source("modules/downstream_module.R")
source("modules/viewer_module.R")

.regenie_env <- .load_shiny_app_from_dir(REGENIE_APP_DIR)
.viewer_env  <- .load_shiny_app_from_dir(VIEWER_APP_DIR)

# Export VIEWER_APP_DIR as bash env variable so ai_cache.R
# can find scripts/ with absolute paths
Sys.setenv(VIEWER_APP_DIR = VIEWER_APP_DIR)

# =============================================================================
# UI
# =============================================================================

regenie_doc_panel <- bsCollapse(
  id = "regenie_doc_collapse",
  bsCollapsePanel(
    title = span(icon("book"), " Documentation — nf-gwas & REGENIE"),
    value = "regenie_doc",
    style = "info",
    div(
      p(style = "color:#555; font-size:14px; line-height:1.7;",
        tags$a("nf-gwas", href = "https://github.com/genepi/nf-gwas", target = "_blank"),
        " is a Nextflow pipeline that wraps ",
        tags$a("REGENIE", href = "https://rgcgithub.github.io/regenie/", target = "_blank"),
        ", a whole-genome regression method designed for large-scale GWAS with
         binary and quantitative traits. It efficiently handles population structure
         and relatedness using a whole-genome regression model fitted in step 1,
         followed by per-variant association testing in step 2."
      ),
      tags$ul(style = "color:#555; font-size:14px; line-height:2; padding-left:20px;",
        tags$li("Supports binary and quantitative traits"),
        tags$li("Handles related individuals and population structure"),
        tags$li("Scalable to biobank-sized datasets (100k+ samples)"),
        tags$li("Runs on SLURM, LSF, and other HPC schedulers via Nextflow")
      ),
      div(style = "margin-top:12px; display:flex; gap:12px; flex-wrap:wrap;",
        tags$a(href = "https://github.com/genepi/nf-gwas", target = "_blank",
               class = "btn btn-sm btn-primary", icon("github"), " nf-gwas on GitHub"),
        tags$a(href = "https://genepi.github.io/nf-gwas/", target = "_blank",
               class = "btn btn-sm btn-default", icon("globe"), " nf-gwas documentation"),
        tags$a(href = "https://rgcgithub.github.io/regenie/", target = "_blank",
               class = "btn btn-sm btn-default", icon("globe"), " REGENIE documentation")
      )
    )
  )
)

downstream_doc_panel <- bsCollapse(
  id = "downstream_doc_collapse",
  bsCollapsePanel(
    title = span(icon("book"), " Documentation — Downstream Pipeline"),
    value = "downstream_doc",
    style = "success",
    div(
      style = "padding:4px 0;",
      p(style = "color:#555; font-size:14px; margin-bottom:16px;",
        "The downstream pipeline runs 9 analytical steps on SLURM after REGENIE completes.
         Configure all paths in the form below and click ", strong("Launch downstream"), " to submit."),
      lapply(list(
        list(n="1", color="#3498db", bg="#eaf4fb", title="GWAS Summary Statistics",
             text="Input GWAS summary statistics from REGENIE are used as the starting point. Manhattan and QQ plots are generated globally and per trait."),
        list(n="2", color="#3498db", bg="#eaf4fb", title="Clumping & Locus Definition",
             text=HTML("Genomic risk loci identified using PLINK following a <a href='https://fuma.ctglab.nl' target='_blank'>FUMA</a>-inspired approach (Watanabe et al., Nat Commun 2017). Independent significant SNPs: P &lt; 5&times;10<sup>-8</sup>, MAF &ge; 0.001, r<sup>2</sup> &lt; 0.6. Lead SNPs identified with secondary clumping (r<sup>2</sup> &lt; 0.1). Regions within 250 kb merged into genomic risk loci.")),
        list(n="3", color="#3498db", bg="#eaf4fb", title="Candidate SNPs",
             text=HTML("Variants in LD with an independent significant SNP (r<sup>2</sup> &ge; 0.6) and nominal association (P &lt; 0.05). Independent significant SNPs included as candidates with r<sup>2</sup> = 1.")),
        list(n="4", color="#3498db", bg="#eaf4fb", title="Novelty / GWAS Catalog",
             text=HTML("Lead SNPs assessed for novelty via GWASLab querying the GWAS Catalog (p &le; 5&times;10<sup>-8</sup>). Novelty defined as no known association within 500 kb in same-trait subset (matched by EFO ID).")),
        list(n="5", color="#3498db", bg="#eaf4fb", title="PubMed Evidence",
             text="PubMed searches retrieve literature evidence for the nearest gene of each lead SNP.
               For known SNPs the search includes the trait name; for novel SNPs it focuses on
               functional terms (expression, pathway, human phenotype)."),
        list(n="6", color="#3498db", bg="#eaf4fb", title="SnpEff Annotation",
             text="Candidate SNPs annotated with SnpEff to summarise predicted functional consequences and annotation proportions across variant categories."),
        list(n="7", color="#3498db", bg="#eaf4fb", title="LocusZoom Regional Plots",
             text="Regional association plots generated for each independent significant SNP, showing LD coloring and gene track within a 500 kb window."),
        list(n="8", color="#3498db", bg="#eaf4fb", title="MAGMA Gene, Gene-Set & Tissue Analysis",
             text="Gene-based association analysis aggregating SNP-level evidence at the gene level. Gene-set analysis (MSigDB) and tissue expression analysis (GTEx v8) identify enriched biological pathways and relevant tissues. Runs in parallel with the main pipeline."),
        list(n="9", color="#3498db", bg="#eaf4fb", title="AI Reports (on-demand)",
             text="AI-powered insights is available on-demand from the Explore Results tab.
                   After reviewing results, select the traits of interest and generate
                   targeted AI reports only for those traits — using a local LLM (Ollama).
                   No data is sent externally.")
      ), function(s) {
        div(
          style = paste0("background:#fff; border:1px solid #e8e8e8; border-radius:8px;
                          padding:12px 16px; margin-bottom:10px; border-left:4px solid ", s$color, ";"),
          div(style = "display:flex; align-items:flex-start; gap:12px;",
              div(style = paste0("background:", s$bg, "; color:", s$color, ";
                                  border-radius:50%; width:28px; height:28px; flex-shrink:0;
                                  display:flex; align-items:center; justify-content:center;
                                  font-weight:700; font-size:13px;"), s$n),
              div(
                strong(style = paste0("color:", s$color, "; font-size:14px;"), s$title),
                p(style = "margin:4px 0 0; color:#555; font-size:13px; line-height:1.6;", s$text)
              )
          )
        )
      })
    )
  )
)

results_doc_panel <- bsCollapse(
  id = "results_doc_collapse",
  bsCollapsePanel(
    title = span(icon("book"), " Documentation — How to use the Results Viewer"),
    value = "results_doc",
    style = "warning",
    tags$ul(style = "color:#555; font-size:14px; line-height:2.2; margin:0; padding-left:20px;",
      tags$li("Enter a results directory in the sidebar and click ", strong("Add dataset"), "; repeat for every run you want to compare."),
      tags$li("Each dataset appears as a tab at the top — click to switch, click ", strong("x"), " to remove."),
      tags$li("Choose a trait from the sidebar."),
      tags$li("Use ", strong("Global view"), " for combined Manhattan plots across traits."),
      tags$li("Use ", strong("Overview"), " for trait-specific Manhattan and QQ plots."),
      tags$li("Click rows in ", strong("Loci"), " to open the corresponding LocusZoom plots."),
      tags$li("Sections for which output files are not yet available will appear empty — the app will not crash.")
    )
  )
)

# =============================================================================
# navbarPage
# =============================================================================
ui <- navbarPage(
  title = tags$img(src = "gwasphere_navbar.png",
             height = "34px",
             style = "margin-top:-4px; vertical-align:middle;"),
  id = "main_nav",
  collapsible = TRUE,

  header = tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }
      .navbar { background-color: #2c3e50 !important; border-color: #1a252f; }
      .navbar .navbar-brand, .navbar-nav > li > a { color: #ecf0f1 !important; }
      .navbar-nav > .active > a { background-color: #3498db !important; color: #fff !important; }
      .status-bar {
        position: fixed; bottom: 0; left: 0; right: 0;
        background: #1a252f; color: #ecf0f1;
        padding: 4px 16px; font-size: 12px; z-index: 9999;
        display: flex; align-items: center; gap: 12px;
      }
      .status-bar .job-badge { background: #27ae60; border-radius: 4px; padding: 2px 8px; font-size: 11px; }
      .status-bar .job-badge.running { background: #e67e22; }
      .status-bar .job-badge.failed  { background: #c0392b; }

      /* Remove ALL padding/margin from the Results tab so the sidebar touches the left edge */
      .tab-content { padding: 0 !important; }
      #tab_viewer { padding: 0 !important; margin: 0 !important; }
      #tab_viewer > .container-fluid { padding: 0 !important; margin: 0 !important; }
      .navbar-default .tab-content > .active { padding: 0 !important; }
      /* Override Bootstrap row negative margins inside viewer */
      #tab_viewer .row { margin-left: 0 !important; margin-right: 0 !important; }
    "))
  ),

  footer = uiOutput("global_status_bar"),

# TAB 1 — ABOUT
  tabPanel(
    title = span(icon("circle-info"), " About"),
    value = "tab_about",
    fluidPage(
      style = "max-width:960px; margin:auto; padding:30px 20px;",
      div(
        style = "text-align:center; margin-bottom:40px;",

        tags$img(src = "gwasphere_about.png",
                 style = "height:100px; margin-bottom:16px;"),
        p(style = "font-size:15px; color:#888; font-style:italic; margin-bottom:8px;",
          "From variants to insights"),
        p(style = "font-size:16px; color:#666; max-width:600px; margin:auto;",
          "A unified, fully automated and privacy-first platform for end-to-end GWAS analysis: easy setup, one-click execution and interactive results exploration.")
      ),
      fluidRow(
        style = "display:flex; flex-wrap:wrap;",
        column(6,
          style = "display:flex;",
          div(style = "background:#fff; border:1px solid #e0e0e0; border-radius:10px;
                       padding:24px; margin-bottom:20px; border-top:4px solid #3498db;
                       box-shadow:0 2px 6px rgba(0,0,0,0.06); width:100%;",
            div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
                div(style = "background:#eaf4fb; border-radius:50%; width:42px; height:42px;
                             display:flex; align-items:center; justify-content:center;",
                    icon("dna", style = "color:#3498db; font-size:18px;")),
                h4(style = "margin:0; color:#2c3e50;", "Run GWAS")),
            p(style = "color:#555; font-size:14px; margin-bottom:12px;",
              "Configure and launch genome-wide association analyses using REGENIE,
               orchestrated via Nextflow (nf-gwas) on your SLURM cluster."),
            tags$ul(style = "color:#555; font-size:13px; padding-left:18px;",
              tags$li("Set genotype and phenotype file paths"),
              tags$li("Configure analysis parameters and covariates"),
              tags$li("Submit and monitor SLURM jobs"),
              tags$li("Supports whole-genome regression (step 1 & 2)"))
          )
        ),
        column(6,
          style = "display:flex;",
          div(style = "background:#fff; border:1px solid #e0e0e0; border-radius:10px;
                       padding:24px; margin-bottom:20px; border-top:4px solid #27ae60;
                       box-shadow:0 2px 6px rgba(0,0,0,0.06); width:100%;",
            div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
                div(style = "background:#eafaf1; border-radius:50%; width:42px; height:42px;
                             display:flex; align-items:center; justify-content:center;",
                    icon("gear", style = "color:#27ae60; font-size:18px;")),
                h4(style = "margin:0; color:#2c3e50;", "Run Post-GWAS analysis")),
            p(style = "color:#555; font-size:14px; margin-bottom:12px;",
              "Post-GWAS analysis pipeline running on SLURM.
               Configure all paths via the web app."),
            tags$ul(style = "color:#555; font-size:13px; padding-left:18px;",
              tags$li("Manhattan & QQ plots"),
              tags$li("LD clumping & locus definition"),
              tags$li("GWAS Catalog novelty check"),
              tags$li("Functional variant annotation (SnpEff)"),
              tags$li("Gene-based analysis (MAGMA): genes, gene-sets, tissues"),
              tags$li("PubMed literature search"),
              tags$li("Regional association plots (LocusZoom)"))
          )
        ),
        column(6,
          style = "display:flex;",
          div(style = "background:#fff; border:1px solid #e0e0e0; border-radius:10px;
                       padding:24px; margin-bottom:20px; border-top:4px solid #e67e22;
                       box-shadow:0 2px 6px rgba(0,0,0,0.06); width:100%;",
            div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
                div(style = "background:#fef5ec; border-radius:50%; width:42px; height:42px;
                             display:flex; align-items:center; justify-content:center;",
                    icon("chart-bar", style = "color:#e67e22; font-size:18px;")),
                h4(style = "margin:0; color:#2c3e50;", "Explore results")),
            p(style = "color:#555; font-size:14px; margin-bottom:12px;",
              "Interactive portal to explore GWAS results across multiple datasets
               and traits, with integrated AI-powered insights."),
            tags$ul(style = "color:#555; font-size:13px; padding-left:18px;",
              tags$li("Interactive viewer"),
              tags$li("Compare datasets and traits"),
              tags$li("AI-generated insights (local LLM)"),
              tags$li("HTML interactive report"))
          )
        ),
        column(6,
          style = "display:flex;",
          div(style = "background:#fff; border:1px solid #e0e0e0; border-radius:10px;
                       padding:24px; margin-bottom:20px; border-top:4px solid #8e44ad;
                       box-shadow:0 2px 6px rgba(0,0,0,0.06); width:100%;",
            div(style = "display:flex; align-items:center; gap:12px; margin-bottom:12px;",
                div(style = "background:#f5eef8; border-radius:50%; width:42px; height:42px;
                             display:flex; align-items:center; justify-content:center;",
                    icon("terminal", style = "color:#8e44ad; font-size:18px;")),
                h4(style = "margin:0; color:#2c3e50;", "Log & Job Monitor")),
            p(style = "color:#555; font-size:14px; margin-bottom:12px;",
              "Real-time monitoring of SLURM jobs and pipeline logs,
               without leaving the browser."),
            tags$ul(style = "color:#555; font-size:13px; padding-left:18px;",
              tags$li("Auto-refreshing log file viewer (every 5 sec)"),
              tags$li("Global job status bar on every tab"))
          )
        )
      ),
      div(
        style = "background:#f8f9fa; border:1px solid #dee2e6; border-radius:8px;
                 padding:16px 20px; margin-top:10px;",
        div(style = "display:flex; align-items:center; gap:10px;",
            icon("shield-halved", style = "color:#27ae60; font-size:20px;"),
            div(
              h5(style = "margin:0 0 4px; color:#2c3e50;", "Privacy by design"),
              p(style = "margin:0; font-size:13px; color:#555;",
                "Privacy-first. Built for HPC.")
            )
        )
      )
    )
  ),

  # TAB 2 — REGENIE
  tabPanel(
    title = span(icon("dna"), " Run GWAS"),
    value = "tab_regenie",
    fluidPage(
      style = "padding-top:16px;",
      regenie_doc_panel,
      if (!is.null(.regenie_env)) .regenie_env$ui
      else div(class="alert alert-danger", style="margin:20px;",
               strong("REGENIE app could not be loaded. "),
               "Check REGENIE_APP_DIR in config.R: ", REGENIE_APP_DIR)
    )
  ),

  # TAB 3 — DOWNSTREAM
  tabPanel(
    title = span(icon("gear"), " Run Post-GWAS analysis"),
    value = "tab_downstream",
    fluidPage(
      style = "padding-top:16px;",
      downstream_doc_panel,
      downstreamUI("downstream")
    )
  ),

  # TAB 4 — RESULTS
  # No fluidPage wrapper — direct UI to avoid any Bootstrap padding
  tabPanel(
    title = span(icon("chart-bar"), " Explore results"),
    value = "tab_viewer",
    div(
      style = "margin: 0; padding: 0;",
      div(style = "padding: 8px 15px 0 15px;", results_doc_panel),
      if (!is.null(.viewer_env)) .viewer_env$ui
      else div(class="alert alert-danger", style="margin:20px;",
               strong("Results viewer could not be loaded. "),
               "Check VIEWER_APP_DIR in config.R: ", VIEWER_APP_DIR)
    )
  ),

  # TAB 5 — LOG
  tabPanel(
    title = span(icon("terminal"), " Log"),
    value = "tab_logs",
    fluidPage(
      h3("Job Monitor"),
      fluidRow(
        column(12,
          wellPanel(
            h5("Latest log file"),
            div(style = "display:flex; gap:8px; margin-bottom:8px;",
              textInput("log_dir_manual", NULL,
                placeholder = "/path/to/your/downstream/output/logs",
                width = "100%"),
              actionButton("log_dir_apply", "Load", icon = icon("folder-open"),
                class = "btn-sm btn-default")
            ),
            selectInput("log_file_select", NULL, choices = NULL, width = "100%"),
            verbatimTextOutput("log_content", placeholder = TRUE)
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  global <- reactiveValues(
    last_regenie_outdir    = NULL,
    last_downstream_outdir = NULL,
    running_jobs           = list()
  )

  if (!is.null(.regenie_env) && exists("server", envir = .regenie_env)) {
    tryCatch({
      old_wd <- setwd(REGENIE_APP_DIR)
      on.exit(setwd(old_wd), add = TRUE)
      .regenie_env$server(input, output, session)
    }, error = function(e) message("[regenie] Server error: ", e$message))
  }

  if (!is.null(.viewer_env) && exists("server", envir = .viewer_env)) {
    tryCatch({
      old_wd <- setwd(VIEWER_APP_DIR)
      on.exit(setwd(old_wd), add = TRUE)
      .viewer_env$server(input, output, session)
    }, error = function(e) message("[viewer] Server error: ", e$message))
  }

  downstreamServer("downstream",
                   pipeline_script = PIPELINE_SCRIPT,
                   global          = global)

  status_timer <- reactiveTimer(5000)

  output$global_status_bar <- renderUI({
    status_timer()
    jobs <- global$running_jobs
    if (length(jobs) == 0) {
      div(class = "status-bar", span("No jobs running"))
    } else {
      badges <- lapply(jobs, function(j) {
        span(class = paste("job-badge", j$status), paste0(j$name, " [", j$id, "]"))
      })
      div(class = "status-bar", span("Active jobs: "), badges)
    }
  })

  # Manual log directory input
  observeEvent(input$log_dir_apply, {
    manual_dir <- trimws(input$log_dir_manual)
    if (nzchar(manual_dir) && dir.exists(manual_dir)) {
      log_files <- list.files(manual_dir,
        pattern = "\\.(out|err|log)$",
        full.names = TRUE, recursive = FALSE)
      if (length(log_files) > 0) {
        updateSelectInput(session, "log_file_select",
          choices = setNames(log_files, basename(log_files)))
      } else {
        showNotification("No .out/.err/.log files found in that directory.", type = "warning")
      }
    } else {
      showNotification("Directory not found.", type = "error")
    }
  })

  observe({
    status_timer()
    dirs <- c(global$last_regenie_outdir, global$last_downstream_outdir)
    dirs <- dirs[!vapply(dirs, is.null, logical(1))]
    dirs <- as.character(dirs)
    dirs <- dirs[nzchar(dirs) & dir.exists(dirs)]
    if (length(dirs) == 0) return()
    log_files <- unlist(lapply(dirs, function(d) {
      list.files(d, pattern = "\\.(out|err|log)$", full.names = TRUE, recursive = FALSE)
    }))
    if (length(log_files) > 0) {
      updateSelectInput(session, "log_file_select",
                        choices = setNames(log_files, basename(log_files)))
    }
  })

  output$log_content <- renderText({
    req(input$log_file_select)
    status_timer()
    path <- input$log_file_select
    if (!file.exists(path)) return("File not found.")
    paste(tail(readLines(path, warn = FALSE), 100), collapse = "\n")
  })
}

shinyApp(ui, server)
