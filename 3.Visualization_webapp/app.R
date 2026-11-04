# GWAS pipeline viewer
# Shiny app skeleton for loci, candidate SNPs, novelty, PubMed, SnpEff, MAGMA and LocusZoom outputs
#
# MULTI-DATASET VERSION:
# Each added directory opens as a closable
# tab in a horizontal bar at the top of the page; clicking a tab makes that
# dataset "active" and everything below (sidebar menu items, trait selector,
# all boxes/plots/tables) reflects that dataset -- exactly as the original
# single-dataset app did, just switchable.

DEFAULT_BASE_DIR <- "path/to/downstream_output"

library(shiny)
library(shinydashboard)
library(DT)
library(data.table)
library(dplyr)
library(stringr)
library(plotly)
library(ggplot2)
library(forcats)
library(jsonlite)
library(digest)

source(file.path("R", "ai_cache.R"), local = TRUE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x


# Helper theme for all plots
plot_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(),
      panel.grid.minor = element_blank()
    )
}

# =========================================================
# HELPERS
# =========================================================
safe_fread <- function(path) {

  if (is.null(path) || is.na(path) || !file.exists(path)) {
    return(data.table())
  }

  tryCatch(

    {
      fread(path, skip = "VARIABLE")
    },

    error = function(e1) {

      tryCatch(

        {
          fread(path)
        },

        error = function(e2) {
          data.table()
        }

      )
    }

  )
}


list_traits_from_mh_plots <- function(mh_dir) {
  if (!dir.exists(mh_dir)) return(character())

  files <- list.files(
    mh_dir,
    pattern = "_mh(_maffilt)?\\.html$",
    full.names = FALSE
  )

  files <- files[!grepl("^combined_mh", files)]

  traits <- files |>
    sub("_mh_maffilt\\.html$", "", x = _) |>
    sub("_mh\\.html$", "", x = _) |>
    unique() |>
    sort()

  traits
}

# Fallback trait discovery: tries mh_plots first, then the loci table,
# then the qqplot folder, then MAGMA gene analysis outputs.
# Never errors even if none of the directories/files exist yet.
derive_traits <- function(p, loci_dt) {

  tr <- list_traits_from_mh_plots(p$mh_dir)

  if (length(tr) == 0 && nrow(loci_dt) > 0 && "trait" %in% names(loci_dt)) {
    tr <- sort(unique(as.character(loci_dt$trait)))
  }

  if (length(tr) == 0 && dir.exists(p$qq_dir)) {
    f <- list.files(p$qq_dir, pattern = "\\.png$", full.names = FALSE)
    if (length(f) > 0) {
      tr <- sort(unique(sub("_qqplot\\.png$", "", f, ignore.case = TRUE)))
    }
  }

  if (length(tr) == 0) {
    gdir <- file.path(p$magma_dir, "gene_analysis")
    if (dir.exists(gdir)) {
      f <- list.files(gdir, pattern = "\\.genes\\.symbol\\.tsv$", full.names = FALSE)
      if (length(f) > 0) {
        tr <- sort(unique(sub("\\.genes\\.symbol\\.tsv$", "", f)))
      }
    }
  }

  tr
}

# Generic fuzzy file finder used for QQ plots / LocusZoom images:
# looks in `dir` for a file whose name contains the (sanitized) trait name,
# optionally further filtered by `must_contain`. Returns NA if dir or file
# is missing -- never errors.
find_trait_asset <- function(dir, trait, must_contain = NULL, ext_pattern = "\\.png$") {

  if (is.null(dir) || !dir.exists(dir) || is.null(trait) || !nzchar(trait)) {
    return(NA_character_)
  }

  safe_trait <- gsub("[^A-Za-z0-9_.-]", "_", trait)
  files <- list.files(dir, pattern = ext_pattern, full.names = FALSE)
  if (length(files) == 0) return(NA_character_)

  hit <- files[grepl(safe_trait, files, fixed = TRUE)]
  if (length(hit) == 0) return(NA_character_)

  if (!is.null(must_contain)) {
    hit2 <- hit[grepl(must_contain, hit, ignore.case = TRUE)]
    if (length(hit2) > 0) hit <- hit2
  }

  hit[1]
}

file_if_exists <- function(...) {
  f <- file.path(...)
  if (file.exists(f)) f else NA_character_
}

read_trait_file <- function(dir, trait, suffix) {
  f <- file.path(dir, paste0(trait, suffix))
  safe_fread(f)
}

format_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA, formatC(x, format = "e", digits = 3))
}

shorten_text <- function(x, n = 280) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  ifelse(nchar(x) > n, paste0(substr(x, 1, n), "..."), x)
}

find_p_col <- function(dt) {
  if (nrow(dt) == 0) return(NA_character_)
  nms <- names(dt)
  clean <- toupper(gsub("[^A-Za-z0-9]+", "", nms))

  priority <- c("P", "PVAL", "PVALUE", "CANDIDATEP", "PJOINT", "PSELF", "PCOND", "PFULL")
  hit <- match(priority, clean, nomatch = 0)
  if (any(hit > 0)) return(nms[hit[hit > 0][1]])

  p_like <- grep("^P", clean)
  if (length(p_like) > 0) return(nms[p_like[1]])

  NA_character_
}

find_label_col <- function(dt, candidates) {
  nms <- names(dt)
  clean <- toupper(gsub("[^A-Za-z0-9]+", "", nms))
  cand_clean <- toupper(gsub("[^A-Za-z0-9]+", "", candidates))
  hit <- match(cand_clean, clean, nomatch = 0)
  if (any(hit > 0)) return(nms[hit[hit > 0][1]])
  nms[1]
}

# Convert data.table columns that htmltools/DT cannot safely serialize.
# In particular, fread() and joins can leave integer64, list, IDate or other
# classed columns that trigger: "is.character(txt) is not TRUE".
normalise_table_for_display <- function(dt) {
  x <- as.data.table(copy(dt))

  clean_names <- as.character(names(x))
  clean_names[is.na(clean_names) | !nzchar(clean_names)] <- "column"
  setnames(x, make.unique(clean_names, sep = "_"))

  collapse_cell <- function(value) {
    if (is.null(value) || length(value) == 0L) return("")
    value <- unlist(value, recursive = TRUE, use.names = FALSE)
    value <- as.character(value)
    value <- value[!is.na(value) & nzchar(value)]
    paste(value, collapse = "; ")
  }

  for (nm in names(x)) {
    value <- x[[nm]]

    if (is.list(value) || is.data.frame(value) || is.matrix(value)) {
      value <- vapply(seq_len(nrow(x)), function(i) {
        if (is.matrix(value)) collapse_cell(value[i, , drop = TRUE])
        else collapse_cell(value[[i]])
      }, character(1))
    } else if (inherits(value, "integer64")) {
      value <- as.character(value)
    } else if (inherits(value, c("POSIXct", "POSIXlt"))) {
      value <- format(value, "%Y-%m-%d %H:%M:%S", tz = "UTC")
    } else if (inherits(value, c("Date", "IDate", "ITime", "difftime"))) {
      value <- as.character(value)
    } else if (is.factor(value)) {
      value <- as.character(value)
    } else if (!is.atomic(value)) {
      value <- as.character(value)
    }

    # Some DT/htmltools version combinations call htmlEscape() directly on
    # numeric/logical vectors when column filters are enabled. htmlEscape()
    # requires character input, so make the display copy unconditionally text.
    # The source reactive data remains unchanged for calculations and plots.
    if (!is.character(value)) value <- as.character(value)
    value[is.na(value)] <- ""

    set(x, j = nm, value = value)
  }

  as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
}

# Wraps datatable() so empty datasets render a friendly message and values are
# normalised before DT/htmltools serialisation.
safe_datatable <- function(dt, ..., empty_msg = "No data available yet.") {
  shiny::validate(shiny::need(!is.null(dt) && nrow(dt) > 0, as.character(empty_msg)))
  DT::datatable(normalise_table_for_display(dt), ...)
}

  build_paths <- function(BASE_DIR) {
  # Read REGENIE input directory:
  # 1. Use manual override from sidebar if provided
  # 2. Fall back to path saved in regenie_input_dir.txt by the pipeline
  regenie_override <- if (exists("input") && !is.null(input$regenie_override_dir) &&
                          nzchar(trimws(input$regenie_override_dir))) {
    trimws(input$regenie_override_dir)
  } else {
    NULL
  }

  regenie_input_file <- file.path(BASE_DIR, "regenie_input_dir.txt")
  regenie_input_dir <- if (!is.null(regenie_override) && dir.exists(regenie_override)) {
    regenie_override
  } else if (file.exists(regenie_input_file)) {
    trimws(readLines(regenie_input_file, warn = FALSE)[1])
  } else {
    NULL
  }

  list(
    loci = file.path(BASE_DIR, "loci", "ALL_TRAITS.loci_table.tsv"),
    loci_rsid = file.path(BASE_DIR, "loci", "ALL_TRAITS.loci_table_with_rsid.tsv"),
    candidate_dir = file.path(BASE_DIR, "ld", "candidate_snps"),
    ld_dir = file.path(BASE_DIR, "ld"),
    gwas_catalog_dir = file.path(BASE_DIR, "novelty_check"),
    pubmed_dir = file.path(BASE_DIR, "novelty_check", "PubMed"),
    snpeff_dir = file.path(BASE_DIR, "SnpEff"),
    snpeff_prop = file.path(BASE_DIR, "SnpEff", "all_trait_prop.tsv"),
    magma_dir = file.path(BASE_DIR, "magma"),
    locuszoom_dir = file.path(BASE_DIR, "LZplot"),
    mh_dir = file.path(BASE_DIR, "mh_plots"),
    gene_mh_dir = file.path(BASE_DIR, "mh_plots_gene"),
    qq_dir = file.path(BASE_DIR, "qqplot"),
    tophits_dir = if (!is.null(regenie_input_dir)) file.path(regenie_input_dir, "tophits") else NULL,
    regenie_input_dir = regenie_input_dir,
    run_params = file.path(BASE_DIR, "run_params.json")
  )
}

find_snp_col <- function(dt) {
  if (nrow(dt) == 0) return(NA_character_)
  nms   <- names(dt)
  clean <- toupper(gsub("[^A-Za-z0-9]+", "", nms))
  priority <- toupper(gsub("[^A-Za-z0-9]+", "", c("candidateSNP", "rsID", "SNP", "ID")))
  hit <- match(priority, clean, nomatch = 0)
  if (any(hit > 0)) return(nms[hit[hit > 0][1]])
  NA_character_
}

# =========================================================
# UI
# =========================================================

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(disable = TRUE),
dashboardSidebar(
  menuItem(HTML("&nbsp; Documentation"), tabName = "documentation", icon = icon("circle-info")),
    tags$div(
      style = "padding: 10px 15px;",
      tags$h5(tags$strong("Add dataset")),
      textInput(
        "new_base_dir",
        "Results directory",
        placeholder = DEFAULT_BASE_DIR,
        width = "100%"
      ),
        helpText(style = "font-size:11px; color:#aaa; margin-top:-8px;",
        icon("circle-info"),
        "Same directory as 'Post-GWAS output directory' in the Run Post-GWAS analysis tab."),
      textInput(
            "new_label",
            "Tab label (optional)",
            value = "",
            width = "100%"
      ),
      actionButton("add_dataset", "Add dataset", icon = icon("plus"), width = "80%"),
      uiOutput("add_dataset_feedback"),
      br(),
      tags$details(
        tags$summary(
          style = "font-size:12px; color:#aaa; cursor:pointer; padding: 4px 0;",
          icon("folder-open"), " Override REGENIE results dir"
        ),
        div(
          style = "padding: 6px 0;",
          textInput(
            "regenie_override_dir",
            NULL,
            value = "",
            placeholder = "/path/to/regenie/results",
            width = "100%"
          ),
          p(style = "font-size:11px; color:#aaa; margin:0;",
            "Optional. Use if the top hits path from the pipeline is broken.
             Leave empty to use the path saved during the pipeline run.")
        )
      )
    ),
    hr(),

    selectInput("trait", "Trait", choices = NULL),
    hr(),

    sidebarMenu(
      id = "tabs",
      menuItem("Global view",            tabName = "global",        icon = icon("chart-area")),
      menuItem("Overview",               tabName = "overview",      icon = icon("chart-line")),
      menuItem("Loci",                   tabName = "loci",          icon = icon("map-marker-alt")),
      menuItem("Candidate SNPs",         tabName = "candidate",     icon = icon("dna")),
      menuItem("Novelty / GWAS Catalog", tabName = "novelty",       icon = icon("book")),
      menuItem("PubMed evidence",        tabName = "pubmed",        icon = icon("file-medical")),
      menuItem("SnpEff annotation",      tabName = "snpeff",        icon = icon("microscope")),
      menuItem("MAGMA",                  tabName = "magma",         icon = icon("project-diagram")),
      menuItem("AI Reports",             tabName = "ai_reports",    icon = icon("robot")),
      menuItem("Cross-trait overlap",    tabName = "compare_trait", icon = icon("code-branch")),
      menuItem("Cross-dataset overlap",  tabName = "compare",       icon = icon("layer-group")),
      menuItem("HTML report",       tabName = "share_report",  icon = icon("file-export"))
    )
  ),
  dashboardBody(

    tags$head(tags$style(HTML("
      th.dt-nowrap { white-space: nowrap !important; }
      /* MAGMA flex alignment */
      .magma-row { display:flex; flex-wrap:wrap; margin:0 -7.5px 15px; }
      .magma-row > .col-sm-6 { padding:0 7.5px; display:flex; flex-direction:column; }
      .magma-row > .col-sm-6 > .box { flex:1; }
      /* Fix sidebar colors to match master navbarPage (#2c3e50) */
      .main-sidebar, .left-side { background-color: #2c3e50 !important; }
      .main-sidebar .sidebar { background-color: #2c3e50 !important; }
      .sidebar-menu > li > a { color: #ecf0f1 !important; }
      .sidebar-menu > li.active > a,
      .sidebar-menu > li:hover > a { background-color: #1a252f !important; color: #fff !important; }
      .sidebar-menu > li > a > .fa { color: #ecf0f1 !important; }
      .treeview-menu > li > a { color: #bdc3c7 !important; background: #243342 !important; }
      .content-wrapper, .right-side { background-color: #ecf0f5 !important; }
      .skin-black .sidebar a { color: #ecf0f1 !important; }
      /* Remove top margin since header is disabled */
      .content-wrapper, .main-sidebar { margin-top: 0 !important; }
      .sidebar { padding-top: 0 !important; }
      /* Inputs inside sidebar */
      .main-sidebar label, .main-sidebar .control-label { color: #ecf0f1 !important; }
      .main-sidebar input, .main-sidebar select { color: #333 !important; background: #fff !important; }
      /* Report button at bottom */
      .sidebar-report-btn { 
        position: absolute; bottom: 0; width: 100%; 
        padding: 10px 15px; background: #1a252f;
        border-top: 1px solid #3d5166;
      }

      .dataset-tab-bar {
        display: flex;
        flex-wrap: wrap;
        gap: 4px;
        padding: 8px 15px;
        background: #ecf0f5;
        border-bottom: 1px solid #d2d6de;
        margin-bottom: 15px;
      }
      .dataset-tab {
        display: inline-flex;
        align-items: center;
        padding: 6px 12px;
        border-radius: 4px 4px 0 0;
        background: #fff;
        border: 1px solid #d2d6de;
        border-bottom: none;
        cursor: pointer;
        font-size: 13px;
        color: #444;
        user-select: none;
      }
      .dataset-tab:hover { background: #f4f4f4; }
      .dataset-tab.active {
        background: #3c8dbc;
        color: #fff;
        font-weight: bold;
      }
      .dataset-tab .close-btn {
        margin-left: 10px;
        font-weight: bold;
        opacity: 0.7;
      }
      .dataset-tab .close-btn:hover { opacity: 1; }
      .dataset-tab-empty { padding: 8px 15px; color: #888; font-style: italic; }

      .sidebar-menu > li:first-child > a > .fa,
      .sidebar-menu > li:first-child > a > svg {
        margin-left: 8px !important;
      }

    "))),

    # Horizontal, closable dataset tab bar -- lets you switch between all
    # results directories loaded via the sidebar "Add dataset" form.
    uiOutput("dataset_tabs_bar"),

      tabItems(
        tabItem(
          tabName = "documentation",

          fluidRow(
            box(
              width = 12,
              title = "Pipeline documentation",
              status = "primary",
              solidHeader = TRUE,

              h3("Overview"),
              p("This Shiny app summarizes the main outputs of the GWAS downstream analysis pipeline. Add one or more results directories from the sidebar; each becomes a tab in the bar at the top so you can switch between and compare different runs."),

              h3("1. GWAS summary statistics"),
              p("Input GWAS summary statistics from REGENIE are used as the starting point. Manhattan and QQ plots are generated globally and per trait."),

              h3("2. LD clumping and locus definition"),
              p("Genomic risk loci identified using a FUMA-inspired pipeline implemented in PLINK (Watanabe et al., Nat Commun 2017). Independent significant SNPs: P < 5x10-8, MAF >= 0.001, r2 < 0.6. Lead SNPs identified with secondary clumping (r2 < 0.1). Regions within 250 kb merged into genomic risk loci."),

              h3("3. Candidate SNPs"),
              p("Variants in LD with an independent significant SNP (r2 >= 0.6) and nominal association (P < 0.05). Independent significant SNPs themselves are included as candidates with r2 = 1."),

              h3("4. Novelty / GWAS Catalog"),
              p("Lead SNPs assessed for novelty via GWASLab querying the GWAS Catalog (p <= 5x10-8). Novelty defined as no known association within 500 kb in the same-trait subset (matched by EFO ID). Cross-trait associations are summarised for each lead SNP."),

              h3("5. PubMed evidence"),
              p("PubMed searches retrieve literature evidence for the nearest gene of each lead SNP. For known SNPs the search includes the trait name; for novel SNPs it focuses on functional terms (expression, pathway, human phenotype)."),

              h3("6. SnpEff annotation"),
              p("Candidate SNPs annotated with SnpEff to summarise predicted functional consequences and annotation proportions across variant categories."),

              h3("7. LocusZoom regional plots"),
              p("Regional association plots generated for each independent significant SNP, showing LD coloring and gene track within a 500 kb window."),

              h3("8. MAGMA gene, gene-set and tissue analysis"),
              p("Gene-based association analysis aggregating SNP-level evidence at the gene level. Gene-set analysis (MSigDB) and tissue expression analysis (GTEx v8) identify enriched biological pathways and relevant tissues. MAGMA runs in parallel with the main pipeline."),

              h3("9. AI reports (on-demand)"),
              p("AI-powered interpretation available on-demand from the Explore Results tab. After reviewing results, select the traits of interest and generate targeted AI reports — using a local LLM (Ollama). No data is sent externally."),

              h3("How to use the app"),
              tags$ul(
                tags$li("Enter a results directory in the sidebar and click 'Add dataset'; repeat for every dataset you want to compare."),
                tags$li("Each dataset appears as a tab at the top of the page -- click a tab to make it the active dataset, click the x to remove it."),
                tags$li("Choose a trait from the sidebar."),
                tags$li("Use Global view for combined Manhattan plots."),
                tags$li("Use Overview for trait-specific Manhattan and QQ plots."),
                tags$li("Click rows in Loci to open corresponding LocusZoom plots."),
                tags$li("Any section for which output files are not yet available will simply appear empty -- the app will not crash.")
              )
            )
          )
        ),
      tabItem(
        tabName = "global",
        fluidRow(
          box(width = 12, background = "light-blue",
            p(style = "margin:0; font-size:13px;",
              icon("circle-info"),
              " Shows a combined Manhattan plot across all loaded traits.
                When multiple traits are present, variants are coloured by trait.
                Use the ", strong("Overview"), " tab for trait-specific Manhattan
                and QQ plots.")
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Combined Manhattan plot",
            status = "primary",
            solidHeader = TRUE,

            uiOutput("global_combined_mh")
          )
        ),

        fluidRow(
          box(
            width = 12,
            title = "Combined Manhattan plot (MAF >= 0.01)",
            status = "info",
            solidHeader = TRUE,

            uiOutput("global_combined_mh_maffilt")
          )
        )
      ),
      tabItem(
        tabName = "overview",
        # 1. Run parameters (collapsible)
        fluidRow(
          box(width = 12, title = "Run parameters", status = "info", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              uiOutput("run_params_ui"))
        ),
        # 2. Value boxes on top
        fluidRow(
          valueBoxOutput("n_loci"),
          valueBoxOutput("n_leads"),
          valueBoxOutput("n_candidates")
        ),
        # 3. QQ plot full width, image centered and small
        fluidRow(
          box(width = 12, title = "QQ plot", status = "warning", solidHeader = TRUE,
              p(style = "font-size:13px; color:#666; margin-bottom:8px;",
                "Observed vs expected p-values. Deviation from diagonal indicates signal or inflation."),
              div(style = "display:flex; justify-content:center;",
                uiOutput("overview_trait_qq")))
        ),
        # 3. Manhattan plots
        fluidRow(
          box(width = 12, title = "Trait-specific Manhattan plot",
              status = "primary", solidHeader = TRUE,
              uiOutput("overview_trait_mh"))
        ),
        fluidRow(
          box(width = 12, title = "Trait-specific Manhattan plot (MAF >= 0.01)",
              status = "info", solidHeader = TRUE,
              uiOutput("overview_trait_mh_maffilt"))
        ),
        fluidRow(
          box(
            width = 12,
            title = "Top hits (REGENIE suggestive associations)",
            status = "success",
            solidHeader = TRUE,
            p(style = "font-size:13px; color:#666;",
              "Suggestive associations from REGENIE top hits (−log10P > 5, ~200kb clumping).
               These are pre-computed by REGENIE and shown here as an exploratory overview
               before genome-wide significant loci are defined by the downstream pipeline."),
            DTOutput("tophits_table")
          )
        )
      ),
      tabItem(
        tabName = "loci",
        fluidRow(
          box(width = 12, background = "light-blue",
            p(style = "margin:0; font-size:13px;",
              icon("circle-info"),
              " Genomic risk loci defined by LD clumping.
                Each locus groups one or more independent significant SNPs
                within a merged LD-based region. Click a row to display the
                corresponding LocusZoom regional plot.")
          )
        ),
        fluidRow(
          box(width = 12, title = "Genomic risk loci", status = "primary", solidHeader = TRUE,
          p("Click one row to display the corresponding LocusZoom plot."),
              DTOutput("loci_table"))
        ),
        fluidRow(
          box(
            width = 12,
            title = "Selected locus - LocusZoom plot",
            status = "info",
            solidHeader = TRUE,
            uiOutput("loci_locuszoom_plot")
          )
        )
      ),
      tabItem(
        tabName = "candidate",
        fluidRow(
          box(width = 12, background = "light-blue",
            p(style = "margin:0; font-size:13px;",
              icon("circle-info"),
              " All variants in LD with an independent significant SNP
                and showing nominal association (P < PTAG). These include the index
                SNPs themselves (r² = 1) and their proxies. Useful for fine-mapping
                and functional annotation of associated regions.")
          )
        ),
        fluidRow(
          box(width = 12, title = "Candidate SNPs", status = "primary", solidHeader = TRUE,
              DTOutput("candidate_table"))
        )
      ),
      tabItem(
        tabName = "novelty",
        fluidRow(
          box(width = 12,
              title = "Trait-specific novelty",
              status = "primary", solidHeader = TRUE,
              uiOutput("novelty_params_ui"),
              DTOutput("novel_table"))
        ),
        fluidRow(
          box(width = 12, title = "Other trait associations", status = "warning", solidHeader = TRUE,
              p(style = "font-size:13px; color:#666; margin-bottom:8px;",
                "GWAS Catalog associations from other traits found within 500 kb of each lead SNP.
                 Useful for identifying pleiotropic effects and cross-trait genetic overlap."),
              DTOutput("other_trait_table"))
        )
      ),
      tabItem(
        tabName = "pubmed",
        fluidRow(
          box(width = 12, background = "light-blue",
            p(style = "margin:0; font-size:13px;",
              icon("circle-info"),
              " PubMed searches retrieve literature evidence for genes near lead SNPs
                to support biological interpretation of associated loci.
                Search strategies differ depending on whether the SNP is novel or known
                (see documentation in the Run Post-GWAS analysis tab for details).")
          )
        ),
        fluidRow(
          box(width = 12, title = "PubMed evidence", status = "primary", solidHeader = TRUE,
              DTOutput("pubmed_table"))
        )
      ),
      tabItem(
        tabName = "snpeff",
        fluidRow(
          box(width = 12, background = "light-blue",
            p(style = "margin:0; font-size:13px;",
              icon("circle-info"),
              " Functional annotation of candidate SNPs using SnpEff.
                The chart shows the proportion of variants falling into each
                functional category (e.g. missense, synonymous, intergenic).
                Helps prioritise potentially functional variants for follow-up.")
          )
        ),
        fluidRow(
          box(width = 12, title = "SnpEff annotation proportions", status = "primary", solidHeader = TRUE,
              plotlyOutput("snpeff_prop_plot", height = "500px"))
        )
      ),
      tabItem(
        tabName = "magma",
        fluidRow(
          box(width = 12, background = "light-blue",
            p(style = "margin:0; font-size:13px;",
              icon("circle-info"),
              " Gene-based, gene-set and tissue enrichment analyses using MAGMA.
                Gene analysis aggregates SNP-level p-values within gene boundaries.
                Gene-set analysis tests enrichment in biological pathways (MSigDB).
                Tissue analysis identifies relevant tissues using GTEx expression data.")
          )
        ),
        fluidRow(
          box(width = 12, title = "Gene-based Manhattan plot", status = "primary", solidHeader = TRUE,
              p(style = "font-size:13px; color:#666; margin-bottom:8px;",
                "Gene-level association p-values from MAGMA for the selected trait."),
              uiOutput("magma_gene_mh"))
        ),
        tags$div(class = "magma-row",
          box(width = 6, title = "Top MAGMA genes", status = "primary", solidHeader = TRUE,
              sliderInput("top_genes_n", "Number of genes", min = 5, max = 50, value = 20, step = 5),
              plotlyOutput("magma_gene_plot", height = "500px")),
          box(width = 6, title = "Gene analysis table", status = "primary", solidHeader = TRUE,
              DTOutput("magma_gene_table"))
        ),
        tags$div(class = "magma-row",
          box(width = 6, title = "MAGMA tissue enrichment", status = "info", solidHeader = TRUE,
              sliderInput("top_tissues_n", "Number of tissues", min = 5, max = 50, value = 20, step = 5),
              plotlyOutput("magma_tissue_plot", height = "500px")),
          box(width = 6, title = "Tissue expression table", status = "info", solidHeader = TRUE,
              DTOutput("magma_tissue_table"))
        ),
        tags$div(class = "magma-row",
          box(width = 6, title = "MAGMA gene-set enrichment", status = "success", solidHeader = TRUE,
              sliderInput("top_genesets_n", "Number of gene sets", min = 5, max = 50, value = 20, step = 5),
              plotlyOutput("magma_geneset_plot", height = "500px")),
          box(width = 6, title = "Gene-set analysis table", status = "success", solidHeader = TRUE,
              DTOutput("magma_geneset_table"))
        )
      ),
      tabItem(
        tabName = "ai_reports",
        fluidRow(
          box(
            width = 12,
            title = "AI-supported interpretation",
            status = "primary",
            solidHeader = TRUE,
            p(
              "Prepare evidence-grounded insights for the active dataset only. ",
              "Generate the selected trait or all traits; SLURM results are cached separately per dataset."
            ),
            uiOutput("ai_request_summary"),
            uiOutput("ai_job_status"),
            fluidRow(
              column(4,
                textInput("ai_partition",
                  HTML('SLURM partition <span title="GPU recommended for speed, but CPU partitions also work (slower)">ⓘ</span>'),
                  value = "gpu",
                  width = "100%")
              )
            ),
            br(),
            actionButton("prepare_ai_insights", "Generate report for the selected trait", icon = icon("brain")),
            actionButton("prepare_all_ai_insights", "Generate report for all traits", icon = icon("layer-group")),
            actionButton("refresh_ai_status", "Refresh status", icon = icon("rotate")),
            actionButton("regenerate_ai_insights", "Regenerate selected trait", icon = icon("arrows-rotate")),
            p(style = "font-size:12px; color:#888; margin-top:8px;",
              icon("circle-info"),
              " AI reports are generated as background SLURM jobs. Once ready, they are automatically included in the HTML report — go to the ",
              strong("HTML report"), " tab to download it."),
            br(),
            checkboxInput(
              "include_ai_insights",
              "Include available AI insights in the HTML report",
              value = TRUE
            ),
            uiOutput("ai_insights_preview")
          )
        )
      ),

      tabItem(
        tabName = "share_report",
        fluidRow(
          box(
            width = 12,
            title = "Download HTML report",
            status = "success",
            solidHeader = TRUE,
            p("Create a standalone HTML report from the datasets currently loaded in the app."),
            textInput(
              "report_title",
              "Report title",
              value = "GWAS Results Report",
              width = "100%"
            ),
            checkboxGroupInput(
              "report_datasets", "Datasets to include",
              choices = NULL, selected = NULL
            ),
            actionLink("report_select_all_datasets", "Select all datasets"),
            tags$hr(),
            radioButtons(
              "report_trait_scope", "Traits to include",
              choices = c("All available traits" = "all", "Currently selected trait" = "selected"),
              selected = "all", inline = TRUE
            ),
            checkboxGroupInput(
              "report_sections", "Sections to include",
              choices = c(
                "Manhattan and QQ plots" = "plots",
                "Loci table and LocusZoom plots" = "loci",
                "Candidate SNP tables" = "candidates",
                "Novelty / GWAS Catalog" = "novelty",
                "PubMed evidence" = "pubmed",
                "SnpEff annotation" = "snpeff",
                "MAGMA gene analysis" = "magma_genes",
                "MAGMA gene-set analysis" = "magma_genesets",
                "MAGMA tissue analysis" = "magma_tissues",
                "Cross-trait candidate SNP overlap" = "trait_overlap",
                "Cross-dataset candidate SNP overlap" = "overlap",
                "AI-supported insights" = "ai"
              ),
              selected = c("plots", "loci", "candidates", "ai")
            ),
            downloadButton("download_quarto_report", "Download report bundle (.zip)")
          )
        )
      ),

      tabItem(
      tabName = "compare_trait",
      fluidRow(
        box(width = 12, background = "light-blue",
          p(style = "margin:0; font-size:13px;",
            icon("circle-info"),
            " Compares candidate SNPs across all traits for the active dataset.",
            " Variants appearing in more than one trait may indicate pleiotropic effects.",
            " Requires at least two traits in the active dataset.")
        )
      ),
      tags$div(class = "magma-row",
        box(width = 6, title = "Candidate SNPs shared between traits",
            status = "primary", solidHeader = TRUE,
            p(style = "font-size:13px; color:#666;",
              "For the active dataset, shows candidate SNPs present in more than one trait."),
            plotlyOutput("trait_overlap_counts", height = "400px")),
        box(width = 6, title = "Detail", status = "info", solidHeader = TRUE,
            DTOutput("trait_overlap_table"))
      )
    ),
      tabItem(
      tabName = "compare",
      fluidRow(
        box(width = 12, background = "light-blue",
          p(style = "margin:0; font-size:13px;",
            icon("circle-info"),
            " Compares candidate SNPs across all loaded datasets for the selected trait.
              Useful when analysing the same trait in different cohorts or with different
              parameters — highlights which variants replicate across runs.
              Load at least two datasets to use this view.")
        )
      ),
      tags$div(class = "magma-row",
        box(width = 6, title = "Candidate SNPs shared between loaded datasets (for the current trait)",
            status = "primary", solidHeader = TRUE,
            p(style = "font-size:13px; color:#666;",
              "For the selected trait, count how many of the loaded datasets contain
               each candidate SNP (matched by SNP ID)."),
            plotlyOutput("candidate_overlap_counts", height = "400px")),
        box(width = 6, title = "Detail for SNP", status = "info", solidHeader = TRUE,
            DTOutput("candidate_overlap_table"))
      )
      )
    )
  )
)

# =========================================================
# SERVER
# =========================================================
server <- function(input, output, session) {

  # ---------------------------------------------------------
  # MULTI-DATASET STATE
  # ---------------------------------------------------------
  # ids_order: character vector, order in which datasets were added (drives
  #            the order of the tab bar).
  # datasets_map: named list, id -> list(label = ..., base_dir = ...).
  # active_id: id of the dataset currently shown in the dashboard.
  ids_order    <- reactiveVal(character())
  datasets_map <- reactiveVal(list())
  active_id    <- reactiveVal(NULL)

  observeEvent(input$add_dataset, {
    req(nzchar(input$new_base_dir))

    bd <- normalizePath(input$new_base_dir, mustWork = FALSE)

    if (!dir.exists(bd)) {
      output$add_dataset_feedback <- renderUI({
        tags$p(style = "color:#dd4b39; font-size:12px;", paste("Directory does not exist:", bd))
      })
      return(invisible(NULL))
    }

    new_id <- paste0("ds_", as.integer(Sys.time()), "_", sample.int(99999, 1))
    label  <- if (nzchar(input$new_label)) input$new_label else basename(bd)

    m <- datasets_map()
    m[[new_id]] <- list(id = new_id, label = label, base_dir = bd)
    datasets_map(m)

    ids_order(c(ids_order(), new_id))
    active_id(new_id)

    output$add_dataset_feedback <- renderUI({
      tags$p(style = "color:#00a65a; font-size:12px;", paste("Added:", label))
    })

    updateTextInput(session, "new_base_dir", value = "")
    updateTextInput(session, "new_label", value = "")
  })

  # Fired by the inline onclick handler on each tab (see dataset_tabs_bar).
  observeEvent(input$tab_click, {
    active_id(input$tab_click)
  })

  # Fired by the inline onclick handler on each tab's "x" button.
  observeEvent(input$remove_click, {
    id_to_remove <- input$remove_click

    m <- datasets_map()
    m[[id_to_remove]] <- NULL
    datasets_map(m)

    remaining <- setdiff(ids_order(), id_to_remove)
    ids_order(remaining)

    if (identical(active_id(), id_to_remove)) {
      active_id(if (length(remaining) > 0) remaining[[length(remaining)]] else NULL)
    }
  })

  output$dataset_tabs_bar <- renderUI({
    ids <- ids_order()
    m <- datasets_map()

    if (length(ids) == 0) {
      return(tags$div(class = "dataset-tab-bar",
        tags$span(class = "dataset-tab-empty", "No dataset loaded yet -- add one from the sidebar.")
      ))
    }

    tabs <- lapply(ids, function(id) {
      is_active <- identical(active_id(), id)
      label <- m[[id]]$label

      tags$div(
        class = paste("dataset-tab", if (is_active) "active" else ""),
        onclick = sprintf("Shiny.setInputValue('tab_click', '%s', {priority: 'event'})", id),
        tags$span(label),
        tags$span(
          class = "close-btn",
          HTML("&times;"),
          onclick = sprintf(
            "event.stopPropagation(); Shiny.setInputValue('remove_click', '%s', {priority: 'event'});",
            id
          )
        )
      )
    })

    tags$div(class = "dataset-tab-bar", tabs)
  })

  # All datasets are selected by default, including newly added datasets.
  # After a manual checkbox change, preserve the user's selection.
  report_dataset_selection <- reactiveVal(NULL)
  report_datasets_initialized <- reactiveVal(FALSE)
  report_datasets_manual <- reactiveVal(FALSE)

  observeEvent(input$report_datasets, {
    if (!isolate(report_datasets_initialized())) return()
    selected <- as.character(input$report_datasets %||% character())
    # Ignore input echoes of our own updateCheckboxGroupInput calls.
    if (setequal(selected, isolate(report_dataset_selection() %||% character()))) return()
    report_dataset_selection(selected)
    report_datasets_manual(TRUE)
  }, ignoreNULL = FALSE)

  observe({
    ids <- ids_order()
    m <- datasets_map()
    if (!length(ids)) return()

    labels <- vapply(ids, function(id) {
      as.character(m[[id]]$label %||% id)
    }, character(1))

    manual <- isolate(report_datasets_manual())
    previous <- isolate(report_dataset_selection())
    selected <- if (!manual) ids else intersect(previous %||% character(), ids)

    report_dataset_selection(selected)
    report_datasets_initialized(TRUE)
    updateCheckboxGroupInput(
      session, "report_datasets",
      choices = stats::setNames(ids, make.unique(labels)),
      selected = selected
    )
  })

  observeEvent(input$report_select_all_datasets, {
    ids <- ids_order()
    report_datasets_manual(FALSE)
    report_dataset_selection(ids)
    updateCheckboxGroupInput(session, "report_datasets", selected = ids)
  })

  # base_dir of the currently active dataset (NULL / req() blocks downstream
  # reactives until at least one dataset is loaded and selected).
    

  active_base_dir <- reactive({
    id <- active_id()
    req(id)
    m <- datasets_map()
    req(m[[id]])
    m[[id]]$base_dir
  })

  paths <- reactive({
    build_paths(active_base_dir())
  })


    # Register resource paths for static assets, but only if the folder
    # actually exists yet -- addResourcePath() would error otherwise.
    # Re-registered every time the active dataset changes, so the static
    # aliases ("mh_plots", "lz_plots", "qq_plots") always point at whichever
    # dataset is currently on screen.
    observeEvent(active_base_dir(), {
      p <- paths()

      if (dir.exists(p$mh_dir))        addResourcePath("mh_plots", p$mh_dir)
      if (dir.exists(p$gene_mh_dir))  addResourcePath("gene_mh_plots", p$gene_mh_dir)
      if (dir.exists(p$locuszoom_dir)) addResourcePath("lz_plots", p$locuszoom_dir)
      if (dir.exists(p$qq_dir))        addResourcePath("qq_plots", p$qq_dir)

      # Set AI cache root to ai_reports/ inside the active downstream folder
      ai_cache <- file.path(active_base_dir(), "ai_reports")
      dir.create(ai_cache, recursive = TRUE, showWarnings = FALSE)
      Sys.setenv(GWAS_AI_CACHE_ROOT = ai_cache)
    })

    #load loci
    loci <- reactive({
    req(active_id())
    p <- paths()
    dt <- safe_fread(p$loci_rsid)
    if (nrow(dt) == 0) dt <- safe_fread(p$loci)
    dt
  })

  observeEvent(active_base_dir(), {
    tr <- derive_traits(paths(), loci())

    current <- input$trait
    keep_current <- !is.null(current) && nzchar(current) && current %in% tr

    updateSelectInput(
      session,
      "trait",
      choices = tr,
      selected = if (keep_current) current
                else if (length(tr) > 0) tr[1]
                else NULL
    )
  })

  trait_loci <- reactive({
    req(active_id())
    req(input$trait)

    dt_all <- loci()

    if (nrow(dt_all) == 0) return(data.table())

    dt <- dt_all[trait == input$trait]
    dt
  })


    output$global_combined_mh <- renderUI({
      req(active_id())

      html_file <- "combined_mh.html"
      f <- file.path(paths()$mh_dir, html_file)

      if (!file.exists(f)) {
        return(tags$p("Combined Manhattan plot not found."))
      }

      tags$iframe(
        # cache-busting query string: without it the browser reuses the
        # cached copy from a previous dataset, even after switching tabs
        src = paste0("mh_plots/", html_file, "?v=", as.integer(Sys.time())),
        width = "100%",
        height = "500px",
        style = "border:none;"
      )
    })

    output$global_combined_mh_maffilt <- renderUI({
      req(active_id())

      html_file <- "combined_mh_maffilt.html"
      f <- file.path(paths()$mh_dir, html_file)

      if (!file.exists(f)) {
        return(tags$p("Combined MAF-filtered Manhattan plot not found."))
      }

      tags$iframe(
        src = paste0("mh_plots/", html_file, "?v=", as.integer(Sys.time())),
        width = "100%",
        height = "500px",
        style = "border:none;"
      )
    })

  candidate_dt <- reactive({
    req(active_id())
    req(input$trait)
    dt <- read_trait_file(paths()$candidate_dir, input$trait, ".candidate_snps.tsv")
    if (nrow(dt) == 0) return(dt)

    nms   <- names(dt)
    clean <- toupper(gsub("[^A-Za-z0-9]+", "", nms))
    hit   <- which(clean == "CANDIDATEP")

    if (length(hit) > 0) {
      old_col <- nms[hit[1]]
      new_col <- "-log10P (candidate SNP)"

      dt[[new_col]] <- round(-log10(suppressWarnings(as.numeric(dt[[old_col]]))), digits=3)
      dt[[old_col]] <- NULL

      # rimette la nuova colonna circa dove stava la vecchia
      new_order <- append(setdiff(names(dt), new_col), new_col, after = hit[1] - 1)
      setcolorder(dt, new_order)
    }

    dt
  })

  # ── Cross-trait overlap ──────────────────────────────────────────────────────
  trait_overlap_raw <- reactive({
    req(active_id())
    p <- paths()
    all_tr <- derive_traits(p, loci())
    shiny::validate(shiny::need(length(all_tr) >= 2,
      "Load a dataset with at least 2 traits to compare."))

    pieces <- lapply(all_tr, function(tr) {
      dt <- read_trait_file(p$candidate_dir, tr, ".candidate_snps.tsv")
      if (!nrow(dt)) return(NULL)
      snp_col <- find_snp_col(dt)
      if (is.na(snp_col)) return(NULL)
      p_col <- find_p_col(dt)
      p_vals <- if (!is.na(p_col)) suppressWarnings(as.numeric(dt[[p_col]])) else rep(NA_real_, nrow(dt))
      data.table(
        snp     = as.character(dt[[snp_col]]),
        trait   = tr,
        log10P  = round(-log10(pmax(p_vals, .Machine$double.xmin)), 3)
      )
    })
    pieces <- Filter(Negate(is.null), pieces)
    shiny::validate(shiny::need(length(pieces) >= 2, "Not enough traits with candidate SNPs."))
    unique(rbindlist(pieces, use.names = TRUE))
  })

  trait_overlap_summary <- reactive({
    raw <- copy(trait_overlap_raw())
    raw[, n_traits := uniqueN(trait), by = snp]
    raw[, traits_list := paste(sort(unique(trait)), collapse = ", "), by = snp]
    raw[, log10P_list := paste(sort(unique(log10P)), collapse = ", "), by = snp]
    out <- unique(raw[, .(snp, n_traits, traits_list, log10P_list)])
    out <- out[n_traits >= 2]
    setorder(out, -n_traits, snp)
    out
  })

  output$trait_overlap_counts <- renderPlotly({
    dt <- trait_overlap_summary()
    cnt <- dt[, .N, by = n_traits]
    setorder(cnt, n_traits)
    plot_ly(data = cnt, x = ~as.character(n_traits), y = ~N,
      type = "bar", marker = list(color = "#8e44ad"),
      text = ~as.character(N), textposition = "outside",
      hovertemplate = "Traits: %{x}<br>SNPs: %{y}<extra></extra>") |>
      layout(
        title = list(text = "Cross-trait candidate SNP overlap"),
        xaxis = list(title = "Number of traits sharing the SNP", type = "category"),
        yaxis = list(title = "Number of candidate SNPs", rangemode = "tozero",
          range = list(0, max(cnt$N) * 1.15)),
        margin = list(t = 60, b = 60)
      )
  })

  output$trait_overlap_table <- renderDT({
    safe_datatable(
      trait_overlap_summary(),
      filter = "top",
      options = list(scrollX = TRUE, pageLength = 10),
      rownames = FALSE,
      empty_msg = "No shared candidate SNPs found between traits."
    )
  })

  # ── Cross-dataset overlap ─────────────────────────────────────────────────
  candidate_overlap_raw <- reactive({
    req(input$trait)
    ids <- ids_order()
    shiny::validate(shiny::need(length(ids) >= 2, "Load at least 2 datasets to compare candidate SNPs."))
    m <- datasets_map()

    pieces <- lapply(ids, function(id) {
      info <- m[[id]]
      p    <- build_paths(info$base_dir)
      dt   <- read_trait_file(p$candidate_dir, input$trait, ".candidate_snps.tsv")
      if (nrow(dt) == 0) return(NULL)

      snp_col <- find_snp_col(dt)
      if (is.na(snp_col)) return(NULL)

      p_col <- find_p_col(dt)
      p_values <- if (!is.na(p_col) && p_col %in% names(dt)) {
        suppressWarnings(as.numeric(dt[[p_col]]))
      } else {
        rep(NA_real_, nrow(dt))
      }

      data.table(
        snp = as.character(dt[[snp_col]]),
        dataset = as.character(info$label),
        candidateP = round(-log10(pmax(p_values, .Machine$double.xmin)), digits = 3)
      )
    })


    pieces <- pieces[!vapply(pieces, is.null, logical(1))]

    shiny::validate(shiny::need(length(pieces) > 0,
                  "No overlapping candidate SNPs."))

    unique(rbindlist(pieces, use.names = TRUE))
  })

  candidate_overlap_summary <- reactive({
    raw <- copy(candidate_overlap_raw())
    raw[, n_datasets := uniqueN(dataset), by = snp]
    raw[, datasets_list := paste(sort(unique(dataset)), collapse = ", "), by = snp]
    raw[, `-log10P_list` := paste(sort(unique(candidateP)), collapse = ", "), by = snp]

    out <- unique(raw[, .(snp, n_datasets, datasets_list,`-log10P_list`)])
    setorder(out, -n_datasets, `-log10P_list`, snp)
    out
  })
  novelty_dt <- reactive({
    req(active_id())
    req(input$trait)
    read_trait_file(paths()$gwas_catalog_dir, input$trait, "_novel.tsv")
  })

  other_traits_dt <- reactive({
    req(active_id())
    req(input$trait)
    read_trait_file(paths()$gwas_catalog_dir, input$trait, "_othertraits_in_window.tsv")
  })

  pubmed_dt <- reactive({
    req(active_id())
    req(input$trait)
    dt <- read_trait_file(paths()$pubmed_dir, input$trait, "_pubmed.tsv")
    if (nrow(dt) == 0) return(dt)

    if ("pmid" %in% names(dt)) {
      dt[, pmid_link := paste0(
        '<a href="https://pubmed.ncbi.nlm.nih.gov/', pmid,
        '/" target="_blank">', pmid, '</a>'
      )]
    }

    if ("abstract" %in% names(dt)) {
      dt[, abstract_full := as.character(abstract)]
      dt[, abstract := paste0(
        '<span title="', htmltools::htmlEscape(as.character(abstract_full)), '">',
        htmltools::htmlEscape(shorten_text(abstract_full, 280)),
        '</span>'
      )]
    }

    # Put clickable PMID in place of raw PMID
    if ("pmid_link" %in% names(dt)) {
      dt[, pmid := pmid_link]
      dt[, pmid_link := NULL]
    }

    if ("abstract_full" %in% names(dt)) {
      dt[, abstract_full := NULL]
    }

    dt
  })

  snpeff_dt <- reactive({
    req(active_id())
    req(input$trait)
    # Adjust suffix to your refined output filename if different
    f1 <- file.path(paths()$snpeff_dir, paste0(input$trait, "_snpeff_refined.tsv"))
    f2 <- file.path(paths()$snpeff_dir, "output", paste0(input$trait, "_snpeff_refined.tsv"))
    if (file.exists(f1)) safe_fread(f1) else safe_fread(f2)
  })

  snpeff_prop_dt <- reactive({
    req(active_id())
    dt <- safe_fread(paths()$snpeff_prop)
    if (nrow(dt) == 0) return(dt)
    if ("trait" %in% names(dt)) dt <- dt[trait == input$trait]
    dt
  })

  magma_gene_dt <- reactive({
    req(active_id())
    req(input$trait)
    safe_fread(file.path(paths()$magma_dir, "gene_analysis", paste0(input$trait, ".genes.symbol.tsv")))
  })

  magma_geneset_dt <- reactive({
    req(active_id())
    req(input$trait)
    safe_fread(file.path(paths()$magma_dir, "gene_set_analysis", paste0(input$trait, ".gsa.out")))
  })

  magma_tissue_dt <- reactive({
    req(active_id())
    req(input$trait)
    safe_fread(file.path(paths()$magma_dir, "tissue_analysis", paste0(input$trait, ".gsa.out")))
  })

  output$n_loci <- renderValueBox({
    valueBox(nrow(trait_loci()), "Risk loci", icon = icon("map-marker-alt"), color = "blue")
  })

  output$n_leads <- renderValueBox({
    dt <- trait_loci()
    n <- if (nrow(dt) == 0 || !"leadSNPs" %in% names(dt)) 0 else length(unique(unlist(strsplit(dt$leadSNPs, ","))))
    valueBox(n, "Lead SNPs", icon = icon("dna"), color = "purple")
  })

  output$n_candidates <- renderValueBox({
    valueBox(nrow(candidate_dt()), "Candidate SNPs", icon = icon("microscope"), color = "green")
  })

  # --- Run parameters (from run_params.json) ---------------------------------
  run_params_data <- reactive({
    req(active_id())
    p <- paths()
    if (!file.exists(p$run_params)) return(NULL)
    tryCatch(jsonlite::read_json(p$run_params), error = function(e) NULL)
  })

  output$run_params_ui <- renderUI({
    params <- run_params_data()
    if (is.null(params)) {
      return(p(style = "color:#888; font-size:13px;",
               "Run parameters not found. Re-run the downstream pipeline to generate them."))
    }

    cl <- params$clumping
    nv <- params$novelty
    mg <- params$magma
    sr <- params$steps_run

    div(
      p(style="font-size:12px; color:#888; margin-bottom:12px;",
        icon("clock"), " Run date: ", strong(params$run_date %||% "N/A")),
      fluidRow(
        column(4,
          h5("Clumping & locus definition"),
          tags$table(class="table table-condensed table-bordered",
            style="font-size:12px;",
            tags$tr(tags$th("Parameter"), tags$th("Value")),
            tags$tr(tags$td("P1 (index SNP threshold)"),   tags$td(cl$P1 %||% "5e-8")),
            tags$tr(tags$td("PTAG (candidate threshold)"), tags$td(cl$PTAG %||% "0.05")),
            tags$tr(tags$td("Step 1 r²"),                  tags$td(cl$STEP1_R2 %||% "0.6")),
            tags$tr(tags$td("Step 2 r²"),                  tags$td(cl$STEP2_R2 %||% "0.1")),
            tags$tr(tags$td("Merge distance (kb)"),        tags$td(cl$MERGE_KB %||% "250")),
            tags$tr(tags$td("MAF filter"),                 tags$td(cl$MAF %||% "0.001"))
          ),
          h5("Steps run"),
          tags$ul(style="font-size:12px;",
            tags$li("MAGMA: ", strong(sr$magma %||% "unknown")),
            tags$li("Novelty + PubMed: ", strong(sr$novelty_pubmed %||% "unknown"))
          )
        ),
        column(4,
          h5("Reference files"),
          tags$table(class="table table-condensed table-bordered",
            style="font-size:12px;",
            tags$tr(tags$th("File"), tags$th("Path")),
            tags$tr(tags$td("PLINK prefix"),
                    tags$td(tags$code(style="font-size:10px;", params$ref_prefix %||% "N/A"))),
            tags$tr(tags$td("GWAS Catalog"),
                    tags$td(tags$code(style="font-size:10px;", basename(nv$catalog_tsv %||% "N/A"))))
          ),
          h5("MAGMA files"),
          tags$table(class="table table-condensed table-bordered",
            style="font-size:12px;",
            tags$tr(tags$th("File"), tags$th("Path")),
            tags$tr(tags$td("Gene locations"),
                    tags$td(tags$code(style="font-size:10px;", basename(mg$gene_loc %||% "N/A")))),
            tags$tr(tags$td("Gene sets"),
                    tags$td(tags$code(style="font-size:10px;", basename(mg$geneset %||% "N/A")))),
            tags$tr(tags$td("GTEx"),
                    tags$td(tags$code(style="font-size:10px;", basename(mg$gtex %||% "N/A")))),
            tags$tr(tags$td("Parallel jobs"),
                    tags$td(mg$parallel %||% "8"))
          )
        ),
        column(4,
          h5("Novelty / GWAS Catalog"),
          tags$table(class="table table-condensed table-bordered",
            style="font-size:12px;",
            tags$tr(tags$th("Parameter"), tags$th("Value")),
            tags$tr(tags$td("Trait map file"),
                    tags$td(tags$code(style="font-size:10px;",
                      if (nzchar(nv$trait_map_xlsx %||% "")) basename(nv$trait_map_xlsx) else "not used"))),
            tags$tr(tags$td("Default trait string"),
                    tags$td(nv$default_trait_string %||% "not specified"))
          ),
          h5("Input / output"),
          tags$table(class="table table-condensed table-bordered",
            style="font-size:12px;",
            tags$tr(tags$th(""), tags$th("Path")),
            tags$tr(tags$td("Input dir"),
                    tags$td(tags$code(style="font-size:10px;", params$input_dir %||% "N/A"))),
            tags$tr(tags$td("Output dir"),
                    tags$td(tags$code(style="font-size:10px;", params$out_dir %||% "N/A")))
          )
        )
      )
    )
  })

  output$novelty_params_ui <- renderUI({
    req(input$trait)
    params <- run_params_data()
    if (is.null(params)) {
      return(p(style="color:#888; font-size:13px;",
               "Run parameters not found. Re-run the downstream pipeline to generate them."))
    }
    nv <- params$novelty
    trait_map <- nv$trait_map_xlsx %||% ""
    def_str   <- nv$default_trait_string %||% ""
    current_trait <- input$trait

    # Try to get the matched term for the current trait
    matched_term <- NULL
    method <- "string"

    if (nzchar(trait_map) && trait_map != "path/to/file/" && file.exists(trait_map)) {
      # User provided a trait map table — read EFO ID for current trait
      tryCatch({
        tm <- readxl::read_excel(trait_map)
        colnames(tm) <- tolower(trimws(colnames(tm)))
        # Find trait column and EFO/catalog trait columns (case-insensitive)
        col_trait   <- names(tm)[tolower(names(tm)) == "trait"][1]
        col_efo     <- names(tm)[tolower(names(tm)) == "efo_id"][1]
        col_catalog <- names(tm)[tolower(names(tm)) %in% c("gwas_catalog_trait", "catalog_trait", "trait_mapped")][1]

        if (!is.na(col_trait)) {
          row <- tm[trimws(as.character(tm[[col_trait]])) == current_trait, ]
          if (nrow(row) > 0) {
            # Prefer GWAS catalog trait name for display
            catalog_label <- if (!is.na(col_catalog)) trimws(as.character(row[[col_catalog]][1])) else NA
            efo           <- if (!is.na(col_efo))     trimws(as.character(row[[col_efo]][1]))     else NA

            display_label <- if (!is.na(catalog_label) && nzchar(catalog_label) && catalog_label != "NA") {
              catalog_label
            } else {
              trimws(as.character(row[[col_trait]][1]))
            }

            if (!is.na(efo) && nzchar(efo) && efo != "NA") {
              matched_term <- paste0(display_label, " (EFO: ", efo, ")")
            } else {
              matched_term <- display_label
            }
            method <- "efo"
          } else {
            matched_term <- if (nzchar(def_str)) def_str else "not specified"
            method <- "string"
          }
        }
      }, error = function(e) NULL)
    }

    if (is.null(matched_term)) {
      # No trait map used — user provided a manual string
      matched_term <- if (nzchar(def_str)) def_str else ""
      method <- "string"
    }

    # If no matching was performed, show warning and hide table
    no_matching <- !nzchar(matched_term)

    if (no_matching) {
      return(tagList(
        div(
          style = "background:#fff3cd; border:1px solid #ffc107; border-radius:6px; padding:12px 16px; font-size:13px; color:#856404;",
          icon("triangle-exclamation"), " ",
          strong("Novelty check was not performed."),
          " No trait string or EFO mapping was provided. ",
          "To enable novelty classification, set a trait string or provide a trait-EFO mapping file when launching the downstream pipeline."
        )
      ))
    }

    icon_col <- if (method == "efo") "#3498db" else "#e67e22"
    icon_name <- if (method == "efo") "table" else "font"
    method_label <- if (method == "efo") "EFO-based matching" else "String-based matching"

    tagList(
      p(
        style = "font-size:13px; color:#666; margin-bottom:8px;",
        paste0(
          "Lead SNPs assessed against the GWAS Catalog (same-trait subset). ",
          "The NOVEL column indicates whether the association is novel: ",
          "TRUE = no known association within 500 kb for this trait in the GWAS Catalog; ",
          "FALSE = already reported in the GWAS Catalog for this trait."
        )
      ),

      p(
        style = "font-size:13px; margin:0;",
        "Novelty assessed for trait ",
        strong(current_trait),
        " matched against the GWAS Catalog using ",
        strong(method_label), ": ",
        tags$code(matched_term)
      )
    )
  })

  output$overview_table <- renderDT({
    dt <- data.table(
      section = c("Risk loci", "Candidate SNPs", "Novelty rows", "Other trait rows", "PubMed rows", "SnpEff rows"),
      n = c(
        nrow(trait_loci()),
        nrow(candidate_dt()),
        nrow(novelty_dt()),
        nrow(other_traits_dt()),
        nrow(pubmed_dt()),
        nrow(snpeff_dt())
      )
    )
    safe_datatable(dt, options = list(pageLength = 10), rownames = FALSE)
  })

  output$overview_trait_mh <- renderUI({
  req(active_id())
  req(input$trait)

  html_file <- paste0(input$trait, "_mh.html")
  f <- file.path(paths()$mh_dir, html_file)

  if (!file.exists(f)) {
    return(tags$p(paste("No Manhattan plot found for trait:", input$trait)))
  }

  tags$iframe(
    src = paste0("mh_plots/", html_file, "?v=", as.integer(Sys.time())),
    width = "100%",
    height = "500px",
    style = "border:none;"
  )
})

  output$overview_trait_mh_maffilt <- renderUI({
    req(active_id())
    req(input$trait)

    html_file <- paste0(input$trait, "_mh_maffilt.html")
    f <- file.path(paths()$mh_dir, html_file)

    if (!file.exists(f)) {
      return(tags$p(paste("No MAF-filtered Manhattan plot found for trait:", input$trait)))
    }

    tags$iframe(
      src = paste0("mh_plots/", html_file, "?v=", as.integer(Sys.time())),
      width = "100%",
      height = "500px",
      style = "border:none;"
    )
  })

  # QQ plot -- PNG, not interactive. File name is matched fuzzily
  # (must contain the trait name and the word "qq"), so it works
  # regardless of the exact naming convention used by the pipeline.
  output$overview_trait_qq <- renderUI({
    req(active_id())
    req(input$trait)

    qdir <- paths()$qq_dir
    fname <- find_trait_asset(qdir, input$trait, must_contain = "qq", ext_pattern = "\\.png$")

    if (is.na(fname)) {
      return(tags$p(paste("No QQ plot found for trait:", input$trait)))
    }

    tags$img(
      src = paste0("qq_plots/", fname, "?v=", as.integer(Sys.time())),
      style = "max-width:420px; width:100%; border:1px solid #ddd; padding:4px;"
    )
  })

  # --- Top hits (from REGENIE tophits folder) --------------------------------
  tophits_dt <- reactive({
    req(active_id())
    req(input$trait)
    p <- paths()

    tophits_dir <- p$tophits_dir
    if (is.null(tophits_dir) || !dir.exists(tophits_dir)) return(data.table())

    # Find file matching trait name
    files <- list.files(tophits_dir,
                        pattern = paste0("^", input$trait, "\\.regenie\\.filtered\\.gz$"),
                        full.names = TRUE)
    if (length(files) == 0) {
      # Try partial match
      files <- list.files(tophits_dir, pattern = input$trait, full.names = TRUE)
    }
    if (length(files) == 0) return(data.table())

    dt <- tryCatch(
      data.table::fread(cmd = paste("zcat", shQuote(files[1])), sep = "\t"),
      error = function(e) data.table()
    )
    if (nrow(dt) == 0) return(dt)

    # Keep only relevant columns
    keep <- intersect(
      c("RSID", "ID", "CHROM", "GENPOS", "ALLELE1", "ALLELE0",
        "A1FREQ", "BETA", "SE", "LOG10P", "GENE_NAME", "GENE_DISTANCE"),
      names(dt)
    )
    dt <- dt[, ..keep]

    # Round numeric columns
    for (col in c("A1FREQ", "BETA", "SE", "LOG10P", "GENE_DISTANCE")) {
      if (col %in% names(dt)) {
        dt[[col]] <- round(suppressWarnings(as.numeric(dt[[col]])), 4)
      }
    }

    # Sort by LOG10P descending
    if ("LOG10P" %in% names(dt)) data.table::setorder(dt, -LOG10P)
    dt
  })

  output$tophits_table <- renderDT({
    safe_datatable(
      tophits_dt(),
      filter = "top",
      options = list(scrollX = TRUE, pageLength = 25),
      rownames = FALSE,
      empty_msg = "No top hits file found for this trait. Make sure the pipeline has run and the REGENIE results directory is correct."
    )
  })

  output$loci_table <- renderDT({
    safe_datatable(
      trait_loci(),
      filter = "top",
      selection = "single",
      options = list(scrollX = TRUE, pageLength = 20),
      rownames = FALSE,
      empty_msg = "No loci table found for this trait yet."
    )
  })

  selected_locus <- reactive({
  dt <- trait_loci()
  row <- input$loci_table_rows_selected

  if (length(row) == 0 || nrow(dt) == 0) return(NULL)

  dt[row[1]]
})

  output$loci_locuszoom_plot <- renderUI({
    sel <- selected_locus()

    if (is.null(sel)) {
      return(tags$p("Select a locus from the table above."))
    }

    trait <- input$trait

    if (!"leadRSIDs" %in% names(sel)) {
      return(tags$p("Column leadRSIDs not found in loci table."))
    }

    lead_snps <- strsplit(as.character(sel$leadRSIDs[1]), ",")[[1]]
    lead_snps <- trimws(lead_snps)
    lead_snps <- lead_snps[nzchar(lead_snps)]

    if (length(lead_snps) == 0) {
      return(tags$p("No lead SNPs found for this locus."))
    }

    lz_dir <- paths()$locuszoom_dir

    if (!dir.exists(lz_dir)) {
      return(tags$p("LocusZoom output folder not found yet."))
    }

    plot_blocks <- lapply(lead_snps, function(lead) {

      safe_lead <- gsub("[^A-Za-z0-9_.-]", "_", lead)

      patterns <- c(
        paste0("^", trait, "_", safe_lead, "\\.html$"),
        paste0("^", trait, "_.*", safe_lead, ".*\\.html$")
      )

      files <- unique(unlist(lapply(patterns, function(pat) {
        list.files(lz_dir, pattern = pat, full.names = TRUE)
      })))

      if (length(files) == 0) {
        return(tags$div(
          style = "margin-bottom:25px;",
          tags$p(strong(paste("Lead SNP:", lead))),
          tags$p("No LocusZoom plot found.")
        ))
      }

      html_file <- basename(files[1])

      tags$div(
        style = "margin-bottom:35px;",
        tags$p(strong(paste("Lead SNP:", lead))),
        tags$p(html_file),
        tags$iframe(
          src = paste0("lz_plots/", html_file, "?v=", as.integer(Sys.time())),
          width = "100%",
          height = "500px",
          style = "border:none;"
        )
      )
    })

    tagList(plot_blocks)
  })

  output$candidate_table <- renderDT({
    safe_datatable(
      candidate_dt(),
      filter = "top",
      selection = "single",
      options = list(scrollX = TRUE, pageLength = 25),
      rownames = FALSE,
      empty_msg = "No candidate SNP file found for this trait yet."
    )
  })

  selected_candidate <- reactive({
    dt <- candidate_dt()
    row <- input$candidate_table_rows_selected
    if (length(row) == 0 || nrow(dt) == 0) return(NULL)
    dt[row[1]]
  })


 output$candidate_locuszoom_plot <- renderUI({

  sel <- selected_candidate()

  if (is.null(sel)) {
    return(tags$p("Select a candidate SNP."))
  }

  trait <- input$trait
  candidate <- as.character(sel$candidateSNP[1])
  index <- if ("bestIndSigSNP" %in% names(sel))
    as.character(sel$bestIndSigSNP[1])
  else
    candidate

  lz_dir <- paths()$locuszoom_dir

  if (!dir.exists(lz_dir)) {
    return(tags$p("LocusZoom output folder not found yet."))
  }

  safe_index <- gsub("[^A-Za-z0-9_.-]", "_", index)
  safe_candidate <- gsub("[^A-Za-z0-9_.-]", "_", candidate)

  patterns <- c(
    paste0("^", trait, "_", safe_index, "\\.png$"),
    paste0("^", trait, "_.*", safe_index, ".*\\.png$"),
    paste0("^", trait, "_", safe_candidate, "\\.png$"),
    paste0("^", trait, "_.*", safe_candidate, ".*\\.png$")
  )

  files <- unique(unlist(lapply(patterns, function(pat) {
    list.files(lz_dir, pattern = pat, full.names = TRUE)
  })))

  if (length(files) == 0) {
    return(tags$div(
      tags$p(strong("No LocusZoom plot found.")),
      tags$p(paste("Selected candidate:", candidate)),
      tags$p(paste("Best independent significant SNP:", index)),
      tags$p("Expected filename pattern similar to: <trait>_<bestIndSigSNP>.png")
    ))
  }

  f <- files[1]

  tags$div(
    tags$p(strong(paste("Candidate:", candidate))),
    tags$p(paste("Best independent significant SNP:", index)),
    tags$p(basename(f)),
    tags$img(
      src = f,
      width = "850px",
      style = "max-width:100%; border:1px solid #ddd; padding:4px;"
    )
  )
})

  output$novel_table <- renderDT({
    dt <- novelty_dt()
    # Hide table if novelty check was not performed (all NOVEL=TRUE with no trait matching)
    params <- run_params_data()
    if (!is.null(params)) {
      nv <- params$novelty
      def_str <- nv$default_trait_string %||% ""
      trait_map <- nv$trait_map_xlsx %||% ""
      no_matching <- !nzchar(def_str) &&
        (!nzchar(trait_map) || trait_map == "path/to/file/" || !file.exists(trait_map))
      if (no_matching) return(safe_datatable(data.table(), rownames = FALSE,
        empty_msg = "Novelty check was not performed — no trait string or EFO mapping was provided."))
    }
    if (nrow(dt) > 0) {
      # Keep only essential columns for readability
      keep <- intersect(
        c("SNPID", "rsID", "CHR", "POS", "GENE_NAME", "GENE_DISTANCE",
          "LOG10P", "BETA", "SE", "EAF", "NOVEL",
          "KNOWN_ID", "KNOWN_AUTHOR", "KNOWN_PUBMED_ID", "DISTANCE_TO_KNOWN"),
        names(dt)
      )
      dt <- dt[, ..keep, with = FALSE]

      # Round numeric columns
      for (col in c("LOG10P", "BETA", "SE", "EAF", "GENE_DISTANCE", "DISTANCE_TO_KNOWN")) {
        if (col %in% names(dt)) dt[[col]] <- round(suppressWarnings(as.numeric(dt[[col]])), 4)
      }
    }
    safe_datatable(
      dt, filter = "top",
      options = list(scrollX = TRUE, pageLength = 20),
      rownames = FALSE,
      empty_msg = "No novelty file found for this trait yet."
    )
  })

  output$other_trait_table <- renderDT({
    safe_datatable(
      other_traits_dt(), filter = "top",
      options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE,
      empty_msg = "No 'other traits in window' file found for this trait yet."
    )
  })

  output$pubmed_table <- renderDT({
    dt <- pubmed_dt()
    shiny::validate(shiny::need(nrow(dt) > 0, "No PubMed file found for this trait yet."))
    dt <- normalise_table_for_display(dt)

    DT::datatable(
      dt,
      filter = "top",
      escape = FALSE,
      options = list(
        scrollX = TRUE,
        pageLength = 20,
        columnDefs = list(
          list(width = '420px', targets = which(names(dt) == "abstract") - 1)
        )
      ),
      rownames = FALSE
    )
  })

  output$snpeff_prop_plot <- renderPlotly({
    dt <- snpeff_prop_dt()
    shiny::validate(shiny::need(nrow(dt) > 0, "No SnpEff proportion file found or no rows for this trait."))
    shiny::validate(shiny::need(all(c("annotation", "proportion") %in% names(dt)), "Expected columns: trait, annotation, N, proportion"))

    dt[, proportion := suppressWarnings(as.numeric(proportion))]
    dt <- dt[is.finite(proportion)]
    dt[, annotation_label := as.character(annotation)]
    setorder(dt, proportion)
    n_display <- if ("N" %in% names(dt)) as.character(dt[["N"]]) else rep("not available", nrow(dt))

    dt[, hover_text := paste0(
      "Annotation: ", annotation_label,
      "<br>N: ", n_display,
      "<br>Proportion: ", round(proportion, 4)
    )]
    dt[, hover_text := as.character(hover_text)]

    plot_ly(
      data = dt,
      x = ~proportion,
      y = ~annotation_label,
      type = "bar",
      orientation = "h",
      marker = list(color = "#966B9D"),
      text = ~hover_text,
      hoverinfo = "text"
    ) |>
      layout(
        title = list(text = paste("SnpEff annotation proportions -", as.character(input$trait))),
        xaxis = list(title = "Proportion"),
        yaxis = list(title = "")
      )
  })

  output$magma_gene_plot <- renderPlotly({
    dt <- magma_gene_dt()
    shiny::validate(shiny::need(nrow(dt) > 0, "No MAGMA gene results found for this trait."))
    p_col <- find_p_col(dt)
    shiny::validate(shiny::need(!is.na(p_col), paste("Could not find a P-value column in MAGMA gene results. Columns found:", paste(names(dt), collapse = ", "))))
    if ("SYMBOL" %in% names(dt)) {
      entrez_col <- if ("ENTREZID" %in% names(dt)) "ENTREZID" else names(dt)[1]
      dt[, gene_display := fifelse(is.na(SYMBOL) | SYMBOL == "", as.character(get(entrez_col)), as.character(SYMBOL))]
    } else {
      dt[, gene_display := as.character(get(names(dt)[1]))]
    }
    dt[, p_numeric := as.numeric(get(p_col))]
    dt <- dt[!is.na(p_numeric) & p_numeric > 0]
    dt[, score := -log10(p_numeric)]
    dt[, gene_display := as.character(gene_display)]
    dt[is.na(gene_display) | gene_display == "NA" | !nzchar(gene_display),
       gene_display := as.character(get(if ("ENTREZID" %in% names(dt)) "ENTREZID" else names(dt)[1]))]
    dt <- dt[!is.na(gene_display) & gene_display != "NA" & nzchar(gene_display)]
    dt <- dt[order(-score)][1:min(.N, input$top_genes_n)]
    setorder(dt, score)
    zstat_display <- if ("ZSTAT" %in% names(dt)) as.character(dt[["ZSTAT"]]) else rep("not available", nrow(dt))
    dt[, hover_text := paste0("Gene: ", gene_display, "<br>ZSTAT: ", zstat_display,
      "<br>P: ", signif(p_numeric, 3), "<br>-log10(P): ", round(score, 3))]
    dt[, hover_text := as.character(hover_text)]
    cat_order_gene <- dt$gene_display
    plot_ly(data = dt, x = ~score, y = ~gene_display, type = "bar", orientation = "h",
      marker = list(color = "#D64550"), text = ~hover_text, hoverinfo = "text", textposition = "none") |>
      layout(title = list(text = as.character(input$trait)),
        xaxis = list(title = "-log10(P)"),
        yaxis = list(title = "", categoryorder = "array", categoryarray = cat_order_gene,
          tickmode = "array", tickvals = cat_order_gene, ticktext = cat_order_gene))
  })

  output$magma_gene_mh <- renderUI({
    req(active_id())
    req(input$trait)
    p <- paths()
    html_file <- paste0(input$trait, "_gene_mh.html")
    f <- file.path(p$gene_mh_dir, html_file)
    if (!file.exists(f)) {
      return(tags$p(style = "color:#888; font-size:13px;",
        "No gene-based Manhattan plot found for this trait. Run the downstream pipeline with MAGMA enabled."))
    }
    tags$iframe(
      src = paste0("gene_mh_plots/", html_file, "?v=", as.integer(Sys.time())),
      width = "100%", height = "500px", style = "border:none;"
    )
  })

  output$magma_gene_table <- renderDT({
    dt <- copy(magma_gene_dt())
    drop_extra <- intersect(c("GENE", "gene_display"), names(dt))
    if (length(drop_extra)) dt[, (drop_extra) := NULL]
    p_col_idx <- which(toupper(gsub("[^A-Za-z0-9]+", "", names(dt))) %in% c("P","PVAL","PVALUE"))
    if (length(p_col_idx)) {
      p_col_name <- names(dt)[p_col_idx[1]]
      dt[, log10P_tmp := round(-log10(suppressWarnings(as.numeric(get(p_col_name)))), 4)]
      setnames(dt, "log10P_tmp", "-log10P")
      new_ord <- append(setdiff(names(dt), "-log10P"), "-log10P", after = p_col_idx[1])
      setcolorder(dt, new_ord)
      dt <- dt[order(suppressWarnings(as.numeric(get(p_col_name))))]
    }
    drop_tmp <- intersect(c("p_numeric", "score", "gene_display", "set_display", "tissue_display"), names(dt))
    if (length(drop_tmp)) dt[, (drop_tmp) := NULL]
    shiny::validate(shiny::need(nrow(dt) > 0, "No MAGMA gene analysis file found for this trait yet."))
    # Convert to data.frame keeping numeric types
    df <- as.data.frame(dt, stringsAsFactors = FALSE)
    DT::datatable(df, filter = "top", rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10,
        columnDefs = list(list(width = "120px", className = "dt-nowrap", targets = which(names(df) == "-log10P") - 1L))))
  })
  output$magma_geneset_plot <- renderPlotly({
    dt <- magma_geneset_dt()
    shiny::validate(shiny::need(nrow(dt) > 0, "No MAGMA gene-set results found for this trait."))

    set_col <- find_label_col(dt, c("VARIABLE", "FULL_NAME", "GENESET", "SET", "NAME"))
    p_col <- find_p_col(dt)
    shiny::validate(shiny::need(!is.na(p_col), paste("Could not find a P-value column in MAGMA gene-set results. Columns found:", paste(names(dt), collapse = ", "))))

    dt[, p_numeric := as.numeric(get(p_col))]
    dt <- dt[!is.na(p_numeric) & p_numeric > 0]
    dt[, score := -log10(p_numeric)]
    dt <- dt[order(-score)][1:min(.N, input$top_genesets_n)]
    dt[, set_display := as.character(get(set_col))]
    setorder(dt, score)
    beta_display <- if ("BETA" %in% names(dt)) round(suppressWarnings(as.numeric(dt[["BETA"]])), 3) else rep(NA_real_, nrow(dt))
    ngenes_display <- if ("NGENES" %in% names(dt)) as.character(dt[["NGENES"]]) else rep("not available", nrow(dt))
    dt[, hover_text := paste0(
      "Gene set: ", set_display,
      "<br>BETA: ", beta_display,
      "<br>NGENES: ", ngenes_display,
      "<br>-log10(P): ", round(score, 3)
    )]
    dt[, hover_text := as.character(hover_text)]
    cat_order_set <- dt$set_display

    plot_ly(data = dt, x = ~score, y = ~set_display, type = "bar", orientation = "h",
      marker = list(color = "#D64550"), text = ~hover_text, hoverinfo = "text", textposition = "none") |>
      layout(title = list(text = as.character(input$trait)),
        xaxis = list(title = "-log10(P)"),
        yaxis = list(title = "", categoryorder = "array", categoryarray = cat_order_set,
          tickmode = "array", tickvals = cat_order_set, ticktext = cat_order_set))
  })

  output$magma_geneset_table <- renderDT({
    dt <- copy(magma_geneset_dt())
    p_col_idx <- which(toupper(gsub("[^A-Za-z0-9]+", "", names(dt))) %in% c("P","PVAL","PVALUE"))
    if (length(p_col_idx)) {
      p_col_name <- names(dt)[p_col_idx[1]]
      dt[, log10P_tmp := round(-log10(suppressWarnings(as.numeric(get(p_col_name)))), 4)]
      setnames(dt, "log10P_tmp", "-log10P")
      new_ord <- append(setdiff(names(dt), "-log10P"), "-log10P", after = p_col_idx[1])
      setcolorder(dt, new_ord)
      dt <- dt[order(suppressWarnings(as.numeric(get(p_col_name))))]
    }
    drop_tmp <- intersect(c("p_numeric", "score", "gene_display", "set_display", "tissue_display"), names(dt))
    if (length(drop_tmp)) dt[, (drop_tmp) := NULL]
    shiny::validate(shiny::need(nrow(dt) > 0, "No MAGMA gene-set analysis file found for this trait yet."))
    # Convert to data.frame keeping numeric types
    df <- as.data.frame(dt, stringsAsFactors = FALSE)
    DT::datatable(df, filter = "top", rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10,
        columnDefs = list(list(width = "120px", className = "dt-nowrap", targets = which(names(df) == "-log10P") - 1L))))
  })
  output$magma_tissue_plot <- renderPlotly({
    dt <- magma_tissue_dt()
    shiny::validate(shiny::need(nrow(dt) > 0, "No MAGMA tissue results found for this trait."))

    tissue_col <- find_label_col(dt, c("VARIABLE", "TISSUE", "TISSUE_NAME", "NAME"))
    p_col <- find_p_col(dt)
    shiny::validate(shiny::need(!is.na(p_col), paste("Could not find a P-value column in MAGMA tissue results. Columns found:", paste(names(dt), collapse = ", "))))

    dt[, p_numeric := as.numeric(get(p_col))]
    dt <- dt[!is.na(p_numeric) & p_numeric > 0]
    dt[, score := -log10(p_numeric)]
    dt <- dt[order(-score)][1:min(.N, input$top_tissues_n)]
    dt[, tissue_display := as.character(get(tissue_col))]
    setorder(dt, score)
    beta_display <- if ("BETA" %in% names(dt)) round(suppressWarnings(as.numeric(dt[["BETA"]])), 3) else rep(NA_real_, nrow(dt))
    dt[, hover_text := paste0(
      "Tissue: ", tissue_display,
      "<br>BETA: ", beta_display,
      "<br>-log10(P): ", round(score, 3)
    )]
    dt[, hover_text := as.character(hover_text)]
    cat_order_tissue <- dt$tissue_display

    plot_ly(data = dt, x = ~score, y = ~tissue_display, type = "bar", orientation = "h",
      marker = list(color = "#D64550"), text = ~hover_text, hoverinfo = "text", textposition = "none") |>
      layout(title = list(text = as.character(input$trait)),
        xaxis = list(title = "-log10(P)"),
        yaxis = list(title = "", categoryorder = "array", categoryarray = cat_order_tissue,
          tickmode = "array", tickvals = cat_order_tissue, ticktext = cat_order_tissue))
  })

  output$magma_tissue_table <- renderDT({
    dt <- copy(magma_tissue_dt())
    p_col_idx <- which(toupper(gsub("[^A-Za-z0-9]+", "", names(dt))) %in% c("P","PVAL","PVALUE"))
    if (length(p_col_idx)) {
      p_col_name <- names(dt)[p_col_idx[1]]
      dt[, log10P_tmp := round(-log10(suppressWarnings(as.numeric(get(p_col_name)))), 4)]
      setnames(dt, "log10P_tmp", "-log10P")
      new_ord <- append(setdiff(names(dt), "-log10P"), "-log10P", after = p_col_idx[1])
      setcolorder(dt, new_ord)
      dt <- dt[order(suppressWarnings(as.numeric(get(p_col_name))))]
    }
    drop_tmp <- intersect(c("p_numeric", "score", "gene_display", "set_display", "tissue_display"), names(dt))
    if (length(drop_tmp)) dt[, (drop_tmp) := NULL]
    shiny::validate(shiny::need(nrow(dt) > 0, "No MAGMA tissue analysis file found for this trait yet."))
    # Convert to data.frame keeping numeric types
    df <- as.data.frame(dt, stringsAsFactors = FALSE)
    DT::datatable(df, filter = "top", rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10,
        columnDefs = list(list(width = "120px", className = "dt-nowrap", targets = which(names(df) == "-log10P") - 1L))))
  })
  output$locuszoom_images <- renderUI({
    req(active_id())
    req(input$trait)

    lz_dir <- paths()$locuszoom_dir
    if (!dir.exists(lz_dir)) return(tags$p("LocusZoom output folder not found yet."))

    files <- list.files(lz_dir, pattern = paste0("^", input$trait, "_.*\\.png$"), full.names = TRUE)
    if (length(files) == 0) return(tags$p("No LocusZoom plots found for this trait."))

    tagList(lapply(files, function(f) {
      tags$div(
        style = "display:inline-block; margin:10px; vertical-align:top;",
        tags$p(basename(f)),
        tags$img(src = f, width = "450px")
      )
    }))
  })

  # ---------------------------------------------------------
  # ASYNCHRONOUS AI INSIGHT CACHE
  # ---------------------------------------------------------
  ai_state <- reactiveValues(last_job_id = NULL)

  selected_dataset_information <- reactive({
    id <- active_id()
    if (is.null(id)) return(list())
    info <- datasets_map()[[id]]
    if (is.null(info)) return(list())
    if (is.null(info$id)) info$id <- id
    list(info)
  })

  ai_target_traits <- reactive({
    info <- selected_dataset_information()
    if (!length(info)) return(character())
    p <- build_paths(info[[1]]$base_dir)
    loci_dt <- safe_fread(if (file.exists(p$loci_rsid)) p$loci_rsid else p$loci)
    traits <- derive_traits(p, loci_dt)
    sort(unique(traits[!is.na(traits) & nzchar(traits)]))
  })

  ai_cache_map <- reactive({
    information <- selected_dataset_information()
    traits <- ai_target_traits()
    if (!length(information) || !length(traits)) return(list())
    tryCatch(
      setNames(lapply(traits, function(trait) ai_cache_dir_for(information, trait)), traits),
      error = function(e) list()
    )
  })

  prepare_ai_traits <- function(traits, force = FALSE) {
    information <- selected_dataset_information()
    req(length(information) > 0)
    traits <- unique(as.character(traits))
    traits <- traits[!is.na(traits) & nzchar(traits)]
    req(length(traits) > 0)

    tryCatch({
      requests <- list()
      for (trait in traits) {
        cache_dir <- ai_cache_dir_for(information, trait)
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        evidence_file <- file.path(cache_dir, "evidence.json")
        assessment_file <- file.path(cache_dir, "assessment.json")
        status_file <- file.path(cache_dir, "status.json")
        status <- read_ai_status(cache_dir)$status %||% "unknown"

        if (!force && (file.exists(assessment_file) || status %in% c("queued", "starting", "running"))) next
        if (force) unlink(c(assessment_file, status_file), force = TRUE)

        write_ai_json(build_ai_evidence(information, trait), evidence_file)
        write_ai_json(list(
          status = "queued", job_id = NULL,
          model = Sys.getenv("OLLAMA_MODEL", "qwen3:8b"),
          trait = trait,
          updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
          error = NULL
        ), status_file)
        requests[[length(requests) + 1L]] <- list(
          evidence_file = evidence_file,
          assessment_file = assessment_file,
          status_file = status_file
        )
      }

      if (!length(requests)) {
        showNotification("All requested traits are already cached or running.", type = "message")
      } else if (length(requests) == 1L) {
        request <- requests[[1]]
        ai_state$last_job_id <- submit_ai_job(
          request$evidence_file, request$assessment_file, request$status_file,
          partition = input$ai_partition
        )
        showNotification(paste("AI job submitted:", ai_state$last_job_id), type = "message", duration = 6)
      } else {
        manifest <- write_ai_array_manifest(requests)
        ai_state$last_job_id <- submit_ai_job_array(manifest, length(requests),
                                                    partition = input$ai_partition)
        showNotification(
          paste("AI job array submitted:", ai_state$last_job_id, "for", length(requests), "traits"),
          type = "message", duration = 8
        )
      }
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  }

  observeEvent(input$prepare_ai_insights, prepare_ai_traits(input$trait, FALSE))
  observeEvent(input$prepare_all_ai_insights, prepare_ai_traits(ai_target_traits(), FALSE))
  observeEvent(input$regenerate_ai_insights, prepare_ai_traits(input$trait, TRUE))

  current_ai_status <- reactive({
    input$refresh_ai_status
    invalidateLater(5000, session)
    trait <- input$trait
    if (is.null(trait) || !nzchar(trait)) return(list(status = "not_generated"))
    cache_dir <- ai_cache_map()[[trait]]
    if (is.null(cache_dir)) return(list(status = "not_generated"))
    read_ai_status(cache_dir)
  })

  output$ai_request_summary <- renderUI({
    req(input$trait)
    information <- selected_dataset_information()
    tags$div(
      tags$p(tags$strong("Active dataset: "), information[[1]]$label),
      tags$p(tags$strong("Selected trait: "), input$trait),
      tags$p(
        tags$strong("All traits: "),
        if (length(ai_target_traits())) paste(ai_target_traits(), collapse = ", ")
        else "No traits available"
      )
    )
  })

  output$ai_job_status <- renderUI({
    status <- current_ai_status()
    state <- status$status %||% "unknown"
    content <- switch(
      state,
      not_generated = list("AI insights have not been prepared.", "#777"),
      queued = list("Waiting in the queue.", "#d58512"),
      starting = list("Ollama is starting.", "#d58512"),
      running = list("The model is generating and validating insights.", "#3c8dbc"),
      completed = list("AI insights are ready and can be included in the report.", "#00a65a"),
      failed = list(paste("AI generation failed:", status$error %||% "unknown error"), "#dd4b39"),
      list("AI status is unknown.", "#777")
    )
    cache_map <- ai_cache_map()
    states <- if (length(cache_map)) vapply(cache_map, function(path) {
      read_ai_status(path)$status %||% "not_generated"
    }, character(1)) else character()
    summary <- if (length(states)) paste(
      sum(states == "completed"), "ready;",
      sum(states %in% c("queued", "starting", "running")), "queued/running;",
      sum(states == "failed"), "failed;",
      sum(states %in% c("not_generated", "unknown")), "not prepared"
    ) else "No report traits found."
    tagList(
      tags$p(style = paste0("font-weight:600;color:", content[[2]], ";"), content[[1]]),
      tags$p(tags$strong("All report traits: "), summary)
    )
  })

  output$ai_insights_preview <- renderUI({
    req(input$trait)
    cache_dir <- ai_cache_map()[[input$trait]]
    if (is.null(cache_dir)) return(NULL)
    assessment <- read_ai_assessment(cache_dir)
    if (is.null(assessment)) return(NULL)

    # Support both old format (insights list) and new format (loci list)
    has_new_format <- !is.null(assessment$loci) && length(assessment$loci) > 0
    has_old_format <- !is.null(assessment$insights) && length(assessment$insights) > 0

    if (!has_new_format && !has_old_format) return(NULL)

    div(
      style = "margin-top:15px;",

      # Summary
      if (!is.null(assessment$summary) && nzchar(assessment$summary)) {
        div(style = "padding:12px 16px; background:#f4f8fb;
                     border-left:4px solid #3c8dbc; border-radius:4px; margin-bottom:16px;",
          tags$b("Summary: "),
          tags$span(assessment$summary %||% "")
        )
      },

      # Key takeaway
      if (!is.null(assessment$key_takeaway) && nzchar(assessment$key_takeaway %||% "")) {
        div(style = "padding:10px 16px; background:#eafaf1;
                     border-left:4px solid #27ae60; border-radius:4px; margin-bottom:16px;",
          tags$b("Key finding: "),
          tags$span(assessment$key_takeaway)
        )
      },

      # Loci (new format)
      if (has_new_format) {
        tagList(
          tags$h5(style="margin-top:16px;", "Locus-by-locus summary"),
          tagList(lapply(assessment$loci, function(locus) {
            novelty_color <- switch(locus$novelty %||% "Unknown",
              "Novel"   = "#27ae60",
              "Known"   = "#e67e22",
              "#888"
            )
            div(style = "background:#fff; border:1px solid #e0e0e0; border-radius:6px;
                         padding:12px 16px; margin-bottom:10px;
                         border-left:4px solid #3498db;",
              fluidRow(
                column(6,
                  tags$b(locus$lead_snp %||% ""),
                  tags$span(style="color:#888; font-size:12px; margin-left:8px;",
                            paste0("-log10P = ", locus$log10p %||% "")), br(),
                  tags$span(style="font-size:13px;",
                    icon("dna"), " Nearest gene: ",
                    tags$b(locus$nearest_gene %||% "N/A"))
                ),
                column(6,
                  tags$span(
                    style = paste0("font-weight:600; color:", novelty_color, ";"),
                    locus$novelty %||% "Unknown"
                  ), br(),
                  tags$span(style="font-size:12px; color:#555;",
                    locus$novelty_detail %||% "")
                )
              ),
              if (!is.null(locus$pubmed_support) && nzchar(locus$pubmed_support %||% "") &&
                  locus$pubmed_support != "No PubMed hits found") {
                p(style="margin:8px 0 0; font-size:13px; color:#555;",
                  icon("book-open"), " ", locus$pubmed_support)
              },
              if (!is.null(locus$interpretation) && nzchar(locus$interpretation %||% "")) {
                p(style="margin:8px 0 0; font-size:13px; font-style:italic; color:#444;",
                  locus$interpretation)
              }
            )
          }))
        )
      },

      # Pathway / tissue enrichment
      if (!is.null(assessment$pathway_enrichment) &&
          nzchar(assessment$pathway_enrichment %||% "")) {
        div(style="margin-top:12px; padding:10px 16px; background:#f5eef8;
                   border-left:4px solid #8e44ad; border-radius:4px;",
          tags$b("Pathway enrichment: "),
          tags$span(assessment$pathway_enrichment)
        )
      },
      if (!is.null(assessment$tissue_enrichment) &&
          nzchar(assessment$tissue_enrichment %||% "")) {
        div(style="margin-top:8px; padding:10px 16px; background:#fef9e7;
                   border-left:4px solid #f39c12; border-radius:4px;",
          tags$b("Tissue enrichment: "),
          tags$span(assessment$tissue_enrichment)
        )
      },

      # Old format fallback
      if (!has_new_format && has_old_format) {
        tagList(
          tags$h5("Insights"),
          tags$ul(lapply(head(assessment$insights, 8), function(item) {
            tags$li(
              tags$strong(paste0(item$section, ": ")),
              htmltools::htmlEscape(as.character(item$text %||% ""))
            )
          }))
        )
      }
    )
  })

  output$download_quarto_report <- downloadHandler(
    filename = function() paste0("GWAS_report_", Sys.Date(), ".zip"),
    content = function(file) {
      withProgress(message = "Creating the HTML report", value = 0, {
      ids <- intersect(ids_order(), input$report_datasets %||% character())
      if (!length(ids)) stop("Select at least one dataset for the HTML report.")
      if (!length(input$report_sections)) stop("Select at least one report section.")
      incProgress(0.10, detail = "Preparing the selected datasets")

      dataset_map <- datasets_map()
      report_title <- input$report_title
      if (is.null(report_title) || !nzchar(trimws(report_title))) {
        report_title <- "GWAS Results Report"
      } else {
        report_title <- trimws(report_title)
      }
      ai_assessments <- list()
      if ("ai" %in% input$report_sections && isTRUE(input$include_ai_insights)) {
        for (id in ids) {
          info <- dataset_map[[id]]
          p <- build_paths(info$base_dir)
          loci_dt <- safe_fread(if (file.exists(p$loci_rsid)) p$loci_rsid else p$loci)
          available <- derive_traits(p, loci_dt)
          wanted <- if (identical(input$report_trait_scope, "selected")) {
            intersect(input$trait %||% character(), available)
          } else available
          for (trait in wanted) {
            cache_dir <- ai_cache_dir_for(list(info), trait)
            if (!identical(read_ai_status(cache_dir)$status, "completed")) next
            assessment <- read_ai_assessment(cache_dir)
            if (is.null(assessment)) next
            assessment$dataset_label <- info$label
            assessment$dataset_id <- id
            assessment$trait <- trait
            ai_assessments[[paste(id, trait, sep = "::")]] <- assessment
          }
        }
      }

      snapshot <- list(
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        report_title = report_title,
        datasets = lapply(ids, function(id) dataset_map[[id]]),
        trait_scope = input$report_trait_scope,
        selected_trait = input$trait,
        sections = input$report_sections,
        ai = list(
          requested = "ai" %in% input$report_sections,
          included = length(ai_assessments) > 0,
          assessments = ai_assessments
        )
      )

      render_dir <- tempfile(pattern = "gwas_report_")
      dir.create(render_dir)
      snapshot_file <- file.path(render_dir, "report_snapshot.rds")
      rendered_file <- file.path(render_dir, "GWAS_report.html")
      template_copy <- file.path(render_dir, "gwas_report.qmd")
      saveRDS(snapshot, snapshot_file)

      if (!file.copy(
        normalizePath(file.path(VIEWER_APP_DIR, "gwas_report.qmd"), mustWork = TRUE),
        template_copy,
        overwrite = TRUE
      )) {
        stop("The report template could not be copied to the temporary directory.")
      }

      incProgress(0.20, detail = "Rendering plots and tables")
      quarto::quarto_render(
        input = template_copy,
        output_file = basename(rendered_file),
        execute_params = list(snapshot_path = snapshot_file),
        quiet = FALSE
      )

      incProgress(0.55, detail = "Creating the ZIP bundle")
      unlink(c(snapshot_file, template_copy))
      zip::zipr(
        zipfile = file,
        files = list.files(render_dir, all.files = FALSE),
        recurse = TRUE,
        compression_level = 9,
        include_directories = TRUE,
        root = render_dir
      )
      incProgress(0.15, detail = "Download ready")
      })
    }
  )


  output$candidate_overlap_counts <- renderPlotly({
    dt <- candidate_overlap_summary()
    cnt <- dt[, .N, by = n_datasets]
    setorder(cnt, n_datasets)
    cnt[, `:=`(
      n_datasets_label = as.character(n_datasets),
      count_label = as.character(N)
    )]

    plot_ly(
      data = cnt,
      x = ~n_datasets_label,
      y = ~N,
      type = "bar",
      marker = list(color = "#3c8dbc"),
      text = ~count_label,
      textposition = "outside",
      hovertemplate = paste0(
        "Datasets: %{x}<br>",
        "Candidate SNPs: %{y}<extra></extra>"
      )
    ) |>
      layout(
        title = list(text = paste("Overlap candidate SNPs -", as.character(input$trait))),
        xaxis = list(title = "Number of datasets in which the SNP appears", type = "category"),
        yaxis = list(title = "Number of candidate SNPs", rangemode = "tozero",
          range = list(0, max(cnt$N) * 1.15)),
        margin = list(t = 70, b = 60)
      )
  })

output$candidate_overlap_table <- renderDT({
  safe_datatable(
    candidate_overlap_summary(),
    filter = "top",
    options = list(scrollX = TRUE, pageLength = 10),
    rownames = FALSE,
    empty_msg = "No overlap data available."
  )
})

}

shinyApp(ui, server)
