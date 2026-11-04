library(rentrez)
library(data.table)
library(readxl)
library(XML)
args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a[1]))) b else a
}

pubmed_fetch <- function(gene,
                         trait = NULL,
                         status = c("known", "novel"),
                         lead_snp = NULL,
                         lead_snp_rsid = NULL,
                         retmax = 20,
                         min_abstract_chars = 400,
                         require_abstract = TRUE,
                         extra_terms_novel = c("function", "expression", "pathway", "human phenotype")) {

  status <- match.arg(status)
  abstract_filter <- if (require_abstract) " AND hasabstract[text]" else ""

  queries <- list()

  if (status == "known") {
    if (is.null(trait) || !nzchar(trait)) {
      stop("For known variants, trait must be provided.")
    }

    queries[["gene_trait"]] <- sprintf(
      "(%s) AND (%s)%s",
      gene,
      trait,
      abstract_filter
    )

    if (!is.null(lead_snp_rsid) && nzchar(lead_snp_rsid) && lead_snp_rsid != "NA") {
      queries[["gene_trait_rsid"]] <- sprintf(
        "(%s) AND (%s) AND (%s)%s",
        gene,
        trait,
        lead_snp_rsid,
        abstract_filter
      )
    }
  }

  if (status == "novel") {
    for (term in extra_terms_novel) {
      label <- gsub("[^A-Za-z0-9]+", "_", term)
      queries[[paste0("novel_", label)]] <- sprintf(
        "(%s) AND (%s)%s",
        gene,
        term,
        abstract_filter
      )
    }
  }

  fetch_one <- function(query, query_type) {
    Sys.sleep(0.5)

    search <- tryCatch(
      rentrez::entrez_search(
        db = "pubmed",
        term = query,
        retmax = retmax,
        sort = "relevance"
      ),
      error = function(e) return(NULL)
    )

    if (is.null(search) || length(search$ids) == 0) return(data.frame())

    xml <- tryCatch(
      rentrez::entrez_fetch(
        db = "pubmed",
        id = search$ids,
        parsed = TRUE,
        rettype = "xml",
        retmode = "xml"
      ),
      error = function(e) return(NULL)
    )

    if (is.null(xml)) return(data.frame())

    get_text <- function(node, path) {
      n <- tryCatch(XML::getNodeSet(node, path), error = function(e) list())
      if (length(n) == 0) return(NA_character_)
      paste(vapply(n, XML::xmlValue, character(1)), collapse = " ")
    }

    articles <- tryCatch(XML::getNodeSet(xml, "//PubmedArticle"), error = function(e) list())
    if (length(articles) == 0) return(data.frame())

    out <- lapply(articles, function(a) {
      data.frame(
        query_type = query_type,
        total_count = as.integer(search$count),
        pmid = get_text(a, ".//PMID"),
        year = get_text(a, ".//ArticleDate/Year"),
        title = get_text(a, ".//ArticleTitle"),
        abstract = get_text(a, ".//Abstract/AbstractText"),
        stringsAsFactors = FALSE
      )
    })

    df <- rbindlist(out, fill = TRUE)

    if (require_abstract) {
      df <- df[!is.na(abstract) & nchar(abstract) >= min_abstract_chars]
    }

    df <- df[!duplicated(pmid)]
    as.data.frame(df)
  }

  res <- rbindlist(
    lapply(names(queries), function(qt) fetch_one(queries[[qt]], qt)),
    fill = TRUE
  )

  if (nrow(res) == 0) return(data.frame())

  res[, lead_snp := lead_snp]
  res[, lead_snp_rsid := lead_snp_rsid]
  res[, mapped_trait := trait %||% ""]
  res[, gene := gene]
  res[, status := status]

  setcolorder(
    res,
    c("lead_snp", "lead_snp_rsid", "mapped_trait", "gene", "status",
      setdiff(names(res), c("lead_snp", "lead_snp_rsid", "mapped_trait", "gene", "status")))
  )

  as.data.frame(res)
}

# =========================================================
# PATHS
# =========================================================
path <- file.path(args[1], "novelty_check")

out_dir <- file.path(args[1], "novelty_check", "PubMed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trait_map_path <- args[2]

use_trait_map <- FALSE
default_trait_string <- args[3]   # fallback string

if (file.exists(trait_map_path)) {

  trait_map <- read_excel(trait_map_path)
  trait_map <- as.data.table(trait_map)

  trait_map <- trait_map[, .(
    TRAIT = trait,
    GWAS_catalog_trait
  )]

  use_trait_map <- TRUE

  message("Using trait mapping file: ", trait_map_path)

} else {

  message("Trait mapping file not found.")
  message("Fallback: using default trait string -> ", default_trait_string)
}
novel_files <- list.files(
  path = path,
  pattern = "_novel.tsv$",
  full.names = TRUE
)

# =========================================================
# MAIN
# =========================================================
for (f in novel_files) {

  trait <- sub("_novel.tsv$", "", basename(f))
  message("Processing trait: ", trait)


  if (use_trait_map) {

  trait_to_search <- trait_map[
    TRAIT == trait,
    GWAS_catalog_trait
  ][1]

  if (is.na(trait_to_search) || !nzchar(trait_to_search)) {
    trait_to_search <- default_trait_string
  }

  } else {

    trait_to_search <- default_trait_string
  }

  if (is.na(trait_to_search) || !nzchar(trait_to_search)) {
    message("  No GWAS_catalog_trait found for ", trait,
            " — will use generic functional search (novel mode) for all loci.")
    trait_to_search <- ""
  }

  loci <- fread(f)

  required_base <- c("SNPID", "rsID", "GENE_NAME")
  missing <- setdiff(required_base, names(loci))
  if (length(missing) > 0) {
    stop("Missing columns in ", f, ": ", paste(missing, collapse = ", "))
  }

  # If NOVEL column is missing or all NA, treat all loci as novel (generic functional search)
  if (!"NOVEL" %in% names(loci) || all(is.na(loci$NOVEL))) {
    message("[INFO] NOVEL column missing or empty in ", basename(f),
            " — no trait string or EFO mapping was provided. Treating all loci as novel (generic literature search).")
    loci[, status := "novel"]
  } else {
    loci[, status := ifelse(as.character(NOVEL) %in% c("TRUE", "True", "true", "1"), "novel", "known")]
  }

  loci <- unique(loci[, .(
    trait = trait,
    lead_snp = SNPID,
    lead_snp_rsid = rsID,
    nearest_gene = GENE_NAME,
    status
  )])

  loci[, nearest_gene := as.character(nearest_gene)]
  loci <- loci[!is.na(nearest_gene) & nzchar(nearest_gene) & nearest_gene != "NA"]

  df_all <- data.frame()

  for (i in seq_len(nrow(loci))) {

    row <- loci[i]
    gene <- as.character(row$nearest_gene)

    docs <- tryCatch({
      if (row$status == "known") {
        pubmed_fetch(
          gene = gene,
          trait = trait_to_search,
          status = "known",
          lead_snp = row$lead_snp,
          lead_snp_rsid = row$lead_snp_rsid
        )
      } else {
        pubmed_fetch(
          gene = gene,
          trait = trait_to_search,
          status = "novel",
          lead_snp = row$lead_snp,
          lead_snp_rsid = row$lead_snp_rsid
        )
      }
    }, error = function(e) {
      message("  PubMed failed: trait=", trait, " gene=", gene, " error=", e$message)
      data.frame()
    })

    if (!is.null(docs) && nrow(docs) > 0) {
      df_all <- rbind(df_all, docs)
    }
  }

  if (nrow(df_all) > 0) {
    df_all$reported_trait <- trait

    df_all <- df_all[, c(
      "reported_trait",
      "mapped_trait",
      "lead_snp",
      "lead_snp_rsid",
      "gene",
      "status",
      "query_type",
      "total_count",
      "pmid",
      "year",
      "title",
      "abstract"
    )]

    write.table(
      df_all,
      file = file.path(out_dir, paste0(trait, "_pubmed.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  } else {
    message("  No PubMed results for ", trait)
  }
}