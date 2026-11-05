AI_PROMPT_VERSION <- "1.0.0"

ai_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  gsub("^_+|_+$", "", x)
}

ai_paths_for <- function(base_dir) {
  list(
    loci = file.path(base_dir, "loci", "ALL_TRAITS.loci_table.tsv"),
    loci_rsid = file.path(base_dir, "loci", "ALL_TRAITS.loci_table_with_rsid.tsv"),
    candidate_dir = file.path(base_dir, "ld", "candidate_snps"),
    gwas_catalog_dir = file.path(base_dir, "novelty_check"),
    pubmed_dir = file.path(base_dir, "novelty_check", "PubMed"),
    snpeff_dir = file.path(base_dir, "SnpEff"),
    magma_dir = file.path(base_dir, "magma")
  )
}

ai_read_trait_loci <- function(base_dir, trait) {
  p <- ai_paths_for(base_dir)
  dt <- safe_fread(p$loci_rsid)
  if (!nrow(dt)) dt <- safe_fread(p$loci)
  if (!nrow(dt) || !"trait" %in% names(dt)) return(data.table::data.table())
  selected_trait <- as.character(trait)
  dt[as.character(get("trait")) == selected_trait]
}

ai_candidate_p_col <- function(dt) {
  if (!nrow(dt)) return(NA_character_)
  clean <- toupper(gsub("[^A-Za-z0-9]+", "", names(dt)))
  hit <- which(clean == "CANDIDATEP")
  if (length(hit)) names(dt)[hit[1]] else NA_character_
}

ai_read_candidates <- function(base_dir, trait) {
  safe_fread(file.path(
    ai_paths_for(base_dir)$candidate_dir,
    paste0(trait, ".candidate_snps.tsv")
  ))
}

ai_cache_root <- function(base_dir = NULL) {
  root <- if (!is.null(base_dir)) file.path(base_dir, "ai_reports") else Sys.getenv("GWAS_AI_CACHE_ROOT", unset = "")
  if (!nzchar(root)) stop("GWAS_AI_CACHE_ROOT is not configured")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  normalizePath(root, mustWork = TRUE)
}

ai_file_signature <- function(path) {
  if (!file.exists(path)) return(list(path = path, exists = FALSE))
  info <- file.info(path)
  list(
    path = normalizePath(path, mustWork = TRUE),
    exists = TRUE,
    size = unname(info$size),
    modified = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z")
  )
}

ai_input_files <- function(base_dir, trait) {
  p <- ai_paths_for(base_dir)
  c(
    if (file.exists(p$loci_rsid)) p$loci_rsid else p$loci,
    file.path(p$candidate_dir, paste0(trait, ".candidate_snps.tsv")),
    file.path(p$gwas_catalog_dir, paste0(trait, "_novel.tsv")),
    file.path(p$gwas_catalog_dir, paste0(trait, "_othertraits_in_window.tsv")),
    file.path(p$pubmed_dir, paste0(trait, "_pubmed.tsv")),
    file.path(p$snpeff_dir, paste0(trait, "_snpeff_refined.tsv")),
    file.path(p$snpeff_dir, "output", paste0(trait, "_snpeff_refined.tsv")),
    file.path(p$magma_dir, "gene_analysis", paste0(trait, ".genes.symbol.tsv")),
    file.path(p$magma_dir, "gene_set_analysis", paste0(trait, ".gsa.out")),
    file.path(p$magma_dir, "tissue_analysis", paste0(trait, ".gsa.out"))
  )
}

make_ai_cache_key <- function(dataset_information, trait) {
  signatures <- lapply(dataset_information, function(info) {
    list(
      dataset_label = info$label,
      base_dir = normalizePath(info$base_dir, mustWork = TRUE),
      files = lapply(ai_input_files(info$base_dir, trait), ai_file_signature)
    )
  })
  digest::digest(
    list(
      trait = trait,
      datasets = signatures,
      model = Sys.getenv("OLLAMA_MODEL", "qwen3:8b"),
      prompt_version = AI_PROMPT_VERSION,
      background_version = Sys.getenv("COHORT_BACKGROUND_VERSION", "none")
    ),
    algo = "sha256"
  )
}

ai_find_p_col <- function(dt) {
  if (!nrow(dt)) return(NA_character_)
  clean <- toupper(gsub("[^A-Za-z0-9]+", "", names(dt)))
  wanted <- c("P", "PVAL", "PVALUE", "PJOINT", "PFULL", "PSELF")
  hit <- match(wanted, clean, nomatch = 0)
  if (any(hit > 0)) names(dt)[hit[hit > 0][1]] else NA_character_
}

ai_top_rows <- function(dt, n = 30L, p_col = NULL) {
  dt <- data.table::as.data.table(dt)
  if (!nrow(dt)) return(dt)
  if (is.null(p_col)) p_col <- ai_find_p_col(dt)
  if (!is.na(p_col) && p_col %in% names(dt)) {
    dt[, ai_p_sort__ := suppressWarnings(as.numeric(get(p_col)))]
    dt[!is.finite(ai_p_sort__), ai_p_sort__ := Inf]
    data.table::setorder(dt, ai_p_sort__)
    dt[, ai_p_sort__ := NULL]
  }
  head(dt, n)
}

ai_records <- function(dt, prefix, n = 30L, p_col = NULL) {
  dt <- ai_top_rows(data.table::copy(dt), n, p_col)
  if (!nrow(dt)) return(list())
  dt[, evidence_id := sprintf("%s_%03d", prefix, seq_len(.N))]
  data.table::setcolorder(dt, c("evidence_id", setdiff(names(dt), "evidence_id")))
  jsonlite::fromJSON(
    jsonlite::toJSON(dt, dataframe = "rows", na = "null", auto_unbox = TRUE, digits = NA),
    simplifyVector = FALSE
  )
}

ai_read_sections <- function(info, trait) {
  p <- ai_paths_for(info$base_dir)
  stable_dataset_id <- paste0(
    "DS",
    toupper(substr(digest::digest(normalizePath(info$base_dir), algo = "xxhash32"), 1L, 8L))
  )
  dataset_prefix <- stable_dataset_id
  snpeff_1 <- file.path(p$snpeff_dir, paste0(trait, "_snpeff_refined.tsv"))
  snpeff_2 <- file.path(p$snpeff_dir, "output", paste0(trait, "_snpeff_refined.tsv"))
  loci <- ai_read_trait_loci(info$base_dir, trait)
  candidates <- ai_read_candidates(info$base_dir, trait)
  pubmed <- safe_fread(file.path(p$pubmed_dir, paste0(trait, "_pubmed.tsv")))
  if (nrow(pubmed) && "abstract" %in% names(pubmed)) {
    pubmed[, abstract := substr(as.character(abstract), 1L, 800L)]
  }
  # Check if novelty classification was performed
  run_params_path <- file.path(info$base_dir, "run_params.json")
  novelty_check_performed <- TRUE
  if (file.exists(run_params_path)) {
    tryCatch({
      rp <- jsonlite::fromJSON(run_params_path)
      def_str   <- if (is.null(rp$novelty$default_trait_string)) "" else rp$novelty$default_trait_string
      trait_map <- if (is.null(rp$novelty$trait_map_xlsx)) "" else rp$novelty$trait_map_xlsx
      if (!nzchar(def_str) && (!nzchar(trait_map) || trait_map == "path/to/file/" || !file.exists(trait_map))) {
        novelty_check_performed <- FALSE
      }
    }, error = function(e) NULL)
  }

  list(
    dataset_id = stable_dataset_id,
    dataset_label = info$label,
    trait = trait,
    novelty_check_performed = novelty_check_performed,
    loci = ai_records(loci, paste0(dataset_prefix, "_LOCUS"), 30L, "minP"),
    candidates = ai_records(candidates, paste0(dataset_prefix, "_CAND"), 60L, ai_candidate_p_col(candidates)),
    novelty = ai_records(
      safe_fread(file.path(p$gwas_catalog_dir, paste0(trait, "_novel.tsv"))),
      paste0(dataset_prefix, "_NOV"), 60L
    ),
    other_traits = ai_records(
      safe_fread(file.path(p$gwas_catalog_dir, paste0(trait, "_othertraits_in_window.tsv"))),
      paste0(dataset_prefix, "_CROSS"), 60L
    ),
    pubmed = ai_records(pubmed, paste0(dataset_prefix, "_PUBMED"), 20L),
    snpeff = ai_records(
      if (file.exists(snpeff_1)) safe_fread(snpeff_1) else safe_fread(snpeff_2),
      paste0(dataset_prefix, "_SNPEFF"), 40L
    ),
    magma_genes = ai_records(
      safe_fread(file.path(p$magma_dir, "gene_analysis", paste0(trait, ".genes.symbol.tsv"))),
      paste0(dataset_prefix, "_MAGMA_GENE"), 20L
    ),
    magma_genesets = ai_records(
      safe_fread(file.path(p$magma_dir, "gene_set_analysis", paste0(trait, ".gsa.out"))),
      paste0(dataset_prefix, "_MAGMA_SET"), 20L
    ),
    magma_tissues = ai_records(
      safe_fread(file.path(p$magma_dir, "tissue_analysis", paste0(trait, ".gsa.out"))),
      paste0(dataset_prefix, "_MAGMA_TISSUE"), 20L
    )
  )
}

ai_record_value <- function(records, candidates) {
  if (!length(records)) return(character())
  values <- unlist(lapply(records, function(record) {
    field <- candidates[candidates %in% names(record)][1]
    if (!length(field) || is.na(field) || is.null(record[[field]])) return(NA_character_)
    as.character(record[[field]])
  }))
  unique(values[!is.na(values) & nzchar(values)])
}

build_ai_comparisons <- function(dataset_evidence, trait) {
  if (length(dataset_evidence) < 2L) return(list())
  candidate_rows <- data.table::rbindlist(lapply(dataset_evidence, function(ds) {
    ids <- ai_record_value(ds$candidates, c("candidateSNP", "rsID", "SNP", "ID"))
    if (!length(ids)) return(NULL)
    data.table::data.table(value = ids, dataset = ds$dataset_label)
  }), fill = TRUE)
  gene_rows <- data.table::rbindlist(lapply(dataset_evidence, function(ds) {
    ids <- ai_record_value(ds$magma_genes, c("SYMBOL", "GENE", "ENTREZID"))
    if (!length(ids)) return(NULL)
    data.table::data.table(value = ids, dataset = ds$dataset_label)
  }), fill = TRUE)
  make_shared <- function(dt, type, prefix) {
    if (!nrow(dt)) return(list())
    shared <- dt[, .(
      n_datasets = data.table::uniqueN(dataset),
      datasets = paste(sort(unique(dataset)), collapse = ", ")
    ), by = value][n_datasets >= 2]
    if (!nrow(shared)) return(list())
    shared[, `:=`(
      comparison_type = type,
      trait = trait,
      evidence_id = sprintf("%s_%03d", prefix, seq_len(.N))
    )]
    data.table::setcolorder(shared, c("evidence_id", "comparison_type", "trait", "value", "n_datasets", "datasets"))
    jsonlite::fromJSON(jsonlite::toJSON(shared, dataframe = "rows", auto_unbox = TRUE), simplifyVector = FALSE)
  }
  c(
    make_shared(candidate_rows, "shared_candidate_variant", "CMP_CAND"),
    make_shared(gene_rows, "shared_magma_gene", "CMP_MAGMA")
  )
}

build_ai_evidence <- function(dataset_information, trait) {
  datasets <- lapply(dataset_information, ai_read_sections, trait = trait)
  list(
    schema_version = "1.0",
    prompt_version = AI_PROMPT_VERSION,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    trait = trait,
    datasets = datasets,
    comparisons = build_ai_comparisons(datasets, trait)
  )
}

write_ai_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  jsonlite::write_json(value, tmp, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  if (!file.rename(tmp, path)) stop("Cannot atomically write ", path)
  invisible(path)
}

submit_ai_job <- function(evidence_file, assessment_file, status_file,
                          partition = Sys.getenv("GWAS_AI_PARTITION", "gpu")) {
  viewer_dir <- Sys.getenv("VIEWER_APP_DIR", unset = getwd())
  runner <- normalizePath(file.path(viewer_dir, "scripts", "run_report_ai.sh"), mustWork = TRUE)
  cache_dir <- dirname(status_file)
  output <- system2(
    "sbatch",
    args = c(
      "--parsable",
      paste0("--partition=", partition),
      paste0("--output=", shQuote(file.path(cache_dir, "slurm_%j.out"))),
      paste0("--error=", shQuote(file.path(cache_dir, "slurm_%j.err"))),
      shQuote(runner), shQuote(evidence_file), shQuote(assessment_file), shQuote(status_file)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(output, "status")
  if (!is.null(exit_status) && exit_status != 0L) stop(paste(output, collapse = "\n"))
  job_id <- sub(";.*$", "", trimws(output[[1]]))
  if (!grepl("^[0-9]+$", job_id)) stop("Unexpected sbatch response: ", paste(output, collapse = " "))
  job_id
}

read_ai_status <- function(cache_dir) {
  path <- file.path(cache_dir, "status.json")
  if (!file.exists(path)) return(list(status = "not_generated"))
  tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) list(status = "unknown", error = e$message))
}

read_ai_assessment <- function(cache_dir) {
  path <- file.path(cache_dir, "assessment.json")
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) NULL)
}

ai_cache_dir_for <- function(dataset_information, trait) {
  if (length(dataset_information) != 1L) stop("AI generation requires exactly one dataset")
  file.path(ai_cache_root(dataset_information[[1]]$base_dir),
            make_ai_cache_key(dataset_information, trait))
}

write_ai_array_manifest <- function(requests) {
  if (!length(requests)) stop("Cannot create an empty AI job-array manifest")
  batch_dir <- file.path(
    ai_cache_root(), "_batches",
    paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  )
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- file.path(batch_dir, "manifest.tsv")
  rows <- data.table::rbindlist(lapply(requests, function(x) data.table::data.table(
    evidence_file = x$evidence_file,
    assessment_file = x$assessment_file,
    status_file = x$status_file
  )))
  data.table::fwrite(rows, manifest, sep = "\t", col.names = FALSE, quote = FALSE)
  manifest
}

submit_ai_job_array <- function(manifest_file, n_tasks,
                                partition = Sys.getenv("GWAS_AI_PARTITION", "gpu")) {
  if (n_tasks < 1L) stop("n_tasks must be positive")
  viewer_dir <- Sys.getenv("VIEWER_APP_DIR", unset = getwd())
  runner <- normalizePath(file.path(viewer_dir, "scripts", "run_report_ai_array.sh"), mustWork = TRUE)
  batch_dir <- dirname(manifest_file)
  max_parallel <- suppressWarnings(as.integer(Sys.getenv("GWAS_AI_ARRAY_PARALLEL", "2")))
  if (!is.finite(max_parallel) || max_parallel < 1L) max_parallel <- 1L
  array_spec <- sprintf("1-%d%%%d", n_tasks, min(max_parallel, n_tasks))
  output <- system2(
    "sbatch",
    args = c(
      "--parsable",
      paste0("--partition=", partition),
      paste0("--array=", array_spec),
      paste0("--output=", shQuote(file.path(batch_dir, "slurm_%A_%a.out"))),
      paste0("--error=", shQuote(file.path(batch_dir, "slurm_%A_%a.err"))),
      shQuote(runner), shQuote(manifest_file)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(output, "status")
  if (!is.null(exit_status) && exit_status != 0L) stop(paste(output, collapse = "\n"))
  job_id <- sub(";.*$", "", trimws(output[[1]]))
  if (!grepl("^[0-9]+$", job_id)) stop("Unexpected sbatch response: ", paste(output, collapse = " "))
  job_id
}
