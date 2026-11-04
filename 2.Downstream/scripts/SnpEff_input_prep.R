#Create vcf file as input to SnpEff
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

input_path <- args[1]
out_dir<- args[2]

regenie_files <- list.files(
  path = input_path,
  pattern = "\\.regenie\\.gz$"
)

dir.create(
  file.path(out_dir, "SnpEff", "input"),
  recursive = TRUE,
  showWarnings = FALSE
)

for (i in seq_along(regenie_files)) {

  tr <- sub("\\.regenie\\.gz$", "", regenie_files[i])
  message("Processing: ", tr)

  candidate_file <- file.path(
    out_dir,
    "ld",
    "candidate_snps",
    paste0(tr, ".candidate_snps.tsv")
  )

  if (!file.exists(candidate_file)) {
    message("  [skip] candidate file not found: ", candidate_file)
    next
  }

  # Read header to check available columns
  header_cols <- names(fread(
    cmd = paste("zcat", shQuote(file.path(input_path, regenie_files[i]))),
    nrows = 0
  ))

  available_cols <- c("CHROM", "GENPOS", "REF", "ALT", "ID")
  if ("INFO" %in% header_cols) {
    available_cols <- c(available_cols, "INFO")
  }

  f <- fread(
    cmd = paste(
      "zcat",
      shQuote(file.path(input_path, regenie_files[i]))
    ),
    select = available_cols
  )

  # Add INFO column as NA if not present
  if (!"INFO" %in% names(f)) {
    f[, INFO := NA_character_]
  }

  candidate_dt <- fread(candidate_file)
  candidate_ids <- unique(candidate_dt$candidateSNP)

  dt <- f[ID %in% candidate_ids]

  if (nrow(dt) == 0) {
    message("  [skip] no candidate SNPs found in REGENIE file")
    next
  }

  setnames(dt, "GENPOS", "POS")

  dt[, QUAL := "."]
  dt[, FILTER := "PASS"]

  # INFO in VCF should be key=value, not raw number
  dt[, VCF_INFO := ifelse(
    is.na(INFO) | INFO == "",
    ".",
    paste0("IMPUTEINFO=", INFO)
  )]

  vcf_dt <- unique(dt[, .(
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    INFO = VCF_INFO
  )])

  setorder(vcf_dt, CHROM, POS)

  vcf_file <- file.path(
    out_dir,
    "SnpEff",
    "input",
    paste0(tr, "_snpeff_input.vcf")
  )
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=IMPUTEINFO,Number=1,Type=Float,Description=\"Imputation INFO score\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ),
    con = vcf_file
  )

  fwrite(
    vcf_dt,
    file = vcf_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    append = TRUE
  )

  message("  wrote: ", vcf_file)
}