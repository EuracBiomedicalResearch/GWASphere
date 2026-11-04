library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
args <- commandArgs(trailingOnly = TRUE)

magma_dir<- file.path(args[1], "magma", "gene_analysis")


files <- list.files(
  magma_dir,
  pattern = "\\.genes\\.out$",
  full.names = TRUE
)

if(length(files) == 0){
  stop("No MAGMA .genes.out files found")
}

for(f in files){

  cat("\n====================================\n")
  cat("Processing:", basename(f), "\n")

  dt <- fread(f)

  if(nrow(dt) == 0){
    cat("  [skip] empty file\n")
    next
  }

  # MAGMA first column = ENTREZ
  gene_col <- names(dt)[1]

  dt[, ENTREZID := as.character(get(gene_col))]

  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(dt$ENTREZID),
    keytype = "ENTREZID",
    columns = c("SYMBOL")
  )

  map <- unique(map)

  dt2 <- merge(
    dt,
    map,
    by = "ENTREZID",
    all.x = TRUE
  )

  # Put SYMBOL first
  cols <- c(
    "SYMBOL",
    setdiff(names(dt2), "SYMBOL")
  )

  dt2 <- dt2[, ..cols]

  out <- sub(
    "\\.genes\\.out$",
    ".genes.symbol.tsv",
    f
  )

  fwrite(
    dt2,
    out,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat("  wrote:", out, "\n")
}

cat("\nDone.\n")