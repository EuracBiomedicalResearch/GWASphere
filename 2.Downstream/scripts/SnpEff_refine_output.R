library(data.table)
args <- commandArgs(trailingOnly = TRUE)

path <- args[1]
snpeff_dir <- file.path(args[1], "SnpEff")

snpeff_files <- list.files(
  path = snpeff_dir,
  pattern = "\\.ann\\.vcf$",
  full.names = TRUE
)

snpeff_table <- data.table()

for (f in snpeff_files) {
  dt <- fread(f, skip = "CHROM")
  print(f)
  tr <- sub("\\.ann\\.vcf$", "", basename(f))
  if(nrow(dt) > 0){
    dt[, trait := tr]
    dt[, annotation := tstrsplit(INFO, "\\|")[[2]] ]
    dt<-dt[,c("trait", "ID", "annotation")]
    snpeff_table <- rbind(snpeff_table, dt, fill = TRUE)
  }
}
prop_dt <- snpeff_table[, .N, by = c("trait","annotation")]
prop_dt[, proportion := N / sum(N), by= trait]

write.table(
  prop_dt,
  file = file.path(path, "SnpEff", "all_trait_prop.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



