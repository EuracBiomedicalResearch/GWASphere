library(ggplot2)
library(tidyr)
library(dplyr)
library(paletteer)
library(data.table)
library(tidyverse)
library(stringr)
library(bslib)

source("scripts/Manhattan_functions_interactive.R")
args <- commandArgs(trailingOnly = TRUE)

out_dir<- file.path(args[1], "mh_plots")
path<- file.path(args[2], "mplot")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


color<-paletteer_d("ggsci::default_igv")
files <- list.files(path, pattern = "_mplot\\.txt$")

if (length(files) == 0) {
  message("No mplot files found in ", path, " — skipping Manhattan plots")
  quit(status = 0)
} else if (length(files) == 1) {
  # Single trait — no merging needed
  f1  <- fread(file.path(path, files[1]))
  tr1 <- gsub("_mplot.txt", "", files[1])
  merged_files <- f1[, c(4, 2, 3, 5, 6, 9)]
  colnames(merged_files) <- c("SNP", "Chromosome", "Position", "GeneName", "A1FREQ", eval(tr1))
} else {
  # Multiple traits — merge loop
  for(i in 2:length(files)){
    if(i==2){
        f1<-fread(file.path(path,files[i-1]))
        tr1<-gsub("_mplot.txt", "", files[i-1])

        f2<-fread(file.path(path, files[i]))
        tr2<-gsub("_mplot.txt", "", files[i])

        merged_files <- merge(f1[,c(4,2,3,5,6,9)],f2[,c(4,2,3,5,6,9)],by="variant",all=TRUE)
        merged_files$chrom.x <- ifelse(is.na(merged_files$chrom.x),merged_files$chrom.y,merged_files$chrom.x)
        merged_files$pos.x <- ifelse(is.na(merged_files$pos.x),merged_files$pos.y,merged_files$pos.x)
        merged_files$gene_name.x <- ifelse(is.na(merged_files$gene_name.x),merged_files$gene_name.y,merged_files$gene_name.x)
        merged_files$freqA1.x <- ifelse(is.na(merged_files$freqA1.x),merged_files$freqA1.y,merged_files$freqA1.x)
        merged_files <- merged_files[,c("variant","chrom.x","pos.x","gene_name.x","freqA1.x","p.x","p.y")]
        colnames(merged_files) <- c("SNP","Chromosome","Position","GeneName","A1FREQ",eval(tr1),eval(tr2))
    } else{
      fn<-fread(file.path(path,files[i]))
      trn<-gsub("_mplot.txt", "", files[i])
      col<-colnames(merged_files)
      merged_files <- merge(merged_files,fn[,c(4,2,3,5,6,9)],by.x="SNP",by.y="variant",all=TRUE)
      merged_files$Chromosome <- ifelse(is.na(merged_files$Chromosome),merged_files$chrom,merged_files$Chromosome)
      merged_files$Position <- ifelse(is.na(merged_files$Position),merged_files$pos,merged_files$Position)
      merged_files$GeneName <- ifelse(is.na(merged_files$GeneName),merged_files$gene_name,merged_files$GeneName)
      merged_files$A1FREQ <- ifelse(is.na(merged_files$A1FREQ),merged_files$freqA1,merged_files$A1FREQ)
      new_col<-c(col,"p")
      merged_files <- merged_files[,..new_col]
      colnames(merged_files) <- c(col, eval(trn))
    }
  }
}

dlong <- pivot_longer(merged_files, cols = -c("SNP","Chromosome","Position","GeneName","A1FREQ"), names_to = "trait", values_to = "pvalue")
dlong<-na.omit(dlong)
dlong<-as.data.table(dlong)[,LOG10P:= -log10(pvalue)]
setnames(dlong, old = "Chromosome", new = "CHROM")
setnames(dlong, old = "Position", new = "GENPOS")
setnames(dlong, old = "GeneName", new = "GENE_NAME")

df_prep <- manhattan_prep(dlong)

df_plot <- downsample_manhattan(
  df_prep,
  top_n_keep = 2000,
  sig_threshold = 7.301,
  bin_size = 250000
)

p_html <- manhattan_plotly_scattergl(
  df_plot,
  combined = TRUE,
  top_n_hover = 50
)

htmlwidgets::saveWidget(
  p_html,
  file = file.path(out_dir, "combined_mh.html"),
  selfcontained = FALSE
)

for (t in unique(dlong$trait)) {
  dt <- dlong[trait == t]
  df_prep <- manhattan_prep(dt)
  
  df_plot <- downsample_manhattan(
    df_prep,
    top_n_keep = 1000,
    sig_threshold = 7.301,
    bin_size = 250000
  )

  p_html <- manhattan_plotly_scattergl(
    df_plot,
    combined = FALSE,
    top_n_hover = 50
  )

  htmlwidgets::saveWidget(
    p_html,
    file = file.path(out_dir, paste0(t, "_mh.html")),
    selfcontained = FALSE
  )
}


dlong[, MAF := pmin(A1FREQ, 1 - A1FREQ)]
dlong_maf<-dlong[MAF >= 0.01]

  df_prep <- manhattan_prep(dlong_maf)
  df_plot <- downsample_manhattan(
    df_prep,
    top_n_keep = 1000,
    sig_threshold = 7.301,
    bin_size = 250000
  )

p_html <- manhattan_plotly_scattergl(
  df_plot,
  combined = TRUE,
  top_n_hover = 50
)

htmlwidgets::saveWidget(
  p_html,
  file = file.path(out_dir, "combined_mh_maffilt.html"),
  selfcontained = FALSE
)

for (t in unique(dlong_maf$trait)){
    dt<-dlong_maf[trait == t]
    df_prep <- manhattan_prep(dt)
    
    df_plot <- downsample_manhattan(
    df_prep,
    top_n_keep = 1000,
    sig_threshold = 7.301,
    bin_size = 250000
  )

    p_html <- manhattan_plotly_scattergl(
      df_plot,
      combined = FALSE,
      top_n_hover = 50
    )

  htmlwidgets::saveWidget(
    p_html,
    file = file.path(out_dir, paste0(t, "_mh_maffilt.html")),
    selfcontained = FALSE
  )

}
  

