## Multiple tracks Q-Q plot
library(data.table)
library(ramwas)

args <- commandArgs(trailingOnly = TRUE)
input_dir<-args[2]

out_dir<- file.path(args[1], "qqplot")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

files <- list.files(
  path = input_dir,
  pattern = "\\.regenie\\.gz$"
)



for(i in 1:length(files)){
  tr1<-gsub(".regenie.gz", "", files[i])
  #if (!file.exists(paste0("QQplot.",tr1,".jpg"))) { 
  f1 <- fread(
    cmd = paste(
      "zcat",
      shQuote(file.path(input_dir, files[i]))
    )
  )
  #Need to remove pvalue with NA values
  f1<-f1[!is.na(LOG10P)]
  qqPrepare <- qqPlotPrepare(pvalues=f1$LOG10P,ismlog10 = T)
  png(
  file.path(out_dir, paste0(tr1, "_qqplot.png")),
  width = 800,
  height = 800,
  res = 150
)
  par(cex.lab = 1.2,  # axis labels
      cex.axis = 1, # axis numbers  
      mar = c(5, 5, 2, 2))  # bottom, left, top, right
  qqPlotFast(qqPrepare,
        xlab = expression("Expected " * -log[10] * "(" * italic(P) * ")"),
        ylab = expression("Observed " * -log[10] * "(" * italic(P) * ")"))
  title(main = tr1)   
  dev.off()
  #s}

}
