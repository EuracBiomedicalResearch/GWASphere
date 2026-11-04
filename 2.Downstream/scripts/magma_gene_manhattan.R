# magma_gene_manhattan.R
# Generates per-trait interactive gene-based Manhattan plots from MAGMA .genes.symbol.tsv
# Output: mh_plots_gene/{trait}_gene_mh.html

library(data.table)
library(plotly)
library(htmlwidgets)

args     <- commandArgs(trailingOnly = TRUE)
OUT_DIR  <- args[1]
gene_dir <- file.path(OUT_DIR, "magma", "gene_analysis")
out_dir  <- file.path(OUT_DIR, "mh_plots_gene")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(gene_dir, pattern = "\\.genes\\.symbol\\.tsv$", full.names = TRUE)

if (length(files) == 0) {
  message("No .genes.symbol.tsv files found — skipping gene Manhattan plots")
  quit(status = 0)
}

SIG_THRESHOLD <- -log10(2.5e-6)  # Bonferroni for ~20k genes

for (f in files) {
  trait <- sub("\\.genes\\.symbol\\.tsv$", "", basename(f))
  message("Processing: ", trait)

  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) next

  # Expected cols: SYMBOL (or GENE), CHR, START, STOP, P
  if (!all(c("CHR", "START", "STOP", "P") %in% names(dt))) {
    message("  [skip] missing required columns")
    next
  }

  # Gene label: SYMBOL if available, else ENTREZID
  if ("SYMBOL" %in% names(dt)) {
    dt[is.na(SYMBOL) | !nzchar(SYMBOL), SYMBOL := as.character(GENE)]
    dt[, gene_label := SYMBOL]
  } else {
    dt[, gene_label := as.character(GENE)]
  }

  dt[, CHR    := as.integer(sub("^chr", "", as.character(CHR)))]
  dt[, MID    := as.integer((as.numeric(START) + as.numeric(STOP)) / 2)]
  dt[, P      := as.numeric(P)]
  dt[, LOG10P := -log10(pmax(P, .Machine$double.eps))]
  dt[, ZSTAT  := if ("ZSTAT" %in% names(dt)) as.numeric(ZSTAT) else NA_real_]

  dt <- dt[!is.na(CHR) & CHR %in% 1:22 & !is.na(MID) & LOG10P > 2]

  if (nrow(dt) == 0) next

  # Cumulative position
  chr_len <- dt[, .(CHR_len = max(MID)), by = CHR]
  chr_len[, tot := cumsum(as.numeric(CHR_len)) - CHR_len]
  dt <- merge(dt, chr_len[, .(CHR, tot)], by = "CHR")
  dt[, MIDcum := MID + tot]

  axisdf <- dt[, .(center = (max(MIDcum) + min(MIDcum)) / 2), by = CHR]

  # Colour by chromosome (alternating)
  dt[, col_group := ifelse(CHR %% 2 == 0, "even", "odd")]

  # Hover text
  dt[, hover := paste0(
    "Gene: ", gene_label,
    "<br>Chr: ", CHR,
    "<br>-log10(P): ", round(LOG10P, 3),
    if ("ZSTAT" %in% names(dt)) paste0("<br>ZSTAT: ", round(ZSTAT, 3)) else ""
  )]

  # Downsample background
  dt_sig <- dt[LOG10P >= SIG_THRESHOLD]
  dt_bg  <- dt[LOG10P < SIG_THRESHOLD]
  dt_bg  <- dt_bg[order(-LOG10P)][1:min(.N, 2000)]
  dt_plot <- rbind(dt_sig, dt_bg)

  p <- plot_ly() %>%
    add_trace(
      data = dt_plot[col_group == "odd"],
      x = ~MIDcum, y = ~LOG10P, type = "scattergl", mode = "markers",
      marker = list(size = 5, color = "#4878CF", opacity = 0.7),
      text = ~hover, hoverinfo = "text", name = "odd", showlegend = FALSE
    ) %>%
    add_trace(
      data = dt_plot[col_group == "even"],
      x = ~MIDcum, y = ~LOG10P, type = "scattergl", mode = "markers",
      marker = list(size = 5, color = "#6ACC65", opacity = 0.7),
      text = ~hover, hoverinfo = "text", name = "even", showlegend = FALSE
    ) %>%
    add_trace(
      data = dt_sig,
      x = ~MIDcum, y = ~LOG10P, type = "scattergl", mode = "markers+text",
      marker = list(size = 7, color = "#D64550"),
      text = ~gene_label, textposition = "top center",
      textfont = list(size = 10),
      hoverinfo = "text", name = "significant", showlegend = FALSE
    ) %>%
    layout(
      title = list(text = paste0("Gene-based Manhattan — ", trait)),
      xaxis = list(
        title = "Chromosome",
        tickmode = "array",
        tickvals = axisdf$center,
        ticktext = axisdf$CHR,
        showgrid = FALSE, zeroline = FALSE
      ),
      yaxis = list(
        title = "-log10(P) gene-based",
        showgrid = FALSE, zeroline = FALSE
      ),
      shapes = list(list(
        type = "line",
        x0 = min(dt$MIDcum), x1 = max(dt$MIDcum),
        y0 = SIG_THRESHOLD, y1 = SIG_THRESHOLD,
        line = list(color = "firebrick", dash = "dash", width = 1.5)
      )),
      plot_bgcolor  = "white",
      paper_bgcolor = "white"
    )

  out_file <- file.path(out_dir, paste0(trait, "_gene_mh.html"))
  htmlwidgets::saveWidget(p, file = out_file, selfcontained = FALSE)
  message("  Written: ", out_file)
}

message("Gene Manhattan plots completed.")
