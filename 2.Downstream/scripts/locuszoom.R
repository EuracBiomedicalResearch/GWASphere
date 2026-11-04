library(data.table)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(patchwork)
library(ggrepel)
library(gggenes)
library(biomaRt)

# =========================================================
# CONFIG
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

out_dir<- file.path(args[1], "LZplot")
cand_path<-args[1]
input_path<-args[2]
GENE_LOC<-args[3]
WINDOW_BP <- 500000

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Load MAGMA gene locations
# =========================================================
genes_dt <- fread(
  GENE_LOC,
  header = FALSE
)

# MAGMA gene.loc format:
# ENTREZID CHR START STOP STRAND

setnames(
  genes_dt,
  c("ENTREZID", "CHROM", "start", "end", "strand", "SYMBOL")
)

genes_dt[, ENTREZID := as.character(ENTREZID)]
genes_dt[, CHROM := as.character(CHROM)]

genes_dt[, gene := fifelse(
  is.na(SYMBOL) | SYMBOL == "",
  ENTREZID,
  SYMBOL
)]

genes_dt <- genes_dt[
  CHROM %in% c(as.character(1:22), "X", "Y")
]

genes_dt <- genes_dt[, .(
  gene,
  CHROM,
  start,
  end,
  strand
)]

# =========================================================
# CUSTOM LOCUSZOOM FUNCTION
# =========================================================
custom_locuszoom_plot <- function(plot_dt,
                                  genes_dt,
                                  lead_id,
                                  lead_label,
                                  trait,
                                  out_html,
                                  out_png = NULL,
                                  window_bp = 500000) {

  lead_row <- plot_dt[ID == lead_id][1]

  if (nrow(lead_row) == 0 || is.na(lead_row$GENPOS)) {
    message("    [skip] lead not found in plot_dt: ", lead_id)
    return(NULL)
  }

  chr <- as.character(lead_row$CHROM)
  pos <- as.numeric(lead_row$GENPOS)

  start_bp <- pos - window_bp
  end_bp <- pos + window_bp

  region <- plot_dt[
    as.character(CHROM) == chr &
      GENPOS >= start_bp &
      GENPOS <= end_bp
  ]


  if (nrow(region) < 2) {
    message("    [skip] too few variants in locus window: ", lead_id)
    return(NULL)
  }

  region[, LOG10P := -log10(P)]
  region[is.na(R2), R2 := 0]
  region[ID == lead_id, R2 := 1]
  region[is.na(RSID) | RSID == "", RSID := ID]

  genes_region <- genes_dt[
    as.character(CHROM) == chr &
      end >= start_bp &
      start <= end_bp
  ]

  if (nrow(genes_region) > 0) {
    genes_region[, strand_label := ifelse(strand %in% c("-", -1), "-", "+")]
    genes_region[, y_gene := seq_len(.N)]
  }

  region[, hover_text := fifelse(
    R2 > 0,
    paste0(
      "RSID: ", RSID,
      "<br>Gene: ", GENE_NAME,
      "<br>Chr: ", CHROM,
      "<br>Position: ", GENPOS,
      "<br>-log10(P): ", round(LOG10P, 3),
      "<br>A1FREQ: ", round(A1FREQ, 3)
    ),
    ""
  )]

  region[, POS_MB := GENPOS / 1e6]
  genes_region[, start_mb := start / 1e6]
  genes_region[, end_mb := end / 1e6]

  p1 <- ggplot(
    region,
    aes(
      x = POS_MB,
      y = LOG10P,
      color = R2,
      text = hover_text
    )
  ) +
    geom_point(size = 1.5, alpha = 0.85) +
    geom_point(
      data = region[ID == lead_id],
      aes(x = POS_MB, y = LOG10P),
      inherit.aes = FALSE,
      color = "red",
      size = 3
    ) +
    ggrepel::geom_label_repel(
      data = region[ID == lead_id],
      aes(x = POS_MB, y = LOG10P, label = lead_label),
      inherit.aes = FALSE,
      color = "red",
      size = 4,
      max.overlaps = Inf
    ) +
    scale_color_gradient(
      low = "grey80",
      high = "#D64550",
      limits = c(0, 1),
      name = "r²"
    ) +
    geom_hline(
      yintercept = -log10(5e-8),
      linetype = "longdash",
      color = "firebrick"
    ) +
    labs(
      title = paste0(trait, " | ", lead_label),
      x = NULL,
      y = "-log(10)(P)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  if (nrow(genes_region) > 0) {
    ## static gene track for PNG
    genes_region[, label := gene]

    p2 <- ggplot(genes_region) +
      geom_segment(
        aes(
          x = start_mb,
          xend = end_mb,
          y = y_gene,
          yend = y_gene
        ),
        color = "navy",
        linewidth = 0.4
      ) +
      geom_text(
        aes(
          x = (start_mb + end_mb) / 2,
          y = y_gene + 0.25,
          label = label
        ),
        size = 3,
        color = "black"
      ) +
      scale_x_continuous(
        limits = c(start_bp / 1e6, end_bp / 1e6),
        name = paste0("Chromosome ", chr, " (Mb)")
      ) +
      theme_classic() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
      )

  } else {

    p2 <- ggplot() +
      annotate(
        "text",
        x = mean(c(start_bp, end_bp)),
        y = 1,
        label = "No genes found in this window"
      ) +
      xlim(start_bp, end_bp) +
      theme_void()
  }

  ## interactive plots
  if (nrow(genes_region) > 0) {

    gene_track_ly <- plot_ly()
    region[, POS_MB := GENPOS / 1e6]
    genes_region[, start_mb := start / 1e6]
    genes_region[, end_mb := end / 1e6]
    genes_region[, label := gene]
    gene_track_ly <- plot_ly()

    for (i in seq_len(nrow(genes_region))) {
      g <- genes_region[i]

      gene_track_ly <- gene_track_ly %>%
        add_segments(
          x = g$start_mb,
          xend = g$end_mb,
          y = g$y_gene,
          yend = g$y_gene,
          line = list(color = "navy", width = 2),
          hoverinfo = "text",
          text = paste0(
            "Gene: ", g$gene,
            "<br>Start: ", round(g$start_mb, 3), " Mb",
            "<br>End: ", round(g$end_mb, 3), " Mb",
            "<br>Strand: ", g$strand
          ),
          showlegend = FALSE
        ) %>%
        add_text(
          x = (g$start_mb + g$end_mb) / 2,
          y = g$y_gene + 0.25,
          text = g$label,
          textfont = list(size = 11, color = "black"),
          hoverinfo = "skip",
          showlegend = FALSE
        )
    }

} else {

  gene_track_ly <- plot_ly()
}

p1ly <- ggplotly(
  p1,
  tooltip = "text"
)

  combined_ly <- subplot(
    p1ly,
    gene_track_ly,
    nrows = 2,
    shareX = TRUE,
    heights = c(0.68, 0.32),
    titleY = TRUE
  ) %>%
    layout(
      yaxis = list(
        title = "-log10(P)",
        zeroline = FALSE
      ),
      yaxis2 = list(
        title = "Genes",
        showticklabels = FALSE,
        zeroline = FALSE
      ),
      xaxis = list(
        title = paste0("Chromosome ", chr, " (Mb)"),
        showticklabels = TRUE,
        ticks = "outside",
        showline = TRUE,
        zeroline = FALSE
      ),
      margin = list(l = 70, r = 80, b = 90, t = 60)
    )


  htmlwidgets::saveWidget(
    combined_ly,
    file = out_html,
    selfcontained = TRUE
  )

  combined <- p1 / p2 + patchwork::plot_layout(heights = c(3, 1.4))

  if (!is.null(out_png)) {
    ggsave(
      filename = out_png,
      plot = combined,
      width = 8,
      height = 7,
      dpi = 150
    )
  }

  return(out_html)
}



regenie_files <- list.files(
  path = input_path,
  pattern = "\\.regenie\\.gz$"
)

for (i in seq_along(regenie_files)) {

  tr <- sub("\\.regenie\\.gz$", "", regenie_files[i])
  message("Processing trait: ", tr)

  cand_file <- file.path(
    cand_path,
    "ld/candidate_snps",
    paste0(tr, ".candidate_snps.tsv")
  )

  if (!file.exists(cand_file)) {
    message("  [skip] candidate SNP file not found: ", cand_file)
    next
  }

  cand <- fread(cand_file)

  if (nrow(cand) == 0) {
    message("  [skip] empty candidate file: ", cand_file)
    next
  }

  f <- fread(
    cmd = paste("zcat", shQuote(file.path(input_path, regenie_files[i]))),
    select = c("ID", "CHROM", "GENPOS", "GENE_NAME", "LOG10P", "RSID","A1FREQ")
  )

  f[, LOG10P := as.numeric(LOG10P)]
  f[, GENPOS := as.numeric(GENPOS)]
  f[, P := 10^(-LOG10P)]
  f <- f[LOG10P > 2.3]

  cand <- cand[!is.na(candidateSNP) & !is.na(bestIndSigSNP)]
  cand[, candidateSNP := as.character(candidateSNP)]
  cand[, bestIndSigSNP := as.character(bestIndSigSNP)]
  cand[, maxR2 := as.numeric(maxR2)]

  for (lead in unique(cand$bestIndSigSNP)) {

    message("  Lead/indSig: ", lead)

    lead_info <- f[ID == lead]

    if (nrow(lead_info) == 0) {
      message("    [skip] lead not found in REGENIE: ", lead)
      next
    }

    rs_id_lead <- lead_info$RSID[1]

    if (is.na(rs_id_lead) || rs_id_lead == "") {
      rs_id_lead <- lead
    }

    ld_for_lead <- cand[
      bestIndSigSNP == lead,
      .(
        ID = candidateSNP,
        R2 = maxR2
      )
    ]

    plot_dt <- merge(
      f,
      ld_for_lead,
      by = "ID",
      all.x = TRUE
    )

    plot_dt[is.na(R2), R2 := 0]
    plot_dt[ID == lead, R2 := 1]
    plot_dt[is.na(RSID) | RSID == "", RSID := ID]

    plot_dt[, GENE_NAME := as.character(GENE_NAME)]
    plot_dt <- unique(plot_dt[, .(
      ID,
      CHROM,
      GENE_NAME,
      GENPOS,
      P,
      RSID,
      A1FREQ,
      R2
    )])

    plot_dt <- plot_dt[
      !is.na(P) &
        P > 0 &
        !is.na(GENPOS) &
        !is.na(CHROM)
    ]

    safe_lead <- gsub("[^A-Za-z0-9_.-]", "_", rs_id_lead)

    out_html <- file.path(
      out_dir, 
      paste0(tr, "_", safe_lead, ".html")
    )

    out_png <- file.path(
      out_dir,
      paste0(tr, "_", safe_lead, ".png")
    )

    custom_locuszoom_plot(
      plot_dt = plot_dt,
      genes_dt = genes_dt,
      lead_id = lead,
      lead_label = rs_id_lead,
      trait = tr,
      out_html = out_html,
      out_png = out_png,
      window_bp = 500000
    )
  }
}