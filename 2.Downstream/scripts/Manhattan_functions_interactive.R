#Manhattan plot fucntions
library(dplyr)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(base64enc)
color<-paletteer_d("ggsci::default_igv")

manhattan_prep <- function(df) {
  df_prep <- df %>%
    group_by(CHROM) %>%
    summarise(CHROM_len = max(GENPOS), .groups = "drop") %>%
    mutate(tot = cumsum(as.numeric(CHROM_len)) - CHROM_len) %>%
    select(-CHROM_len) %>%
    left_join(df, ., by = "CHROM") %>%
    arrange(CHROM, GENPOS) %>%
    mutate(GENPOScum = GENPOS + tot)

  return(df_prep)
}

manhattan_plot <- function(df_prep,
                           ymax = 0,
                           hlines = TRUE,
                           top_n_hover = 20,
                           top_n_labels = 20) {

  df_top20 <- df_prep %>%
    slice_max(order_by = LOG10P, n = top_n_labels, with_ties = FALSE) %>%
    group_by(GENE_NAME) %>%
    slice_max(order_by = LOG10P, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    slice_max(order_by = LOG10P, n = top_n_labels, with_ties = FALSE)

  df_prep <- df_prep %>%
    mutate(
      hover_text = paste0(
        "SNP: ", SNP,
        "<br>Gene: ", GENE_NAME,
        "<br>Chr: ", CHROM,
        "<br>Position: ", GENPOS,
        "<br>-log10(P): ", round(LOG10P, 3),
        "<br>A1FREQ: ", round(A1FREQ, 3)
      )
    )

  axisdf <- df_prep %>%
    group_by(CHROM) %>%
    summarise(
      center = (max(GENPOScum) + min(GENPOScum)) / 2,
      .groups = "drop"
    )

  ymin <- 2.3
  ydim <- if (ymax == 0) {
    c(ymin, max(df_prep$LOG10P, na.rm = TRUE) + 0.5)
  } else {
    c(ymin, ymax)
  }

  df_top_hover <- df_prep %>%
    slice_max(order_by = LOG10P, n = top_n_hover, with_ties = FALSE)

  df_other <- df_prep %>%
    anti_join(df_top_hover, by = "SNP")

  plot <- ggplot() +
    geom_point(
      data = df_other,
      aes(
        x = GENPOScum,
        y = LOG10P,
        color = as.factor(CHROM)
      ),
      size = 1
    ) +
    geom_point(
      data = df_top_hover,
      aes(
        x = GENPOScum,
        y = LOG10P,
        color = as.factor(CHROM),
        text = hover_text
      ),
      size = 1.4
    ) +
    scale_color_manual(values = rep(c("#779ECB", "#03254c"), 22)) +
    scale_x_continuous(
      labels = axisdf$CHROM,
      breaks = axisdf$center,
      name = "Chromosome",
      expand = c(0.01, 0.01),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.01)),
      name = "-log10(p-value)",
      limits = ydim
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14),
      axis.ticks = element_line(color = "black")
    )

  if (hlines) {
    plot <- plot +
      geom_hline(
        yintercept = 7.301,
        linetype = "longdash",
        color = "firebrick"
      )
  }

  plot <- plot +
    geom_label(
      data = df_top20,
      aes(
        x = GENPOScum,
        y = LOG10P,
        label = GENE_NAME
      ),
      inherit.aes = FALSE,
      fontface = "italic",
      size = 3.5,
      label.padding = unit(0.15, "lines"),
      label.size = 0.25
    )

  return(plot)
}

combo_manhattan_plot <- function(df_prep,
                                 ymax = 0,
                                 hlines = TRUE,
                                 top_n_hover = 20,
                                 top_n_labels = 20) {

  df_top20 <- df_prep %>%
    slice_max(order_by = LOG10P, n = top_n_labels, with_ties = FALSE) %>%
    group_by(GENE_NAME) %>%
    slice_max(order_by = LOG10P, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    slice_max(order_by = LOG10P, n = top_n_labels, with_ties = FALSE)

  df_prep <- df_prep %>%
    mutate(
      hover_text = paste0(
        "Trait: ", trait,
        "<br>SNP: ", SNP,
        "<br>Gene: ", GENE_NAME,
        "<br>Chr: ", CHROM,
        "<br>Position: ", GENPOS,
        "<br>-log10(P): ", round(LOG10P, 3),
        "<br>A1FREQ: ", round(A1FREQ, 3)
      )
    )

  axisdf <- df_prep %>%
    group_by(CHROM) %>%
    summarise(
      center = (max(GENPOScum) + min(GENPOScum)) / 2,
      .groups = "drop"
    )

  ymin <- 2.3
  ydim <- if (ymax == 0) {
    c(ymin, max(df_prep$LOG10P, na.rm = TRUE) + 0.5)
  } else {
    c(ymin, ymax)
  }

  df_rand <- df_prep %>%
    slice_sample(prop = 1)

  df_top_hover <- df_rand %>%
    slice_max(order_by = LOG10P, n = top_n_hover, with_ties = FALSE)

  df_other <- df_rand %>%
    anti_join(df_top_hover, by = "SNP")

  plot <- ggplot() +
    geom_point(
      data = df_other,
      aes(
        x = GENPOScum,
        y = LOG10P,
        color = as.factor(trait)
      ),
      size = 1,
      alpha = 0.7
    ) +
    geom_point(
      data = df_top_hover,
      aes(
        x = GENPOScum,
        y = LOG10P,
        color = as.factor(trait),
        text = hover_text
      ),
      size = 1.4,
      alpha = 0.95
    ) +
    scale_color_manual(values = as.character(color)) +
    scale_x_continuous(
      labels = axisdf$CHROM,
      breaks = axisdf$center,
      name = "Chromosome",
      expand = c(0.01, 0.01),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.01)),
      name = "-log10(p-value)",
      limits = ydim
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14),
      axis.ticks = element_line(color = "black")
    )

  if (hlines) {
    plot <- plot +
      geom_hline(
        yintercept = 7.301,
        linetype = "longdash",
        color = "firebrick"
      )
  }

  plot <- plot +
    geom_label(
      data = df_top20,
      aes(
        x = GENPOScum,
        y = LOG10P,
        label = GENE_NAME
      ),
      inherit.aes = FALSE,
      fontface = "italic",
      size = 3.5,
      label.padding = unit(0.15, "lines"),
      label.size = 0.25
    )

  return(plot)
}


manhattan_plotly_scattergl <- function(df_prep,
                                       combined = FALSE,
                                       top_n_hover = 50,
                                       ymax = 0,
                                       hlines = TRUE,
                                       sig_threshold = 7.301) {

  axisdf <- df_prep %>%
    group_by(CHROM) %>%
    summarise(
      center = (max(GENPOScum, na.rm = TRUE) + min(GENPOScum, na.rm = TRUE)) / 2,
      .groups = "drop"
    )

  ymin <- 2.3
  y_max <- if (ymax == 0) max(df_prep$LOG10P, na.rm = TRUE) + 0.5 else ymax

  df_prep <- df_prep %>%
    mutate(
      is_significant = LOG10P >= sig_threshold,
      hover_text = paste0(
        if (combined) paste0("Trait: ", trait, "<br>") else "",
        "SNP: ", SNP,
        "<br>Gene: ", GENE_NAME,
        "<br>Chr: ", CHROM,
        "<br>Position: ", GENPOS,
        "<br>-log10(P): ", round(LOG10P, 3),
        "<br>A1FREQ: ", round(A1FREQ, 6)
      )
    )

  df_top_hover <- df_prep %>%
    slice_max(order_by = LOG10P, n = top_n_hover, with_ties = FALSE)

  df_other <- df_prep %>%
    anti_join(df_top_hover, by = "SNP")

  df_other_nonsig <- df_other %>% filter(!is_significant)
  df_other_sig <- df_other %>% filter(is_significant)

  color_formula <- if (combined) ~as.factor(trait) else ~as.factor(CHROM)

  palette <- as.character(color)


  p <- plotly::plot_ly() %>%
    plotly::add_trace(
      data = df_other_nonsig,
      x = ~GENPOScum,
      y = ~LOG10P,
      type = "scattergl",
      mode = "markers",
      color = color_formula,
      colors = palette,
      marker = list(size = 3, opacity = 0.55),
      hoverinfo = "skip",
      showlegend = combined
    ) %>%
    plotly::add_trace(
      data = df_other_sig,
      x = ~GENPOScum,
      y = ~LOG10P,
      type = "scattergl",
      mode = "markers",
      color = color_formula,
      colors = palette,
      marker = list(size = 6, opacity = 0.9),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    plotly::add_trace(
      data = df_top_hover,
      x = ~GENPOScum,
      y = ~LOG10P,
      type = "scattergl",
      mode = "markers",
      color = color_formula,
      colors = palette,
      marker = list(size = 7, opacity = 1),
      text = ~hover_text,
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    plotly::layout(
      hovermode = "closest",
      plot_bgcolor = "white",
      paper_bgcolor = "white",
      legend = list(
        itemsizing = "constant",
        font = list(size = 12)),
      xaxis = list(
        title = "Chromosome",
        tickmode = "array",
        tickvals = axisdf$center,
        ticktext = axisdf$CHROM,
        tickangle = 0,
        showgrid = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "-log10(p-value)",
        range = c(ymin, y_max),
        showgrid = FALSE,
        zeroline = FALSE
      )
    )

 if (hlines) {

  p <- p %>%
    plotly::layout(
      shapes = c(
        p$x$layout$shapes,

        list(
          type = "line",
          x0 = min(df_prep$GENPOScum, na.rm = TRUE),
          x1 = max(df_prep$GENPOScum, na.rm = TRUE),
          y0 = sig_threshold,
          y1 = sig_threshold,
          line = list(
            color = "firebrick",
            dash = "dash",
            width = 1.5
          )
        )
      )
    )
}
  return(p)
}

downsample_manhattan <- function(df_prep,
                                 top_n_keep = 2000,
                                 sig_threshold = 7.301,
                                 bin_size = 250000) {

  df_top <- df_prep %>%
    dplyr::filter(LOG10P >= sig_threshold) %>%
    dplyr::bind_rows(
      df_prep %>%
        dplyr::slice_max(order_by = LOG10P, n = top_n_keep, with_ties = FALSE)
    ) %>%
    dplyr::distinct(SNP, trait, .keep_all = TRUE)

  df_background <- df_prep %>%
    dplyr::filter(!paste(SNP, trait) %in% paste(df_top$SNP, df_top$trait)) %>%
    dplyr::mutate(bin = floor(GENPOScum / bin_size)) %>%
    dplyr::group_by(CHROM, trait, bin) %>%
    dplyr::slice_max(order_by = LOG10P, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-bin)

  dplyr::bind_rows(df_background, df_top) %>%
    dplyr::distinct(SNP, trait, .keep_all = TRUE)
}