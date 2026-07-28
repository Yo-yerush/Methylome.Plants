delta_metaplot <- function(list_type, var1, var2, is_geneFeature = FALSE, context_legend = TRUE, gene_feature_binSize = 10) {

  cat(paste0("delta >> calculate within each bin [", var2, " - ", var1,"]\n"))

  # read csv files
  delta_by_cntx <- function(cntx) {
    if (!is_geneFeature) {
      suffixes <- c("up.stream", "gene.body", "down.stream")
    } else {
      suffixes <- c(
        "Promoters.features", "fiveUTR.features", "introns.features",
        "CDS.features", "threeUTR.features"
      )
    }

    var1_cntx <- do.call(rbind, lapply(
      paste0(list_type, "/metaPlot_tables/", var1, ".", cntx, ".", suffixes, ".csv"),
      read.csv
    ))
    var2_cntx <- do.call(rbind, lapply(
      paste0(list_type, "/metaPlot_tables/", var2, ".", cntx, ".", suffixes, ".csv"),
      read.csv
    ))

    var2_cntx$Proportion - var1_cntx$Proportion
  }

  total_bins <- ifelse(is_geneFeature, gene_feature_binSize * 5, 60)
  cntx_bind <- cbind(
    data.frame(pos = 1:total_bins),
    CG = delta_by_cntx("CG"),
    CHG = delta_by_cntx("CHG"),
    CHH = delta_by_cntx("CHH")
  ) %>%
    pivot_longer(cols = c(CG, CHG, CHH), names_to = "context", values_to = "delta")

  max_value <- ifelse(max(cntx_bind$delta) < 0, 0, max(cntx_bind$delta))
  min_value <- ifelse(min(cntx_bind$delta) > 0, 0, min(cntx_bind$delta))
  # mid_value <- max_value / 2
  # mid_value_label <- ifelse(nchar(mid_value) == (nchar(max_value) + 1), paste0("  ", mid_value), mid_value)

  legend_labels <- c("CG", "\nCHG", "\n\nCHH")

  plot_colors <- c("#3d53b4", "#3b8f3e", "#bb4949")

  if (!is_geneFeature) {
    breaks_and_labels <- list(breaks = c(1.35, 20, 40, 59.65), labels = c("  -2kb", "TSS", "TTS", "+2kb   "))
    breaks_vline <- c(20, 40)
  } else {
    breaks_x <- seq(0, total_bins, by = gene_feature_binSize)
    breaks_vline <- breaks_x[-c(1, length(breaks_x))]
    breaks_and_labels <- list(
      breaks = breaks_x[1:length(breaks_x)] + (gene_feature_binSize / 2),
      labels = c("Promoter", "5'UTR", "CDS", "Intron", "3'UTR")
    )
  }

  plot_out <- ggplot(data = cntx_bind, aes(x = pos, y = delta, color = context, group = context)) +
    geom_vline(xintercept = breaks_vline, colour = "gray", linetype = "solid", linewidth = 0.5) +
    geom_hline(yintercept = 0, colour = "gray30", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.4) +
    theme_classic() +
    labs(
      title = list_type,
      x = "",
      y = "Δ methylation"
    ) +
    theme(
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9)
    ) +
    scale_x_continuous(breaks = breaks_and_labels$breaks, labels = breaks_and_labels$labels, expand = expansion(add = c(0, 0))) + # minor_breaks = breaks_x
    scale_y_continuous(
      breaks = c(min_value, 0, max_value), limits = c(min_value, max_value),
      labels = c(
      format(min_value, scientific = T, digits = 2),
      0,
      format(max_value, scientific = T, digits = 2)
      )
    ) +
    scale_color_manual(values = plot_colors) +
    {
      if (context_legend) {
        annotate("text",
          x = 3,
          y = max_value - max_value * 0.01,
          label = legend_labels,
          hjust = 0, vjust = 0.75, size = 2.35,
          color = plot_colors, fontface = "bold"
        )
      }
    }


  img_device(paste0(list_type, "/", list_type, "_delta_metaPlot_", var2, "_vs_", var1), w = ifelse(is_geneFeature, 3.6, 1.95), h = 1.94)
  print(plot_out)
  dev.off()
}
