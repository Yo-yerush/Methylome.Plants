DMRs_ann_plots <- function(var1, var2, context, dH_analysis = F) {
  ###### prepare plot df
  plot_levels_loop <- c("Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs") # ,"Transposable_Elements")
  ann_plot_final_df <- data.frame(ann = plot_levels_loop, total = NA, gain = NA)

  for (ann.loop in plot_levels_loop) {
    if (!dH_analysis) {
      y_title <- "# of DMRs"
      DMRsReplicates_loop_path <- paste0(context, "/DMRs_", ann.loop, "_", context, "_genom_annotations.csv")
      up_col <- "#d97777"
      down_col <- "#7676d6"
      # up_col_leg = "#d96868"
      # down_col_leg = "#6969db"
    } else {
      y_title <- "# of SurpMRs"
      DMRsReplicates_loop_path <- paste0(context, "/SurpMRs_", ann.loop, "_", context, "_genom_annotations.csv")
      up_col <- "#e2b25d"
      down_col <- "#47cf79"
      # up_col_leg = "#c58741"
      # down_col_leg = "#31a35b"
    }

    if (file.exists(DMRsReplicates_loop_path)) {
      ann_plot_vec <- read.csv(DMRsReplicates_loop_path)$regionType
    } else {
      ann_plot_vec <- data.frame()
    }

    ann_plot_final_df[grep(ann.loop, ann_plot_final_df$ann), 2] <- length(ann_plot_vec)
    ann_plot_final_df[grep(ann.loop, ann_plot_final_df$ann), 3] <- length(ann_plot_vec[ann_plot_vec == "gain"])
  }

  ann_plot_final_df$percent <- round((ann_plot_final_df$gain / ann_plot_final_df$total) * 100, digits = 1)
  # ann_plot_final_df$percent[ann_plot_final_df$percent == 0] = 100
  ann_plot_final_df$percent_loss <- round((100 - ann_plot_final_df$percent), digits = 1)

  plot_levels <- c("Promoters", "CDS", "Introns", "5'UTRs", "3'UTRs") # , "TEG")
  ann_plot_final_df$ann <- plot_levels

  y_lim_total <- ann_plot_final_df$total
  y_lim_max <- max(y_lim_total)
  #########################
  # plot
  # img_device(paste0(context, "_genom_annotations"), w = 2.45, h = 2)
  ann_plot <- ann_plot_final_df %>% ggplot() +
    geom_col(
      aes(x = factor(ann, level = plot_levels), y = total, fill = -percent), # -percent beacouse of the color scale
      colour = "black",
      position = "dodge2"
      # alpha = .5
    ) +
    # scale_fill_gradient(low="mediumblue", high="tomato", limits=c(0,100), breaks=seq(0,100,by=10)) +
    scale_fill_gradient2("Regions type\ndistribution",
      midpoint = -50, low = up_col, mid = "#FFFFFF", high = down_col,
      limits = c(-100, 0), breaks = seq(-100, 0, by = 50),
      labels = c("100% Gain", "50%", "100% Loss")
    ) +
    # guides(fill = guide_colourbar()) +
    guides(fill = guide_legend()) +
    # guides(fill = guide_legend(override.aes = list(alpha = .5))) +
    ylim(0, max(ann_plot_final_df$total) * 1.15) +
    # coord_polar() +
    coord_curvedpolar(start = 0) +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(size = 10, face = "bold"), # , margin = margin(t = 0, r = 20, b = 0, l = 0)),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 6, face = "bold"), # element_blank(),
      title = element_text(size = 8, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      # legend.title=element_text(size = 8, face = "bold"),
      # legend.text = element_text(size = 8, face = "bold"),
      legend.position = "none"
    ) +
    geom_text(aes(x = ann, y = y_lim_max * 0.8, label = total), size = 2.65) +
    labs(
      title = context,
      x = element_blank(),
      y = y_title
    )
  # plot(ann_plot)
  # dev.off()
  return(ann_plot)
}

#########################
### legend
gene_ann_legend <- function(up_col = "#d96868", down_col = "#6969db") {
  legend_df <- data.frame(x = 0.85, y = 31:100, fill_val = 31:100)
  labels_left_df <- data.frame(y = c(30, 80), label = c("Loss", "Gain"))
  labels_right_df <- data.frame(y = c(0, 33, 65, 99), label = c("", "100", "50", "100"))

  legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = fill_val)) +
    geom_tile(width = 0.7, height = 1) +
    scale_fill_gradient2(
      midpoint = 65,
      low = down_col,
      mid = "#FFFFFF",
      high = up_col
    ) +
    geom_text(
      data = labels_left_df, aes(x = 0.4, y = y, label = label),
      hjust = 0, vjust = 0, size = 3, inherit.aes = FALSE, angle = 90, fontface = "bold"
    ) +
    geom_text(
      data = labels_right_df, aes(x = 1.3, y = y, label = label),
      hjust = 0, size = 2.5, inherit.aes = FALSE
    ) +
    geom_rect(aes(xmin = 0.5, xmax = 1.2, ymin = 30.5, ymax = 100.5),
      fill = NA, color = "black", linewidth = 0.15, inherit.aes = FALSE
    ) +
    coord_cartesian(xlim = c(0, 4), clip = "off") +
    labs(title = "Regions type\ndistribution (%)") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0, face = "bold", size = 9.5),
      plot.margin = margin(5, 5, 5, 5)
    )

  return(legend_plot)
}
