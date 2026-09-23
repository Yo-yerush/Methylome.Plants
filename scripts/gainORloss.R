gainORloss <- function(DMRsReplicates, context, comparison_name, up_col = "#d96c6c", down_col = "#6c96d9", add_count = FALSE) {
  gainORloss_vec <- DMRsReplicates$regionType

  gain_DMRs <- length(grep("gain", gainORloss_vec))
  loss_DMRs <- length(grep("loss", gainORloss_vec))
  total <- gain_DMRs + loss_DMRs
  pres_gain <- round((gain_DMRs / total) * 100, 1)
  pres_loss <- round((loss_DMRs / total) * 100, 1)

  count_gain <- ifelse(add_count, paste0(" (", gain_DMRs, ")"), "")
  count_loss <- ifelse(add_count, paste0(" (", loss_DMRs, ")"), "")

  pie_data <- data.frame(
    group = c("gain", "loss"),
    text = c(paste0(pres_gain, "%", count_gain), paste0(pres_loss, "%", count_loss)),
    value = c(gain_DMRs, loss_DMRs)
  )

  # pie plot
  pie_plot <- ggplot(pie_data, aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white", lwd = 0.5) +
    coord_polar(theta = "y", start = pi / 2, direction = 1) +
    scale_fill_manual(values = c(gain = up_col, loss = down_col)) +
    theme_void() +
    theme(
      axis.line = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      text = element_text(family = "serif")
    )

  # text of the plot from the rught side
  text_plot <- ggplot() +
    geom_text(aes(x = -0.25, y = 1), label = pie_data$text[1], size = 4) +
    geom_text(aes(x = -0.25, y = -1), label = pie_data$text[2], size = 4) +
    xlim(-1, 1) +
    ylim(-2, 2) +
    theme_void() +
    theme(text = element_text(family = "serif"))

  comb_plot <- cowplot::plot_grid(pie_plot, text_plot,
    nrow = 1,
    rel_widths = c(2, 2)
  )

  if (add_count) {
    x_pos <- c(0.375, 0.51)
  } else {
    x_pos <- c(0.410, 0.575)
  }
  final_plot <- ggdraw(comb_plot) +
    # hyper line
    draw_line(
      x = x_pos,
      y = c(0.725, 0.725),
      color = "black",
      size = 0.25
    ) +
    # hypo line
    draw_line(
      x = x_pos,
      y = c(0.275, 0.275),
      color = "black",
      size = 0.25
    )

  img_device(paste0("pie_", context, "_gainORloss_", comparison_name), w = ifelse(add_count, 2.5, 2), h = 1)
  par(mar = rep(0, 4)) # bottom, left, top, right
  print(final_plot)
  dev.off()
}

########################################################

ratio.distribution <- function(DMRsReplicates, context, comparison_name, up_col = "#b36b74", down_col = "#9396c2") {
  dmrs_ratio <- DMRsReplicates$proportion2 / DMRsReplicates$proportion1
  n_breaks <- if (length(dmrs_ratio) >= 1800) {
     150
  } else if (length(dmrs_ratio) <= 100) {
     length(dmrs_ratio) / 3
  } else {
     length(dmrs_ratio) / 12
  }

  img_device(paste0("ratio.distribution_", context, "_gainORloss"), w = 3.58, h = 3.3)

  h <- hist(log2(dmrs_ratio), breaks = n_breaks, plot = FALSE)
  hist(log2(dmrs_ratio),
    main = paste(context, "-", length(DMRsReplicates$regionType), "DMRs"),
    xlab = "log2(fold-change)",
    ylab = "Frequency",
    # xlim = c(-10,10),
    # col = c(rep("blue",length(h$counts)/2),rep("red",length(h$counts)/2)),
    # border = c(rep("blue",length(h$counts)/2),rep("red",length(h$counts)/2)),
    border = ifelse((h$mids > 0), up_col, down_col),
    breaks = n_breaks
  )

  dev.off()
}
