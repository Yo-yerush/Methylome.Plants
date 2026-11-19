total_meth_levels <- function(rep_var1, rep_var2, var1, var2,
                              heterochromatin_ranges = GRanges()) {
  has_structural_regions <- length(heterochromatin_ranges) > 0L
  if (has_structural_regions) {
    hetero.chromatin_ranges_var1 <- subsetByOverlaps(rep_var1, heterochromatin_ranges)
    hetero.chromatin_ranges_var2 <- subsetByOverlaps(rep_var2, heterochromatin_ranges)
    eu.chromatin_ranges_var1 <- subsetByOverlaps(rep_var1, heterochromatin_ranges, invert = TRUE)
    eu.chromatin_ranges_var2 <- subsetByOverlaps(rep_var2, heterochromatin_ranges, invert = TRUE)
  }

  #############################################################

  total_meth_levels_fun <- function(rep_var1_f, rep_var2_f, var1_f, var2_f, plot_title) {
    ############ edit and add ratio (methylated / total) column
    # var1
    meth_ratio_var1 <- rep_var1_f
    meth_ratio_var1$readsT1 <- (meth_ratio_var1$readsM1 / meth_ratio_var1$readsN1) * 100
    try(
      {
        meth_ratio_var1$readsT2 <- (meth_ratio_var1$readsM2 / meth_ratio_var1$readsN2) * 100
      },
      silent = TRUE
    )
    try(
      {
        meth_ratio_var1$readsT3 <- (meth_ratio_var1$readsM3 / meth_ratio_var1$readsN3) * 100
      },
      silent = TRUE
    )

    # var2
    meth_ratio_var2 <- rep_var2_f
    meth_ratio_var2$readsT1 <- (meth_ratio_var2$readsM1 / meth_ratio_var2$readsN1) * 100
    try(
      {
        meth_ratio_var2$readsT2 <- (meth_ratio_var2$readsM2 / meth_ratio_var2$readsN2) * 100
      },
      silent = TRUE
    )
    try(
      {
        meth_ratio_var2$readsT3 <- (meth_ratio_var2$readsM3 / meth_ratio_var2$readsN3) * 100
      },
      silent = TRUE
    )


    ############ average for each context separately
    var_context <- function(x, cntx.loop) {
      x_total <- x[which(x$context == cntx.loop)]
      if (any(grepl("readsT3", names(mcols(x))))) {
        x_cntx <- c(
          mean(x_total@elementMetadata$readsT1, na.rm = T),
          mean(x_total@elementMetadata$readsT2, na.rm = T),
          mean(x_total@elementMetadata$readsT3, na.rm = T)
        )
      } else if (any(grepl("readsT2", names(mcols(x))))) {
        x_cntx <- c(
          mean(x_total@elementMetadata$readsT1, na.rm = T),
          mean(x_total@elementMetadata$readsT2, na.rm = T)
        )
      } else {
        x_cntx <- mean(x_total@elementMetadata$readsT1, na.rm = T)
      }

      return(x_cntx)
    }

    var1_CG <- var_context(meth_ratio_var1, "CG")
    var1_CHG <- var_context(meth_ratio_var1, "CHG")
    var1_CHH <- var_context(meth_ratio_var1, "CHH")

    var2_CG <- var_context(meth_ratio_var2, "CG")
    var2_CHG <- var_context(meth_ratio_var2, "CHG")
    var2_CHH <- var_context(meth_ratio_var2, "CHH")

    ############ plot
    # meth_plot_df = data.frame(type = rep(c("CG", "CHG", "CHH"),2),
    #                           treatment = c(rep(var1,3), rep(var2,3)),
    #                           levels = c(mean(var1_CG),mean(var1_CHG),mean(var1_CHH),
    #                                      mean(var2_CG),mean(var2_CHG),mean(var2_CHH)),
    #                           SD = c(sd(var1_CG),sd(var1_CHG),sd(var1_CHH),
    #                                  sd(var2_CG),sd(var2_CHG),sd(var2_CHH)))

    meth_plot_df <- data.frame(
      type = rep(
        c("CG", "CHG", "CHH", "CG", "CHG", "CHH"),
        times = c(
          length(var1_CG), length(var1_CHG), length(var1_CHH),
          length(var2_CG), length(var2_CHG), length(var2_CHH)
        )
      ),
      treatment = c(
        rep(var1, sum(lengths(list(var1_CG, var1_CHG, var1_CHH)))),
        rep(var2, sum(lengths(list(var2_CG, var2_CHG, var2_CHH))))
      ),
      levels = c(
        var1_CG, var1_CHG, var1_CHH,
        var2_CG, var2_CHG, var2_CHH
      )
    )

    summary_df <- meth_plot_df %>%
      dplyr::group_by(type, treatment) %>%
      dplyr::summarise(
        mean = mean(levels),
        SD = sd(levels),
        n = dplyr::n(),
        .groups = "drop"
      )

    level_order <- c(var1_f, var2_f)
    max_pos <- which.max(summary_df$mean)
    y_max_plot <- (summary_df$mean[max_pos] + summary_df$SD[max_pos]) * 1.1

    if (nchar(var1) > 6 | nchar(var2) > 6) {
      leg_horiz <- 0.75
    } else {
      leg_horiz <- 0.9
    } # legend position

    g1 <- ggplot(data = summary_df, aes(x = type, y = mean, fill = factor(treatment, level = level_order))) +
      geom_bar(stat = "identity", position = position_dodge(), colour = "black") +

      # geom_point(data = meth_plot_df,
      #          aes(x = type, y = levels, color = factor(treatment,level=level_order)),
      #          position = position_jitterdodge(
      #            jitter.width = 0.01, dodge.width = 0.9),
      #          size = 1.75, inherit.aes = FALSE) +
      geom_point(
        data = meth_plot_df %>%
          dplyr::group_by(type, treatment) %>%
          dplyr::mutate(show_points = ifelse(dplyr::n() > 1, 0.6, 0)) %>%
          dplyr::ungroup(),
        aes(
          x = type, y = levels,
          color = factor(treatment, levels = level_order),
          alpha = show_points
        ),
        position = position_jitterdodge(
          jitter.width = 0.01, dodge.width = 0.9
        ),
        size = 1.75,
        inherit.aes = FALSE
      ) +
      scale_alpha_identity(guide = "none") +

      # scale_fill_manual(values=alpha(c("gray88","gray22"),0.2)) +
      scale_fill_manual(values = alpha(c("gray95", "gray35"), 0.35)) +
      scale_color_manual(values = c("steelblue4", "steelblue4"), guide = "none") +
      geom_errorbar(
        # data = summary_df %>% filter(n > 2),
        aes(ymax = mean + SD, ymin = mean - SD),
        position = position_dodge(width = 0.9),
        width = 0.2
      ) +
      theme_classic() +
      theme(
        plot.margin = unit(c(1, 1, 4, 1), "lines"),
        # title = element_text(size = 9, face="bold"),
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.line.x = element_line(linewidth = 1.1),
        axis.line.y = element_line(linewidth = 1.1),
        axis.ticks = element_line(linewidth = 1.1),
        axis.ticks.length = unit(0.01, "cm"),
        legend.position = c(leg_horiz, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, face = "bold")
      ) +
      labs(y = "5-mC%", title = gsub("_", " ", plot_title)) +
      # geom_jitter(color="steelblue4", size=1, alpha=0.9, width = 0.18) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, y_max_plot))

    img_device(paste0(plot_title, "_", var2_f, "_vs_", var1_f), w = 2.5, h = 3.65)
    plot(g1)
    dev.off()

    write.csv(meth_plot_df, paste0(plot_title, "_", var2_f, "_vs_", var1_f, ".csv"), row.names = F)
  }
  #############################################################

  total_meth_levels_fun(rep_var1, rep_var2, var1, var2, "Whole_Genome") # all genom
  cat(".")

  if (has_structural_regions) {
    total_meth_levels_fun(hetero.chromatin_ranges_var1, hetero.chromatin_ranges_var2, var1, var2, "Heterochromatin_region")
    cat(".")
    total_meth_levels_fun(eu.chromatin_ranges_var1, eu.chromatin_ranges_var2, var1, var2, "Euchromatin_region")
    cat(".")
  }
}
