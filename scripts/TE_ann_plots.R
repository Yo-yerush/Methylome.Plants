TE_ann_plots <- function(context, TE.gr, sum_dH = F) {
  region_analysis <- ifelse(!sum_dH, "DMRs", "SurpMRs")

  dir.create("TEs_addiotionnal_results", showWarnings = F)

  DMRsReplicates_TE_file.0 <- paste0(context, "/", region_analysis, "_Transposable_Elements_", context, "_genom_annotations.csv")
  DMRsReplicates_Genes_file.0 <- paste0(context, "/", region_analysis, "_Genes_", context, "_genom_annotations.csv")

  if (file.exists(DMRsReplicates_TE_file.0)) {
    DMRsReplicates_TE_file <- read.csv(DMRsReplicates_TE_file.0)
    DMRsReplicates_Genes_file <- read.csv(DMRsReplicates_Genes_file.0)
    DMRsReplicates_Genes_file <- DMRsReplicates_Genes_file[DMRsReplicates_Genes_file$gene_model_type == "protein_coding", ]

    ######################### TE pie plot
    dir.create("TEs_addiotionnal_results/super_family_frequency", showWarnings = F)
    color_vec <- c(brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Dark2"))
    color_vec <- paste0(color_vec, "90")
    col_indx <- data.frame(
      SF_name = sort(unique(TE.gr$Transposon_Super_Family)),
      col = color_vec[1:length(unique(TE.gr$Transposon_Super_Family))]
    )

    TE_Freq <- as.data.frame(table(DMRsReplicates_TE_file$Transposon_Super_Family))
    names(TE_Freq)[1] <- "SF_name"
    TE_Freq <- merge.data.frame(TE_Freq, col_indx, by = "SF_name")

    # display top TE-SuperFamily labels
    top_n <- 6
    labels <- rep("", nrow(TE_Freq))
    indices_top_n <- order(TE_Freq$Freq, decreasing = TRUE)[1:top_n]
    labels[indices_top_n] <- as.character(TE_Freq$SF_name[indices_top_n])

    img_device(paste0("TEs_addiotionnal_results/super_family_frequency/", region_analysis, "_", context, "_TE_super_family"), w = 8, h = 5.5)
    par(mar = c(1, 4, 1, 4))
    par(fig = c(0, 6, 0, 10) / 10)
    par(lwd = 2)
    # pie plot
    pie(TE_Freq$Freq,
      labels = labels,
      main = "",
      border = "white",
      col = TE_Freq$col
    )
    # outer circle
    symbols(0, 0,
      circles = 0.8, inches = FALSE, add = TRUE,
      fg = "#575652"
    )
    # title
    title(
      main = paste0(context, " context"),
      line = -3, cex.main = 1.25
    )

    # legend
    par(new = T)
    par(mar = c(1, 4, 1, 4))
    par(fig = c(6, 10, 0, 10) / 10)
    par(lwd = 1.25)
    legend("center", legend = as.character(TE_Freq$SF_name), fill = TE_Freq$col, bty = "n")

    dev.off()
  }

  ######################### TE vs Genes plot
  # dir.create("TEs_addiotionnal_results/coding_genes_vs_TEs", showWarnings = F)
  GvsT_df <- data.frame(
    x = c("Genes", "TE"),
    y = c(nrow(DMRsReplicates_Genes_file), nrow(DMRsReplicates_TE_file))
  )
  GvsT_plot <- ggplot(GvsT_df, aes(x = x, y = y, fill = x)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_col(color = "gray25") +
    scale_fill_manual(values = c("Genes" = "#80808020", "TE" = "#86671220")) +
    xlab("") +
    ylab(paste("# of", region_analysis)) +
    ggtitle(context) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 11),
          legend.position = "none")
  # img_device(paste0("TEs_addiotionnal_results/coding_genes_vs_TEs/", context,"_coding_genes_vs_TEs"), w = 1.6, h = 2.1)
  # plot(GvsT_plot)
  # dev.off()
  return(GvsT_plot)
}
