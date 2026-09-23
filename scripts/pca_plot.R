pca_plot <- function(replicates_joints, var1, var2, var1_path, var2_path, cntx = "all_contexts") {
    cntx_title <- ifelse(cntx == "all_contexts", "Total", cntx)
    
    condition <- c(
        rep(var1, length(var1_path)),
        rep(var2, length(var2_path))
    )

    ratios_df <- as.data.frame(replicates_joints)
    if (cntx != "all_contexts") {
        ratios_df <- dplyr::filter(ratios_df, context == cntx)
    }
    ratios_df <- dplyr::select(ratios_df, matches("readsM|readsN"))

    # add and keep ratio columns
    n.col.0 <- ncol(ratios_df)
    i.suffix = 1
    i.col <- 1
    while (i.col < n.col.0) {
        ratios_df[, ncol(ratios_df) + 1] <- ratios_df[, i.col] / ratios_df[, i.col + 1]
        i.col <- i.col + 2
    }

    ratios_df <- ratios_df %>%
        dplyr::select(-matches("readsM|readsN")) %>%
        na.omit()

    # PCA
    pca <- prcomp(t(ratios_df))

    # PCA results
    pca_df <- as.data.frame(pca$x) %>%
        mutate(Genotype = condition) %>%
        dplyr::relocate(Genotype, .before = PC1)

    # PCA table
    write.csv(pca_df, paste0(cntx, "_PCA_table_", var2, "_vs_", var1, ".csv"), row.names = F)

    # PCA plot
    pca_df$col[pca_df$Genotype == var1] = "gray40"
    pca_df$col[pca_df$Genotype == var2] = "#bb5e1b"
    
    pca_p <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
        geom_point(aes(color = Genotype), size = 3.5, alpha = 0.85) +
        xlab(paste0("PC1 (", round(summary(pca)$importance[2, 1], 3) * 100, "%)")) +
        ylab(paste0("PC2 (", round(summary(pca)$importance[2, 2], 3) * 100, "%)")) +
        ggtitle(cntx_title) +
        theme_classic() +
        theme(
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text = element_text(size = 9),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            legend.title = element_blank()
        ) +
        geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
        geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
        scale_color_manual(values = pca_df$col, breaks = pca_df$Genotype)

    img_device(paste0(cntx, "_PCA_plot_", var2, "_vs_", var1), w = 3.25, h = 2.25)
    print(pca_p)
    dev.off()

    cat(".")
}
