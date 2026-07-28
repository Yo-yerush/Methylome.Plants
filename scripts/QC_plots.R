###########################################################################
# QC_plots.R
# Sample-level quality control plots for methylation data:
#   1. Per-region methylation distributions (gene body, promoter, TE)
#   2. Sample-vs-sample scatter plots on genomic tiles
#   3. Sample correlation heatmap
#   4. Variance/dispersion test (treatment vs control)
###########################################################################


###########################################################################
# Helper: per-tile methylation matrix from per-condition joined replicates
###########################################################################
qc_tile_methylation <- function(meth_var1_rep, meth_var2_rep,
                                var1, var2, var1_path, var2_path,
                                tile_width = 100000, minReadsPerTile = 6000,
                                context = "CG") {

  n_var1 <- length(var1_path)
  n_var2 <- length(var2_path)

  # TAIR10 chromosome lengths
  chr_lengths <- tapply(end(meth_var1_rep), seqnames(meth_var1_rep), max)
  genome_gr <- GRanges(seqnames = names(chr_lengths),
                       ranges = IRanges(start = 1, end = chr_lengths))

  # adjust to available seqlevels
  avail_chr <- intersect(names(chr_lengths), seqlevels(meth_var1_rep))
  genome_gr <- genome_gr[seqnames(genome_gr) %in% avail_chr]

  tiles <- unlist(tile(genome_gr, width = tile_width))

  # compute per-tile weighted methylation for each replicate
  tile_ratios <- function(meth_rep, n_reps, cntx) {
    meth_cx <- meth_rep[meth_rep$context == cntx]
    hits <- findOverlaps(meth_cx, tiles)
    q_idx <- queryHits(hits)
    s_idx <- subjectHits(hits)

    result_list <- list()
    for (i in seq_len(n_reps)) {
      mCol <- paste0("readsM", i)
      nCol <- paste0("readsN", i)
      dt <- data.table(
        tile_idx = s_idx,
        rM = mcols(meth_cx)[[mCol]][q_idx],
        rN = mcols(meth_cx)[[nCol]][q_idx]
      )
      tile_stats <- dt[, .(sumM = sum(rM, na.rm = TRUE),
                           sumN = sum(rN, na.rm = TRUE)), by = tile_idx]
      tile_stats <- tile_stats[sumN >= minReadsPerTile]
      tile_stats[, ratio := sumM / sumN * 100]
      result_list[[i]] <- tile_stats[, .(tile_idx, ratio)]
      setnames(result_list[[i]], "ratio", paste0("ratio_", i))
    }

    merged <- result_list[[1]]
    if (n_reps > 1) {
      for (i in 2:n_reps) {
        merged <- merge(merged, result_list[[i]], by = "tile_idx", all = FALSE)
      }
    }
    return(merged)
  }

  var1_tiles <- tile_ratios(meth_var1_rep, n_var1, context)
  var2_tiles <- tile_ratios(meth_var2_rep, n_var2, context)

  # sample labels
  var1_sample_names <- if (n_var1 == 1) var1 else paste0(var1, "_rep", seq_len(n_var1))
  var2_sample_names <- if (n_var2 == 1) var2 else paste0(var2, "_rep", seq_len(n_var2))

  setnames(var1_tiles, paste0("ratio_", seq_len(n_var1)), var1_sample_names)
  setnames(var2_tiles, paste0("ratio_", seq_len(n_var2)), var2_sample_names)

  merged_all <- merge(var1_tiles, var2_tiles, by = "tile_idx", all = FALSE)

  sample_names <- c(var1_sample_names, var2_sample_names)
  condition_labels <- c(rep(var1, n_var1), rep(var2, n_var2))

  return(list(
    tiles_df = merged_all,
    sample_names = sample_names,
    condition_labels = condition_labels,
    n_var1 = n_var1,
    n_var2 = n_var2
  ))
}


###########################################################################
# 1. Per-region CG/CHG/CHH methylation distribution (violin + boxplot)
#    for coding-gene bodies, promoters, TEs and gbM list
###########################################################################
qc_meth_distribution <- function(meth_var1, meth_var2, annotation.gr, TE_gr,
                                 var1, var2, Methylome.At_path = ".", minReads = 6) {

  Genes_0 <- annotation.gr[annotation.gr$type == "gene"]
  Genes <- Genes_0[grepl("protein_coding", Genes_0$gene_model_type)]
  Promoters <- promoters(Genes, upstream = 2000, downstream = 0, use.names = TRUE)

  # load stable/dynamic gbM gene lists (Williams et al. 2023)
  stable_gbM_file <- file.path(Methylome.At_path, "annotation_files", "At_stable_gbM_Williams23.txt")
  dynamic_gbM_file <- file.path(Methylome.At_path, "annotation_files", "At_dynamic_gbM_Williams23.txt")

  features <- list(
    list(name = "Gene body", gr = Genes),
    list(name = "Promoter", gr = Promoters),
    list(name = "TE", gr = TE_gr)
  )
  feature_levels <- c("Gene body", "Promoter", "TE")

  if (file.exists(stable_gbM_file)) {
    stable_ids <- readLines(stable_gbM_file)
    stable_ids <- trimws(stable_ids[nchar(stable_ids) > 0])
    stable_gr <- Genes_0[Genes_0$gene_id %in% stable_ids]
    if (length(stable_gr) > 0) {
      features[[length(features) + 1]] <- list(name = "Stable gbM", gr = stable_gr)
      feature_levels <- c(feature_levels, "Stable gbM")
    }
  }

  if (file.exists(dynamic_gbM_file)) {
    dynamic_ids <- readLines(dynamic_gbM_file)
    dynamic_ids <- trimws(dynamic_ids[nchar(dynamic_ids) > 0])
    dynamic_gr <- Genes_0[Genes_0$gene_id %in% dynamic_ids]
    if (length(dynamic_gr) > 0) {
      features[[length(features) + 1]] <- list(name = "Dynamic gbM", gr = dynamic_gr)
      feature_levels <- c(feature_levels, "Dynamic gbM")
    }
  }

  # weighted methylation per genomic region
  region_meth <- function(meth_data, regions, context) {
    meth_cx <- meth_data[meth_data$context == context & meth_data$readsN >= minReads]
    hits <- findOverlaps(meth_cx, regions)
    if (length(hits) == 0) return(numeric(0))
    dt <- data.table(
      region_idx = subjectHits(hits),
      readsM = meth_cx$readsM[queryHits(hits)],
      readsN = meth_cx$readsN[queryHits(hits)]
    )
    region_stats <- dt[, .(meth_pct = sum(readsM) / sum(readsN) * 100), by = region_idx]
    return(region_stats$meth_pct)
  }

  plot_data_list <- list()
  for (cntx in c("CG", "CHG", "CHH")) {
    for (feat in features) {
      vals1 <- region_meth(meth_var1, feat$gr, cntx)
      vals2 <- region_meth(meth_var2, feat$gr, cntx)
      if (length(vals1) > 0 | length(vals2) > 0) {
        plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
          methylation = c(vals1, vals2),
          condition = c(rep(var1, length(vals1)), rep(var2, length(vals2))),
          context = cntx,
          feature = feat$name,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  plot_df <- do.call(rbind, plot_data_list)
  plot_df$condition <- factor(plot_df$condition, levels = c(var1, var2))
  plot_df$context <- factor(plot_df$context, levels = c("CG", "CHG", "CHH"))
  plot_df$feature <- factor(plot_df$feature, levels = feature_levels)

  n_feat <- length(feature_levels)
  p <- ggplot(plot_df, aes(x = condition, y = methylation, fill = condition)) +
    geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
    facet_grid(feature ~ context, scales = "free_y") +
    scale_fill_manual(values = c("gray50", "#bb5e1b")) +
    labs(y = "5-mC%", x = NULL) +
    theme_classic() +
    theme(
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )

  img_device(paste0("QC_meth_distribution_", var2, "_vs_", var1), w = 7, h = 6)
  print(p)
  dev.off()

  # summary table
  summary_df <- plot_df %>%
    group_by(feature, context, condition) %>%
    summarise(
      n_regions = dplyr::n(),
      median_meth = round(median(methylation, na.rm = TRUE), 2),
      mean_meth = round(mean(methylation, na.rm = TRUE), 2),
      sd_meth = round(sd(methylation, na.rm = TRUE), 2),
      .groups = "drop"
    )
  write.csv(summary_df,
            paste0("QC_meth_distribution_summary_", var2, "_vs_", var1, ".csv"),
            row.names = FALSE)
}


###########################################################################
# 2. Sample-vs-sample scatter plots (genomic tiles)
###########################################################################
qc_sample_scatter <- function(tile_data, var1, var2, context) {

  tiles_df <- tile_data$tiles_df
  sample_names <- tile_data$sample_names
  condition_labels <- tile_data$condition_labels
  n_samples <- length(sample_names)

  if (n_samples < 2) return(invisible(NULL))

  # condition colors (same as distribution plot)
  cond_colors <- c(setNames("gray50", var1), setNames("#bb5e1b", var2))

  pairs <- combn(seq_len(n_samples), 2, simplify = FALSE)

  for (pr in pairs) {
    s1 <- sample_names[pr[1]]
    s2 <- sample_names[pr[2]]
    col_x <- cond_colors[condition_labels[pr[1]]]
    col_y <- cond_colors[condition_labels[pr[2]]]

    sub_df <- data.frame(
      x = tiles_df[[s1]],
      y = tiles_df[[s2]]
    )
    sub_df <- sub_df[complete.cases(sub_df), ]
    if (nrow(sub_df) < 10) next

    r_val <- round(cor(sub_df$x, sub_df$y, use = "complete.obs"), 3)

    p <- ggplot(sub_df, aes(x = x, y = y)) +
      geom_bin2d(bins = 80) +
      scale_fill_viridis_c(option = "inferno", trans = "log10") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
      labs(
        x = paste0(gsub("_", " ", s1), " (5-mC%)"),
        y = paste0(gsub("_", " ", s2), " (5-mC%)"),
        title = paste0("r = ", r_val)
      ) +
      coord_fixed() +
      theme_classic() +
      theme(
        axis.title.x = element_text(size = 10, color = col_x, face = "bold"),
        axis.title.y = element_text(size = 10, color = col_y, face = "bold"),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7)
      )

    fname <- paste0("QC_scatter_", context, "_", gsub(" ", "_", s2), "_vs_", gsub(" ", "_", s1))
    img_device(fname, w = 3.25, h = 3)
    print(p)
    dev.off()
  }

  return(invisible(NULL))
}


###########################################################################
# 3. Sample correlation heatmap
###########################################################################
qc_correlation_heatmap <- function(tile_data, var1, var2, context) {

  tiles_df <- tile_data$tiles_df
  sample_names <- tile_data$sample_names

  if (length(sample_names) < 2) return(invisible(NULL))

  mat <- as.matrix(tiles_df[, ..sample_names])
  cor_mat <- cor(mat, use = "pairwise.complete.obs")

  write.csv(round(cor_mat, 4),
            paste0("QC_correlation_matrix_", context, "_", var2, "_vs_", var1, ".csv"))

  # melt correlation matrix for ggplot (no reshape2 dependency)
  cor_long <- data.frame(
    Var1 = rep(rownames(cor_mat), each = ncol(cor_mat)),
    Var2 = rep(colnames(cor_mat), times = nrow(cor_mat)),
    value = as.vector(t(cor_mat)),
    stringsAsFactors = FALSE
  )
  cor_long$Var1 <- factor(cor_long$Var1, levels = sample_names)
  cor_long$Var2 <- factor(cor_long$Var2, levels = rev(sample_names))

  p <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 3)), size = 3) +
    scale_fill_gradient(
      low = "#ffffdba8", high = "#5e0059a8",
      limits = c(min(cor_long$value), 1),
      name = "Pearson r"
    ) +
    labs(title = paste0(context, " sample correlation"), x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    ) +
    coord_fixed()

  n_s <- length(sample_names)
  hw <- max(2.35, n_s * 0.75 + 1)

  img_device(paste0("QC_correlation_heatmap_", context, "_", var2, "_vs_", var1),
             w = hw, h = hw)
  print(p)
  dev.off()
}


###########################################################################
# 4. Variance / dispersion test (treatment vs control per tile)
###########################################################################
qc_variance_test <- function(tile_data, var1, var2, context) {

  n_var1 <- tile_data$n_var1
  n_var2 <- tile_data$n_var2

  # need >= 2 replicates per condition
  if (n_var1 < 2 | n_var2 < 2) {
    cat(" (skipped - requires >= 2 replicates per condition)")
    return(invisible(NULL))
  }

  tiles_df <- tile_data$tiles_df
  sample_names <- tile_data$sample_names
  condition_labels <- tile_data$condition_labels

  var1_cols <- sample_names[condition_labels == var1]
  var2_cols <- sample_names[condition_labels == var2]

  var1_mat <- as.matrix(tiles_df[, ..var1_cols])
  var2_mat <- as.matrix(tiles_df[, ..var2_cols])

  var1_variance <- apply(var1_mat, 1, var, na.rm = TRUE)
  var2_variance <- apply(var2_mat, 1, var, na.rm = TRUE)

  keep <- !is.na(var1_variance) & !is.na(var2_variance) &
    var1_variance > 0 & var2_variance > 0
  var1_v <- var1_variance[keep]
  var2_v <- var2_variance[keep]

  if (length(var1_v) < 10) {
    cat(" (skipped - too few informative tiles)")
    return(invisible(NULL))
  }

  f_ratio <- mean(var2_v) / mean(var1_v)

  # paired Wilcoxon: is treatment variance systematically higher?
  wt <- wilcox.test(var2_v, var1_v, paired = TRUE, alternative = "greater")

  result_df <- data.frame(
    context = context,
    n_tiles = sum(keep),
    mean_var_control = round(mean(var1_v), 4),
    mean_var_treatment = round(mean(var2_v), 4),
    median_var_control = round(median(var1_v), 4),
    median_var_treatment = round(median(var2_v), 4),
    variance_ratio = round(f_ratio, 4),
    wilcox_p = signif(wt$p.value, 4),
    stringsAsFactors = FALSE
  )

  # scatter of per-tile variances (log10 scale)
  var_df <- data.frame(ctrl = var1_v, trt = var2_v)

  p <- ggplot(var_df, aes(x = log10(ctrl + 1), y = log10(trt + 1))) +
    geom_bin2d(bins = 60) +
    scale_fill_viridis_c(option = "mako", trans = "log10") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red3") +
    labs(
      x = paste0("log10(", var1, " var + 1)"),
      y = paste0("log10(", var2, " var + 1)"),
      title = paste0(context, "   F-ratio=", round(f_ratio, 2),
                     "   Wilcox p=", signif(wt$p.value, 3))
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 9, hjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
    )

  img_device(paste0("QC_variance_", context, "_", var2, "_vs_", var1), w = 4.15, h = 3.6)
  print(p)
  dev.off()

  return(result_df)
}


###########################################################################
# Main QC runner
###########################################################################
run_QC_plots <- function(meth_var1, meth_var2,
                         meth_var1_replicates, meth_var2_replicates,
                         var1, var2, var1_path, var2_path,
                         annotation.gr, TE_gr,
                         Methylome.At_path = ".",
                         tile_width = 1500, minReadsPerTile = 10) {

  # create subdirectories
  qc_base <- getwd()
  dist_dir <- file.path(qc_base, "distribution")
  scatter_dir <- file.path(qc_base, "scatter")
  corr_dir <- file.path(qc_base, "correlation")
  var_dir <- file.path(qc_base, "variance")
  for (d in c(dist_dir, scatter_dir, corr_dir, var_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  ##### 1. methylation distribution plots #####
  cat("\nmethylation distribution (gene body / promoter / TE)...")
  message(time_msg(), "QC: methylation distribution plots: ", appendLF = FALSE)
  tryCatch(
    {
      setwd(dist_dir)
      qc_meth_distribution(meth_var1, meth_var2, annotation.gr, TE_gr, var1, var2,
                            Methylome.At_path = Methylome.At_path)
      setwd(qc_base)
      message("done")
      cat(" done\n")
    },
    error = function(cond) {
      setwd(qc_base)
      cat("\n*\n QC methylation distribution:\n", as.character(cond), "*\n")
      message("fail")
    }
  )

  ##### 2-4. tile-based QC (scatter, correlation, variance) per context #####
  variance_results <- list()
  for (cntx in c("CG", "CHG", "CHH")) {
    cat(paste0("\n", cntx, " tiled QC:"))
    message(time_msg(), cntx, " tiled QC (scatter / correlation / variance): ", appendLF = FALSE)
    tryCatch(
      {
        tile_data <- qc_tile_methylation(
          meth_var1_replicates, meth_var2_replicates,
          var1, var2, var1_path, var2_path,
          tile_width = tile_width,
          minReadsPerTile = minReadsPerTile,
          context = cntx
        )

        # scatter plots (individual files)
        setwd(scatter_dir)
        qc_sample_scatter(tile_data, var1, var2, cntx)
        setwd(qc_base)
        cat(" scatter.")

        # correlation heatmap
        setwd(corr_dir)
        qc_correlation_heatmap(tile_data, var1, var2, cntx)
        setwd(qc_base)
        cat(" correlation.")

        # variance test
        setwd(var_dir)
        var_result <- qc_variance_test(tile_data, var1, var2, cntx)
        setwd(qc_base)
        if (!is.null(var_result)) variance_results[[cntx]] <- var_result
        cat(" variance.")

        message("done")
        cat(" done")
      },
      error = function(cond) {
        setwd(qc_base)
        cat("\n*\n ", cntx, " tiled QC:\n", as.character(cond), "*\n")
        message("fail")
      }
    )
  }

  # write combined variance test results
  if (length(variance_results) > 0) {
    setwd(var_dir)
    var_results_df <- do.call(rbind, variance_results)
    write.csv(var_results_df,
              paste0("QC_variance_test_", var2, "_vs_", var1, ".csv"),
              row.names = FALSE)
    setwd(qc_base)
  }

  cat("\n")
}
