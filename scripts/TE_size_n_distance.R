# suppressMessages({
#     library(dplyr)
#     library(tidyr)
#     library(ggplot2)
#     library(DMRcaller)
#     library(GenomicFeatures)
#     library(plyranges)
#     library(parallel)
#     library(data.table)
# })

##################################################################################

TE_delta_meth <- function(pooled_list, TE_gr, context = NULL) {
    # use pooled and trimmed object
    meth_ctrl <- pooled_list[[1]]
    meth_trnt <- pooled_list[[2]]

    # delta
    meth_delta <- meth_ctrl[, 1]
    meth_delta$Proportion <- (meth_trnt$readsM / meth_trnt$readsN) - (meth_ctrl$readsM / meth_ctrl$readsN)
    meth_delta$Proportion[is.nan(meth_delta$Proportion)] <- 0

    if (is.null(context) | length(context) != 1) {
        context <- c("CG", "CHG", "CHH")
    }

    # plots
    te_meth_delta <- list()
    for (cntx in context) {
        # delta
        te_meth_delta[[cntx]] <- calculate_te_methylation(meth_delta, TE_gr, cntx, T)
        te_meth_delta[[cntx]]$sample <- "delta"
    }

    return(te_meth_delta)
}

##################################################################################

calculate_te_methylation <- function(meth_data, TE_gr, context, is.delta = F) {
    meth_data <- meth_data[meth_data$context == context]
    if (!is.delta) {
        meth_data$Proportion <- meth_data$readsM / meth_data$readsN
    }
    # Find overlaps between methylation data and TEs
    overlaps <- findOverlaps(TE_gr, meth_data)
    overlaps_df <- data.frame(
        te_id = queryHits(overlaps),
        meth_idx = subjectHits(overlaps)
    )

    # Get TE IDs and methylation levels for overlapping regions
    overlaps_df$te_id <- TE_gr$gene_id[overlaps_df$te_id]
    overlaps_df$meth_level <- meth_data$Proportion[overlaps_df$meth_idx]

    # Calculate average methylation level for each TE
    te_avg_meth <- overlaps_df %>%
        group_by(te_id) %>%
        summarise(avg_meth = mean(meth_level, na.rm = TRUE), .groups = "drop")
    te_avg_meth$context <- context

    te_width <- data.frame(
        te_id = TE_gr$gene_id,
        width = width(TE_gr),
        superfamily = TE_gr$Transposon_Super_Family
    )

    # te_width <- te_width %>% filter(width < 20000)

    as.data.frame(te_avg_meth) %>% merge(., te_width, by = "te_id")
}

##################################################################################

te_size_plot <- function(x_list, cntx, line_col = "red4", point_col = "gray20") {
    x_max <- max(x_list[[cntx]]$width, na.rm = TRUE)
    ggplot(x_list[[cntx]], aes(x = width, y = avg_meth, color = sample)) +
        geom_vline(xintercept = 4000, linetype = "dashed", color = "gray60") +
        geom_point(alpha = 0.6, size = 0.3, shape = 20) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        geom_smooth(
            method = lm, formula = y ~ splines::bs(x, 3),
            se = TRUE,
            linewidth = 0.5,
            color = line_col
        ) +
        scale_color_manual(values = c("delta" = point_col)) +
        labs(
            title = paste(cntx, "context"),
            x = "TE size (Kbp)",
            y = paste0("Δ ", "Methylation")
        ) +
        theme_classic() +
        theme(
            text = element_text(family = "serif"),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
        ) +
        scale_x_continuous(
            breaks = pretty(c(0, x_max), n = 5),
            limits = c(0, x_max), expand = c(0, 0),
            labels = function(x) round(x / 1000, 1)
        ) +
        scale_y_continuous(limits = c(-0.3, 0.3), expand = c(0, 0), breaks = c(-0.296, 0, 0.296), labels = c(-0.3, 0, 0.3))
}

##################################################################################

distance_from_centromer <- function(TE_meth_delta_list, TE_gr = TE_gr,
                                    centromeres = GRanges(), window_size = 1e6,
                                    lines_col = c("#3d53b4", "#3b8f3e", "#bb4949"),
                                    y_max = NULL, y_min = NULL, y_breaks = NULL) {
    if (length(centromeres) == 0L) {
        return(list(
            df = data.frame(),
            plot = ggplot() +
                annotate('text', x = 0, y = 0, label = 'Centromere coordinates not supplied') +
                theme_void()
        ))
    }
    cen_df <- as.data.frame(centromeres)[, c('seqnames', 'start', 'end')]
    cen_df$centromere <- (cen_df$start + cen_df$end) / 2
    cen_pos <- stats::setNames(cen_df$centromere, cen_df$seqnames)
    te_distance <- data.frame(
        chr = as.character(seqnames(TE_gr)),
        pos = (as.numeric(start(TE_gr)) + as.numeric(end(TE_gr))) / 2,
        te_id = TE_gr$gene_id,
        centromere = NA
    )
    te_distance$centromere <- unname(cen_pos[te_distance$chr])
    te_distance <- te_distance[!is.na(te_distance$centromere), , drop = FALSE]
    te_distance$distance <- abs(te_distance$centromere - te_distance$pos)
    all_cx_dis <- rbind(
        TE_meth_delta_list[["CG"]],
        TE_meth_delta_list[["CHG"]],
        TE_meth_delta_list[["CHH"]]
    )
    keep_col <- grep("te_id|distance", names(te_distance))
    te_distance_merged_out <- merge(te_distance[, keep_col], all_cx_dis, by = "te_id")
    te_distance_merged <- te_distance_merged_out
    te_distance_merged$distance <- te_distance_merged$distance / window_size

    ########

    te_distance_cntx <- rbind(
        te_distance_merged %>%
            filter(context == "CG") %>%
            mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
            group_by(window, context) %>%
            summarise(
                avg_meth = mean(avg_meth, na.rm = TRUE),
                distance = mean(distance, na.rm = TRUE),
                .groups = "drop"
            ),
        te_distance_merged %>%
            filter(context == "CHG") %>%
            mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
            group_by(window, context) %>%
            summarise(
                avg_meth = mean(avg_meth, na.rm = TRUE),
                distance = mean(distance, na.rm = TRUE),
                .groups = "drop"
            ),
        te_distance_merged %>%
            filter(context == "CHH") %>%
            mutate(window = floor(distance)) %>% # floor(distance / 1e6) * 1e6) %>%
            group_by(window, context) %>%
            summarise(
                avg_meth = mean(avg_meth, na.rm = TRUE),
                distance = mean(distance, na.rm = TRUE),
                .groups = "drop"
            )
    ) %>%
        as.data.frame() %>%
        filter(avg_meth != max(avg_meth, na.rm = TRUE))

    ########

    y_axis_limits <- if (!is.null(y_max) | !is.null(y_min)) {
        if (!is.null(y_breaks)) {
            y_breaks <- c(y_min, 0, y_max)
        }
        scale_y_continuous(
            limits = c(y_min, y_max),
            breaks = y_breaks
        )
    } else {
        geom_blank()
    }

    ########

    max_distance <- max(te_distance_cntx$distance, na.rm = TRUE)
    te_distance_plot <- ggplot(data = te_distance_cntx, aes(x = distance, y = avg_meth, color = context, group = context)) +
        geom_line(linewidth = 0.85) +
        theme_bw() + # theme_classic() +
        labs(
            title = " ",
            x = "Distance from centromer (Mbp)",
            y = "Δ Methylation"
        ) +
        theme(
            text = element_text(family = "serif"),
            legend.position = "none",
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(hjust = 0.5, size = 10),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
        ) +
        scale_x_continuous(
            limits = c(0, max_distance),
            breaks = pretty(c(0, max_distance), n = 4),
            expand = c(0, 0)
        ) +
        y_axis_limits +
        annotate("text",
            x = max_distance * 0.85,
            y = max(te_distance_cntx$avg_meth) * 0.98,
            label = c("CG", "\nCHG", "\n\nCHH"),
            hjust = 0, vjust = 0.75, size = 3.25,
            color = lines_col, fontface = "bold", family = "serif"
        )

    return(list(
        df = select(te_distance_merged_out, -sample),
        plot = te_distance_plot
    ))
}
