run_sum_deltaH_CX <- function(ctrl_name, trnt_name, ctrl_pool, trnt_pool, annotation.gr, TE.gr, description_df, num_cores, fdr = 0.95) {
    # genome_ann_path <- getwd()

    ctrl_pool <- as.data.frame(ctrl_pool) %>%
        select(seqnames, start, context, m_ctrl = readsM, n_ctrl = readsN)
    trnt_pool <- as.data.frame(trnt_pool) %>%
        select(seqnames, start, context, m_trnt = readsM, n_trnt = readsN)

    merged_pool <- inner_join(
        ctrl_pool, trnt_pool,
        by = c("seqnames", "start", "context")
    )

    seqlengths_merged <- merged_pool %>%
        group_by(seqnames) %>%
        summarise(seqlength = max(start)) %>%
        arrange(match(seqnames, unique(merged_pool$seqnames)))

    # main loop
    cat("\nsum dH in 100bp windowsize:")
    hd_list <- mclapply(
        c("CG", "CHG", "CHH"), function(cntx) {
            dH_bins(merged_pool, cntx)
        },
        mc.cores = ifelse(num_cores >= 3, 3, num_cores)
    )
    cg_hd <- hd_list[[1]]
    chg_hd <- hd_list[[2]]
    chh_hd <- hd_list[[3]]

    # filter by FDR threshold
    q_cg <- quantile(cg_hd$sum_surprisal, fdr) %>% as.numeric()
    q_chg <- quantile(chg_hd$sum_surprisal, fdr) %>% as.numeric()
    q_chh <- quantile(chh_hd$sum_surprisal, fdr) %>% as.numeric()
    cg_filtered <- cg_hd %>%
        filter(sum_surprisal > q_cg) %>%
        mutate(regionType = ifelse(pi_log2FC > 0, "gain", "loss"))
    chg_filtered <- chg_hd %>%
        filter(sum_surprisal > q_chg) %>%
        mutate(regionType = ifelse(pi_log2FC > 0, "gain", "loss"))
    chh_filtered <- chh_hd %>%
        filter(sum_surprisal > q_chh) %>%
        mutate(regionType = ifelse(pi_log2FC > 0, "gain", "loss"))

    # print df (both to terminal and .log file)
    df_2_print <- data.frame(
        context = c("CG", "CHG", "CHH"),
        total = c(nrow(cg_hd), nrow(chg_hd), nrow(chh_hd)),
        top_5 = c(nrow(cg_filtered), nrow(chg_filtered), nrow(chh_filtered)),
        threshold = paste0("S>", c(round(q_cg, 1), round(q_chg, 1), round(q_chh, 1)))
    )
    print(knitr::kable(df_2_print))
    message("\n", paste(
        c(
            paste(names(df_2_print), collapse = "\t\t"),
            paste(sapply(names(df_2_print), function(nm) paste(rep("=", nchar(nm)), collapse = "")), collapse = "\t\t"),
            apply(df_2_print, 1, function(row) paste(row, collapse = "\t\t"))
        ),
        collapse = "\n"
    ), "\n")

    # genome density plot
    dHRs_circular_plot(cg_filtered, chg_filtered, chh_filtered, annotation.gr, TE.gr, paste0(trnt_name, "_vs_", ctrl_name))

    # manhattan plots
    manH_par <- mclapply(
        list(list(cg_hd, "CG", fdr), list(chg_hd, "CHG", fdr), list(chh_hd, "CHH", fdr)),
        function(x) {
            mannh_plots(x[[1]], x[[2]], x[[3]])
        },
        mc.cores = ifelse(num_cores >= 3, 3, num_cores)
    )
    png(paste0("sum_dH_manhattan_plot_", trnt_name, "_vs_", ctrl_name, ".png"), width = 10, height = 2.75, units = "in", res = 300)
    multiplot(
        print(manH_par[[1]]), print(manH_par[[2]]), print(manH_par[[3]]),
        cols = 3
    )
    dev.off()

    # write bigWig file function
    write_bw <- function(x, cntx) {
        x <- x %>%
            select(seqnames, start, end, sum_surprisal) %>%
            filter(!is.na(sum_surprisal)) %>%
            mutate(score = as.numeric(sum_surprisal)) %>%
            makeGRangesFromDataFrame(., keep.extra.columns = TRUE, )
        seqlengths(x) <- setNames(seqlengths_merged$seqlength, seqlengths_merged$seqnames)
        rtracklayer::export.bw(x, con = paste0("surprisal_", cntx, "_", trnt_name, "_vs_", ctrl_name, ".bw"))
    }

    # annotate regions and write files
    ann_list <- genome_ann(annotation.gr, TE.gr)
    write_sum_dH <- function(filtered_df, cntx) {
        filtered_df <- makeGRangesFromDataFrame(filtered_df, keep.extra.columns = T) %>% as.data.frame() # save as gr object configuration
        write.csv(filtered_df, paste0("SurpMRs_", cntx, "_", trnt_name, "_vs_", ctrl_name, ".csv"), row.names = F)
        write_bw(filtered_df, cntx)
        # setwd(genome_ann_path)
        setwd("genome_annotation")
        SurpMRs_bins <- makeGRangesFromDataFrame(filtered_df, keep.extra.columns = T)
        suppressMessages(DMRs_ann(ann_list, SurpMRs_bins, cntx, description_df, sum_dH = T))
        DMRs_ann_plots(ctrl_name, trnt_name, cntx, sum_dH = T)
        TE_ann_plots(cntx, TE.gr, sum_dH = T)
        TE_Super_Family_Frequency(cntx, TE.gr, sum_dH = T)
        setwd("../")
    }

    write_sum_dH(cg_filtered, "CG")
    write_sum_dH(chg_filtered, "CHG")
    write_sum_dH(chh_filtered, "CHH")
}

########################################################################

dH_bins <- function(joint, cntx, min_coverage = 6) {
    thresh <- c(CG = 0.4, CHG = 0.2, CHH = 0.1)[cntx]
    min_C_sites <- c(CG = 4, CHG = 4, CHH = 4)[cntx]

    cat(paste0("\n[", cntx, "]\tmin coverage: ", min_coverage, "; min proportion value: ", thresh, "; min significant sites: ", min_C_sites))
    ### per-site surprisal value
    joint <- joint %>%
        filter(
            context == cntx,
            n_ctrl >= min_coverage,
            n_trnt >= min_coverage
        ) %>%
        # mutate(
        #     delta = abs(m_trnt / n_trnt - m_ctrl / n_ctrl),
        # ) %>%
        # filter(delta > thresh) %>%
        mutate(
            fc_ctrl = m_ctrl / n_ctrl,
            fc_trnt = m_trnt / n_trnt
        ) %>%
        filter(
            fc_ctrl > thresh,
            fc_trnt > thresh
        ) %>%
        calculate_surp_1()
        # calculate_surp_2()

    ### 100bp windowSize
    gr_sites <- GRanges(joint$seqnames, IRanges(joint$start, width = 1))
    seqlevels <- unique(joint$seqnames)
    seqlengths_per_chr <- joint %>%
        group_by(seqnames) %>%
        summarise(seqlength = max(start)) %>%
        arrange(match(seqnames, seqlevels))

    seqlengths(gr_sites) <- seqlengths_per_chr$seqlength
    names(seqlengths(gr_sites)) <- seqlevels
    tiles100 <- unlist(tileGenome(seqlengths(gr_sites), tilewidth = 100))

    ### sum surprisal within each window                  ###
    ov <- findOverlaps(gr_sites, tiles100, ignore.strand = TRUE)
    joint$tile_id <- subjectHits(ov) # map each site → window

    window_surprisal <- joint %>%
        group_by(tile_id) %>%
        summarise_surp_1()
        # summarise_surp_2(., min_C_sites)

    # final df
    out <- cbind(as.data.frame(tiles100)[window_surprisal$tile_id, 1:3], window_surprisal) %>%
        mutate(context = cntx) %>%
        dplyr::relocate(context, .after = tile_id) %>%
        select(-tile_id)
    out
}

########################################################################
# option 1 [S-Value = −log P_Bin(c | n, π0)]
calculate_surp_1 <- function(x) {
    mutate(x,
        pi0 = pmin(pmax(m_ctrl / n_ctrl, 1e-6), 1 - 1e-6), # ctrl expectation π0
        surprisal = -(lchoose(n_trnt, m_trnt) +
            m_trnt * log(pi0) +
            (n_trnt - m_trnt) * log(1 - pi0))
    )
}

# option 2 [H = p * log(1/p)]
calculate_surp_2 <- function(x) {
    mutate(x,
        surprisal_ctrl = fc_ctrl * log2(1 / (fc_ctrl + 1e-6)),
        surprisal_trnt = fc_trnt * log2(1 / (fc_trnt + 1e-6)),
        surprisal_site = surprisal_trnt + surprisal_ctrl,
        surprisal_diff = surprisal_trnt - surprisal_ctrl
    )
}

# option 3 [H = p * log(1/q)]
calculate_surp_3 <- function(x) {
    mutate(x,
        surprisal_cros = fc_trnt * log2(1 / (fc_ctrl + 1e-6)),
        surprisal_site = surprisal_trnt + surprisal_ctrl,
        surprisal_diff = surprisal_trnt - surprisal_ctrl
    )
}

# option 4 - KL-divergence [D = (p*log(1/q)) - (p*log(1/p))]
calculate_surp_4 <- function(x) {
    mutate(x,
        surprisal_cros = fc_trnt * log2(1 / (fc_ctrl + 1e-6)),
        surprisal_trnt = fc_trnt * log2(1 / (fc_trnt + 1e-6)),
        surprisal_kldv = surprisal_cros - surprisal_trnt,
        surprisal_site = surprisal_trnt + surprisal_ctrl,
        surprisal_diff = surprisal_trnt - surprisal_ctrl
    )
}


########################################################################

summarise_surp_1 <- function(x) {
    summarise(x,
        sum_surprisal = sum(surprisal, na.rm = TRUE),
        m_ctrl = sum(m_ctrl),
        n_ctrl = sum(n_ctrl),
        m_trnt = sum(m_trnt),
        n_trnt = sum(n_trnt),
        # adding 1 m_reads and 2 n_reads (total) if the count is 0
        pi_ctrl = (m_ctrl + (m_ctrl == 0)) / (n_ctrl + 2 * (m_ctrl == 0)),
        pi_trnt = (m_trnt + (m_trnt == 0)) / (n_trnt + 2 * (m_trnt == 0)),
        pi_log2FC = log2(pi_trnt / pi_ctrl),
        n_sites = n()
    ) %>%
        ungroup()
}

summarise_surp_2 <- function(x, min_c) {
    summarise(x,
        sum_surprisal_ctrl = sum(surprisal_ctrl),
        sum_surprisal_trnt = sum(surprisal_trnt),
        sum_surprisal_site = sum(surprisal_site),
        sum_surprisal_diff = sum(surprisal_diff),
        p_wilcoxon_diff = wilcox.test(surprisal_diff)$p.value,
        n_sites = n()
    ) %>%
        ungroup() %>%
        filter(n_sites >= min_c) %>%
        mutate(
            # FDR_diff = p.adjust(p_wilcoxon_diff, "BH"),
            call = case_when(
                p_wilcoxon_diff < 0.05 & sum_surprisal_diff > 0 ~ "hyperchaotic",
                p_wilcoxon_diff < 0.05 & sum_surprisal_diff < 0 ~ "hypochaotic",
                TRUE ~ "NS"
            )
        )
}


########################################################################

dHRs_circular_plot <- function(CG_df, CHG_df, CHH_df, ann.file, TE_4_dens, comparison_name) {
    chromosomes <- as.character(seqlevels(ann.file))
    chr_amount <- length(chromosomes)

    genes_type <- ann.file[which(ann.file$type == "gene")]

    ############# the plot #############
    img_device(paste0("dHR_Density_", comparison_name), w = 3.25, h = 3.25)

    circos.par(start.degree = 90)
    circos.genomicInitialize(as.data.frame(ann.file)[, 1:3], sector.names = chromosomes, axis.labels.cex = 0.325, labels.cex = 1.35)

    circos.genomicDensity(
        list(
            CG_df[CG_df$pi_log2FC > 0, 1:3],
            CG_df[CG_df$pi_log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            CHG_df[CHG_df$pi_log2FC > 0, 1:3],
            CHG_df[CHG_df$pi_log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            CHH_df[CHH_df$pi_log2FC > 0, 1:3],
            CHH_df[CHH_df$pi_log2FC < 0, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("#FF000080", "#304ed180"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.genomicDensity(
        list(
            as.data.frame(genes_type)[, 1:3],
            as.data.frame(TE_4_dens)[, 1:3]
        ),
        bg.col = "#fafcff", bg.border = NA, count_by = "number",
        col = c("gray80", "#fcba0320"), border = T, track.height = 0.165, track.margin = c(0, 0)
    )

    circos.clear()
    dev.off()
}

########################################################################

mannh_plots <- function(df, cntx, fdr = 0.95) {
    S_threshold <- as.numeric(quantile(df$sum_surprisal, fdr))

    # add cumulative genome coordinate
    seq_names <- as.character(df$seqnames)
    chr_order <- if (is.factor(df$seqnames)) {
        levels(droplevels(df$seqnames))
    } else {
        unique(seq_names)
    }
    chr_lengths <- vapply(chr_order, function(chr) max(df$end[seq_names == chr]), numeric(1))
    chr_offset <- c(0, cumsum(chr_lengths[-length(chr_lengths)]))
    names(chr_offset) <- chr_order

    df$coord <- df$start + chr_offset[seq_names]
    df$chr_index <- match(seq_names, chr_order)

    manH_p <- ggplot(df, aes(coord, sum_surprisal,
        colour = as.factor(chr_index %% 2)
    )) +
        geom_point(size = 0.025) +
        geom_hline(yintercept = S_threshold, color = "#bf6828", linetype = "dashed", size = 0.7) +
        scale_colour_manual(values = c("grey40", "#7ca182"), guide = "none") +
        scale_x_continuous(
            breaks = chr_offset + chr_lengths / 2,
            labels = chr_order
        ) +
        expand_limits(y = 0) +
        labs(
            x = "", # "Chromosome",
            y = "Window surprisal"
        ) +
        theme_bw()

    manH_p
}
