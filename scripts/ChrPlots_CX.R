chromosome_plot <- function(profile_vars, chr.n, y_max, y_mid, y_min, col, cntx) {
    ### filter by chromosome
    profile_vars <- lapply(profile_vars, function(x) {
        x[x@seqnames == chr.n]
    })

    ### if the seqlevels are not the same, align them
    ### by taking the max length of each seqlevel
    maxLens <- Reduce(
        function(a, b) pmax(a, b, na.rm = TRUE),
        lapply(profile_vars, seqlengths)
    )
    common_si <- Seqinfo(names(maxLens), maxLens)
    profile_vars <- lapply(profile_vars, function(gr) {
        seqlevels(gr) <- seqlevels(common_si) # ensure same ordering
        seqinfo(gr) <- common_si # overwrite lengths
        gr
    })

    ### make it as gr list
    methylationProfiles <- GRangesList()
    for (i in seq_along(profile_vars)) {
        methylationProfiles[[i]] <- profile_vars[[i]]
    }

    lwd_f <- ifelse(length(profile_vars) > 2, 1, 1.25)

    pos <- (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]])) / 2
    plot(pos, methylationProfiles[[1]]$mean_value,
        type = "o",
        ylim = c(y_min, y_max), xlab = "", ylab = "",
        col = col[1], pch = 26, lty = 1,
        yaxt = "n", xaxt = "n", lwd = lwd_f,
        main = "", frame.plot = FALSE
    )

    if (length(profile_vars) > 1) {
        for (i_p in 1:length(profile_vars)) {
            lines(pos, methylationProfiles[[i_p]]$mean_value,
                type = "o", col = col[i_p], lty = 1, lwd = lwd_f, pch = 26
            )
        }
    }

    if (cntx != "TE") {
        axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray20", tck = 0.025, pos = y_min) ###
    } else {
        # color bellow the lines
        polygon(c(pos, rev(pos)),
            c(methylationProfiles[[1]]$mean_value, rep(y_min, length(pos))),
            col = adjustcolor(col[1], alpha.f = 0.25), border = NA
        )
        axis(1, at = c(min(pos), max(pos)), labels = FALSE, col = "gray10", tck = 0.1, pos = 0) ###
    }

    if (cntx != "TE" & y_min != 0) {
        lines(pos, rep(0, length(pos)), type = "o", col = "gray50", lty = 2, pch = 26)
    }
}

###################################################################

ChrPlots_CX_all <- function(
    meth_var_list,
    y_max_cg = 1,
    y_max_chg = 0.5,
    y_max_chh = 0.2,
    y_mid_cg = NULL,
    y_mid_chg = NULL,
    y_mid_chh = NULL,
    y_min_cg = 0,
    y_min_chg = 0,
    y_min_chh = 0,
    legend_names = NULL,
    italic_legend_names = TRUE,
    ylab_suffix = NULL,
    y_title_cex = 1,
    chr_amount,
    chr_length,
    is_subCX,
    chr_labels = as.character(seq_len(chr_amount)),
    TE_as_gr = NULL) {
    ### color palette
    if (!is_subCX) {
        if (length(meth_var_list) == 2) {
            col_vec <- c("gray40", "#bf6828")
        } else {
            col_vec <- "gray20"
        }
    } else {
        col_vec <- brewer.pal(n = 9, name = "Set1")
    }

    ### y-mid edit
    y_mid_cg <- ifelse(is.null(y_mid_cg), (y_max_cg + y_min_cg) / 2, y_mid_cg)
    y_mid_chg <- ifelse(is.null(y_mid_chg), (y_max_chg + y_min_chg) / 2, y_mid_chg)
    y_mid_chh <- ifelse(is.null(y_mid_chh), (y_max_chh + y_min_chh) / 2, y_mid_chh)

    ### ylab suffix - in addition to 'CNTX methylation'
    ### add 'ylab_suffix=(delta)' to get 'CG methylation (delta)'
    if (!is.null(ylab_suffix)) {
        ylab_suffix <- paste0(" ", ylab_suffix)
    }

    ### Low resolution profiles plot ###
    ## column 1: y-axis strip
    ## columns 2-6: chromosomes 1-5
    ## column 7: legend

    if (is.null(TE_as_gr)) {
        panel_count <- 4 * (chr_amount + 1)
        lay <- cbind(
            matrix(seq_len(panel_count), nrow = 4, byrow = TRUE),
            rep(panel_count + 1, 4)
        )

        layout(lay,
            widths  = c(0.375, chr_length, 0.5),
            heights = c(1, 1, 1, 0.25) # last row for Chrs
        )
    } else {
        panel_count <- 5 * (chr_amount + 1)
        lay <- cbind(
            matrix(seq_len(panel_count), nrow = 5, byrow = TRUE),
            rep(panel_count + 1, 5)
        )

        layout(lay,
            widths  = c(0.375, chr_length, 0.5),
            heights = c(1, 1, 1, 0.35, 0.30)
        )
    }


    for (cntx in c("CG", "CHG", "CHH", "TE")) {
        meth_vars_context <- lapply(meth_var_list, function(inner_list) inner_list[[tolower(cntx)]])

        if (is_subCX) {
            meth_vars_context <- meth_vars_context[[1]]
            legend_col <- "gray20"
        } else {
            legend_col <- col_vec
        }

        if (cntx != "TE") {
            ## y-axis
            y_title <- paste0(cntx, " methylation", ylab_suffix)
            y_max_cntx <- ifelse(cntx == "CG", y_max_cg, ifelse(cntx == "CHG", y_max_chg, y_max_chh))
            y_mid_cntx <- ifelse(cntx == "CG", y_mid_cg, ifelse(cntx == "CHG", y_mid_chg, y_mid_chh))
            y_min_cntx <- ifelse(cntx == "CG", y_min_cg, ifelse(cntx == "CHG", y_min_chg, y_min_chh))

            par(mar = c(0, 4, 0, 0)) # c(1, 4, 2, 0))
            plot(NA, NA,
                xlim = c(0, 1), ylim = c(y_min_cntx, y_max_cntx),
                axes = FALSE, xlab = "",
                ylab = y_title,
                cex.lab = y_title_cex
            )
            axis(2,
                at = c(y_min_cntx, y_mid_cntx, y_max_cntx),
                labels = FALSE # ,  c(paste0("       ", y_min_cntx), y_mid_cntx, paste0(y_max_cntx, "      ")),
                # col = "gray35", cex.axis = 1
            )
            mtext(
                side = 2,
                text = y_min_cntx,
                at = y_min_cntx,
                adj = ifelse(y_min_cntx == 0, 0.5, 0), #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )
            mtext(
                side = 2,
                text = y_mid_cntx,
                at = y_mid_cntx,
                adj = 0.5, #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )
            mtext(
                side = 2,
                text = y_max_cntx,
                at = y_max_cntx,
                adj = ifelse(y_max_cntx == 0 | y_max_cntx == 1, 0.5, 0.9), #
                line = 0.65,
                col = "gray25",
                cex = 0.55
            )

            ## chromosome
            par(mar = c(0, 0, 0, 0)) # c(1, 0, 2, 0))
            for (chr in seq(chr_amount)) {
                suppressWarnings(chromosome_plot(meth_vars_context, chr, y_max_cntx, y_mid_cntx, y_min_cntx, col_vec, cntx))

                if (chr == 1 & cntx == "CG") {
                    ## legend
                    text(
                        x = 3, y = y_max_cntx,
                        labels = paste0("\n\n",legend_names[1]),
                        pos = 4, # right of point
                        col = legend_col[1],
                        font = 2, # bold
                        cex = 1.25
                    )
                    if (length(legend_names) > 1) {
                        text(
                            x = 3, y = y_max_cntx,
                            labels = paste0("\n\n\n\n",legend_names[2]),
                            pos = 4,
                            col = legend_col[2],
                            font = 2,
                            cex = 1.25
                        )
                    }
                }
            }
        } else if (!is.null(TE_as_gr)) {
            te_plot_conf(TE_as_gr, chr_amount)
        }
    }

    ## x-axis - chromosome
    par(mar = c(1, 0, 0, 0))
    plot.new()
    for (chr in seq(chr_amount)) {
        plot.new()
        mtext(
            side = 1,
            text = chr_labels[chr],
            line = -0.5,
            at = 0.5,
            adj = 0.5,
            col = "gray25"
        )
    }

    if (is_subCX) {
        ## legend
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend(x = 0, y = 1, # 1 - 0.0125,
            legend = toupper(names(meth_var_list[[1]][["cg"]])),
            # title = "CG",
            col = brewer.pal(n = 4, name = "Set1"),
            lty = 1, bty = "n", cex = 1, lwd = 2
        )
        legend(x = 0, y = 0.75, # 1 - 0.25 - (0.0125 * 2),
            legend = toupper(names(meth_var_list[[1]][["chg"]])),
            # title = "CHG",
            col = brewer.pal(n = 3, name = "Set1"),
            lty = 1, bty = "n", cex = 1, lwd = 2
        )
        legend(x = 0, y = 0.468, # 1 - 0.5 - (0.0125 * 5),
            legend = toupper(names(meth_var_list[[1]][["chh"]])),
            # title = "CHH",
            col = brewer.pal(n = 9, name = "Set1"),
            lty = 1, bty = "n", cex = 1, lwd = 2
        )
    }
}

###################################################################

te_plot_conf <- function(x, chr_amount) {
    te_gr <- x %>%
        circlize::genomicDensity(window.size = 0.5e6) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
        identity()
    names(mcols(te_gr)) <- "mean_value"

    y_max_te <- 1 # max(te_gr$mean_value)
    y_mid_te <- 0.5 # y_max_te / 2
    y_min_te <- 0

    par(mar = c(0, 4, 0, 0))
    plot(NA, NA,
        xlim = c(0, 1), ylim = c(y_min_te, y_max_te),
        axes = FALSE, xlab = "",
        ylab = "TEs", cex.lab = 1.2
    )

    axis(2, at = c(y_min_te, y_max_te), labels = FALSE)
    mtext(
        side = 2, text = c("", ""), # c(y_min_te, y_mid_te, round(y_max_te, 2)),
        at = c(y_min_te, y_max_te) # ,
        # line = 0.65, col = "gray25", cex = 0.7, adj = c(0, .5, .9)
    )

    par(mar = c(0, 0, 0, 0))
    te_list <- list(te_gr)
    te_names <- "TE"

    for (chr in seq(chr_amount)) {
        suppressWarnings(chromosome_plot(te_list, chr, y_max_te, y_mid_te, y_min_te, col = "#55555590", cntx = "TE"))
    }
}

###################################################################

windowSize_mcol <- function(x, mcol_name, windowSize = 2.5e5) {
    # chromosome lengths
    seqlens <- vapply(split(end(x), seqnames(x)), max, numeric(1))
    seqlengths(x) <- seqlens[seqlevels(x)]

    # windows by window size
    windows <- tileGenome(seqlens,
        tilewidth = windowSize,
        cut.last.tile.in.chrom = TRUE
    )

    # map to windows and calculate mean value
    hits <- findOverlaps(windows, x, ignore.strand = TRUE)
    mValue <- tapply(mcols(x)[[mcol_name]][subjectHits(hits)],
        queryHits(hits),
        mean,
        na.rm = TRUE
    )
    mcols(windows)$mean_value <- NA_real_
    mcols(windows)$mean_value[as.integer(names(mValue))] <- mValue

    windows
}

###################################################################

var_sep <- function(a, subCX = F, num_cores) {
    subContext_list <- list(
        CG = c("CGA", "CGT", "CGC", "CGG"),
        CHG = c("CAG", "CTG", "CCG"),
        CHH = c("CAA", "CAT", "CAC", "CTA", "CTT", "CTC", "CCA", "CCT", "CCC")
    )

    if (num_cores >= 27) {
        num_cores_1 <- 9
        num_cores_2 <- 3
    } else if (num_cores >= 18) {
        num_cores_1 <- 9
        num_cores_2 <- 2
    } else if (num_cores >= 9) {
        num_cores_1 <- 9
        num_cores_2 <- 1
    } else if (num_cores >= 3) {
        num_cores_1 <- 3
        num_cores_2 <- 1
    } else {
        num_cores_1 <- 1
        num_cores_2 <- 1
    }


    if (!subCX) {
        var_subSep <- function(gr, cntx) {
            gr[gr$context == cntx, ] %>%
                makeGRangesFromDataFrame(keep.extra.columns = T) %>%
                windowSize_mcol("Proportion")
        }
        mclapply(c("CG", "CHG", "CHH"), function(cntx) {
            var_subSep(a, cntx)
        }, mc.cores = num_cores) %>%
            setNames(c("cg", "chg", "chh"))
    } else {
        var_subSep <- function(gr, cntx) {
            gr[gr$trinucleotide_context == cntx, ] %>%
                makeGRangesFromDataFrame(keep.extra.columns = T) %>%
                windowSize_mcol("Proportion")
        }
        mclapply(names(subContext_list), function(cx) {
            setNames(
                mclapply(subContext_list[[cx]], function(subcx) var_subSep(a, subcx),
                    mc.cores = num_cores_1
                ),
                tolower(subContext_list[[cx]])
            )
        }, mc.cores = num_cores_2) %>%
            setNames(tolower(names(subContext_list)))
    }
}

###################################################################

run_ChrPlots_CX <- function(ctrl_name, trnt_name, ctrl_pool, trnt_pool, TE.gr, num_cores, CX_plot = T, subCX_plot = T) {
    cat("\n")
    cat(paste0("\rCalculate methylated/unmethylated C's ratio... [", ctrl_name, "]          "))
    ctrl_pool <- as.data.frame(ctrl_pool) %>%
        mutate(Proportion = readsM / readsN) %>%
        select(seqnames, start, end, Proportion, context, trinucleotide_context) %>%
        mutate(Proportion = ifelse(is.nan(Proportion), 0, Proportion)) # %>%
    # filter(!is.nan(Proportion)) ###

    cat(paste0("\rCalculate methylated/unmethylated C's ratio... [", trnt_name, "]          "))
    trnt_pool <- as.data.frame(trnt_pool) %>%
        mutate(Proportion = readsM / readsN) %>%
        select(seqnames, start, end, Proportion, context, trinucleotide_context) %>%
        mutate(Proportion = ifelse(is.nan(Proportion), 0, Proportion)) # %>%
    # filter(!is.nan(Proportion)) ###

    # Map arbitrary sequence names to stable internal panel indices while
    # keeping the original names for plot labels.
    chr_labels <- unique(c(as.character(ctrl_pool$seqnames), as.character(trnt_pool$seqnames)))
    chr_map <- stats::setNames(as.character(seq_along(chr_labels)), chr_labels)
    internal_levels <- as.character(seq_along(chr_labels))
    ctrl_pool$seqnames <- factor(unname(chr_map[as.character(ctrl_pool$seqnames)]), levels = internal_levels)
    trnt_pool$seqnames <- factor(unname(chr_map[as.character(trnt_pool$seqnames)]), levels = internal_levels)
    if (length(TE.gr)) {
        te_current <- as.character(seqlevels(TE.gr))
        te_present <- intersect(te_current, names(chr_map))
        TE.gr <- keepSeqlevels(TE.gr, te_present, pruning.mode = 'coarse')
        TE.gr <- renameSeqlevels(TE.gr, stats::setNames(unname(chr_map[te_present]), te_present))
    }

    # delta df
    cat(paste0("\rCalculate methylated/unmethylated C's ratio... [delta]          "))
    delta_pool <- ctrl_pool %>%
        mutate(Proportion = trnt_pool$Proportion - ctrl_pool$Proportion) %>%
        select(seqnames, start, end, Proportion, context, trinucleotide_context) %>%
        filter(!is.nan(Proportion))
    cat(paste0("\rCalculate methylated/unmethylated C's ratio:   done!            "))
    cat("\n")
    
    # normelize Chr panel size to its length
    cat(paste0("\rNormelize chromosome panel size to its length...    "))
    chr_length <- rbind(ctrl_pool, trnt_pool) %>%
        group_by(seqnames) %>%
        summarise(max_start = max(start), .groups = "drop") %>%
        pull(max_start)
    max_chr_length <- chr_length / max(chr_length)
    chr_amount <- length(chr_labels)
    cat(paste0("\rNormelize chromosome panel size to its length: done!"))
    cat("\n-------------\n")

    ########################################################################
    ### ChrPlot
    if (CX_plot) {
        cat(paste0("\rChrPlots... [", trnt_name, " & ", ctrl_name, "]"))
        img_device(paste0("ChrPlot_", trnt_name, "_vs_", ctrl_name), w = 7, h = 4)
        try({
            ChrPlots_CX_all(
                meth_var_list = list(var_sep(ctrl_pool, F, num_cores), var_sep(trnt_pool, F, num_cores)),
                y_max_cg = 1,
                y_max_chg = 0.6,
                y_max_chh = 0.25,
                y_mid_cg = NULL,
                y_mid_chg = NULL,
                y_mid_chh = NULL,
                y_min_cg = 0,
                y_min_chg = 0,
                y_min_chh = 0,
                legend_names = c(ctrl_name, trnt_name),
                italic_legend_names = FALSE,
                ylab_suffix = NULL,
                y_title_cex = 1,
                chr_amount = chr_amount,
                chr_labels = chr_labels,
                chr_length = max_chr_length,
                is_subCX = FALSE,
                TE_as_gr = TE.gr
            )
        })
        dev.off()

        # delta
        cat(paste0("\rChrPlots... [delta]            "))
        img_device(paste0("ChrPlot_difference_", trnt_name, "_vs_", ctrl_name), w = 7, h = 4)
        try({
            ChrPlots_CX_all(
                meth_var_list = list(var_sep(delta_pool, F, num_cores)),
                y_max_cg = 0.1,
                y_max_chg = 0.1,
                y_max_chh = 0.1,
                y_mid_cg = 0,
                y_mid_chg = 0,
                y_mid_chh = 0,
                y_min_cg = -0.1,
                y_min_chg = -0.1,
                y_min_chh = -0.1,
                legend_names = paste0(trnt_name, " vs ", ctrl_name),
                italic_legend_names = FALSE,
                ylab_suffix = "(Δ)",
                y_title_cex = 1,
                chr_amount = chr_amount,
                chr_labels = chr_labels,
                chr_length = max_chr_length,
                is_subCX = FALSE,
                TE_as_gr = TE.gr
            )
        })
        dev.off()
        cat(paste0("\rChrPlots:   done!              "))
        cat("\n")
    }
    ########################################################################
    ### ChrPlot sub-CX
    if (subCX_plot) {
        # var1
        cat(paste0("\rChrPlots for sub-contexts... [", ctrl_name, "]          "))
        img_device(paste0("subCX/ChrPlot_subCX_", ctrl_name), w = 7, h = 4)
        try({
            ChrPlots_CX_all(
                meth_var_list = list(var_sep(ctrl_pool, T, num_cores), var_sep(trnt_pool, T, num_cores)),
                y_max_cg = 1,
                y_max_chg = 0.6,
                y_max_chh = 0.25,
                y_mid_cg = NULL,
                y_mid_chg = NULL,
                y_mid_chh = NULL,
                y_min_cg = 0,
                y_min_chg = 0,
                y_min_chh = 0,
                legend_names = ctrl_name,
                italic_legend_names = FALSE,
                ylab_suffix = NULL,
                y_title_cex = 1,
                chr_amount = chr_amount,
                chr_labels = chr_labels,
                chr_length = max_chr_length,
                is_subCX = TRUE,
                TE_as_gr = TE.gr
            )
        })
        dev.off()

        # var2
        cat(paste0("\rChrPlots for sub-contexts... [", trnt_name, "]          "))
        img_device(paste0("subCX/ChrPlot_subCX_", trnt_name), w = 7, h = 4)
        try({
            ChrPlots_CX_all(
                meth_var_list = list(var_sep(trnt_pool, T, num_cores)),
                y_max_cg = 1,
                y_max_chg = 0.6,
                y_max_chh = 0.25,
                y_mid_cg = NULL,
                y_mid_chg = NULL,
                y_mid_chh = NULL,
                y_min_cg = 0,
                y_min_chg = 0,
                y_min_chh = 0,
                legend_names = trnt_name,
                italic_legend_names = FALSE,
                ylab_suffix = NULL,
                y_title_cex = 1,
                chr_amount = chr_amount,
                chr_labels = chr_labels,
                chr_length = max_chr_length,
                is_subCX = TRUE,
                TE_as_gr = TE.gr
            )
        })
        dev.off()

        # delta
        cat("\rChrPlots for sub-contexts... [delta]          ")
        img_device(paste0("subCX/ChrPlot_difference_subCX_", trnt_name, "_vs_", ctrl_name), w = 7, h = 4)
        try({
            ChrPlots_CX_all(
                meth_var_list = list(var_sep(delta_pool, T, num_cores)),
                y_max_cg = 0.1,
                y_max_chg = 0.1,
                y_max_chh = 0.1,
                y_mid_cg = 0,
                y_mid_chg = 0,
                y_mid_chh = 0,
                y_min_cg = -0.1,
                y_min_chg = -0.1,
                y_min_chh = -0.1,
                legend_names = paste0(trnt_name, " vs ", ctrl_name),
                italic_legend_names = FALSE,
                ylab_suffix = "(Δ)",
                y_title_cex = 1,
                chr_amount = chr_amount,
                chr_labels = chr_labels,
                chr_length = max_chr_length,
                is_subCX = TRUE,
                TE_as_gr = TE.gr
            )
        })
        dev.off()
        cat("\rChrPlots for sub-contexts:   done!            ")
    }
    cat("\n")
}
