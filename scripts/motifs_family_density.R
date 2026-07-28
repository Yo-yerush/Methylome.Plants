TF_motifs <- function(jointed_gr, context = "all", windowSize = 1e6, DMP_fdr = 0.05, ann.gr, tfbs_path = "https://github.com/Yo-yerush/Methylome.At/raw/refs/heads/main/annotation_files/TAIR10_compressed_TFBSs.bed.gz") {
    if (context != "all") {
        jointed_gr <- jointed_gr[which(jointed_gr$context == context)]
        out_file_name <- paste0(context, "_TFs_superfamilies_DMPs_density")
    } else {
        context <- ""
        out_file_name <- "TFs_superfamilies_DMPs"
    }

    ##################### read TFBS file
    cat("read TFBS file\n")
    tfbs_data <- fread(tfbs_path, showProgress = F)

    colnames(tfbs_data) <- c("seqnames", "start", "end", "motif", "phase", "strand", "thickStart", "thickEnd", "itemRgb")

    # extract moif column information
    tfbs_data$name_parsed <- strsplit(as.character(tfbs_data$motif), "_")
    tfbs_data$Experiment <- sapply(tfbs_data$name_parsed, `[`, 1)
    tfbs_data$Condition <- sapply(tfbs_data$name_parsed, `[`, 2)
    tfbs_data$TF <- sapply(tfbs_data$name_parsed, `[`, 3)
    tfbs_data$MotifID <- sapply(tfbs_data$name_parsed, `[`, 4)

    ##################### group TF families
    tf_family <- c(
        # MADS-box
        "AG" = "MADS-box", "AP1" = "MADS-box", "AP3" = "MADS-box", "PI" = "MADS-box",
        "SEP3" = "MADS-box", "SOC1" = "MADS-box", "SVP" = "MADS-box", "AGL27" = "MADS-box",
        "FLC" = "MADS-box",

        # LFY
        "LFY" = "LFY/FLO",

        # bHLH
        "PIF1" = "bHLH", "PIF3" = "bHLH", "PIF4" = "bHLH", "PIF5" = "bHLH",

        # WRKY
        "WRKY18" = "WRKY", "WRKY33" = "WRKY", "WRKY40" = "WRKY",

        # bZIP
        "ABF1" = "bZIP", "ABF3" = "bZIP", "ABF4" = "bZIP",
        "GBF2" = "bZIP", "GBF3" = "bZIP", "HY5" = "bZIP",

        # Homeobox
        "ATHB-5" = "HD-ZIP", "ATHB-6" = "HD-ZIP",
        "ATHB-7" = "HD-ZIP", "HAT22" = "HD-ZIP",

        # GARP
        "ARR1" = "GARP", "ARR10" = "GARP", "ARR14" = "GARP", "KAN1" = "GARP",

        # MYB
        "MYB3" = "MYB", "CCA1" = "MYB",

        # Other TF families
        "SPL7" = "Others", # SBP/SPL
        "IDD4" = "Others", # IDD/BIRD
        "JKD" = "Others", # IDD/BIRD
        "TCP4" = "Others",
        "BPC1" = "Others",
        "E2FA" = "Others",
        "FHY3" = "Others",

        # Chromatin regulator (keep separate)
        "REF6" = "REF6"
    )

    tfbs_data$TF_family <- tf_family[tfbs_data$TF]

    tfbs_gr <- tfbs_data %>%
        mutate(
            seqnames = gsub("chr", "Chr", tfbs_data$seqnames),
            hex = sapply(strsplit(as.character(tfbs_data$itemRgb), ","), function(x) {
                rgb(as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]), maxColorValue = 255)
            })
        ) %>%
        filter(seqnames != "ChrMt") %>%
        filter(seqnames != "ChrPt") %>%
        dplyr::select(-c(motif, name_parsed, itemRgb, thickStart, thickEnd, phase)) %>%
        dplyr::relocate(strand, TF, TF_family, .after = end) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = T)


    # m1 <- findOverlaps(meth_var1, tfbs_gr)
    # var1_meth_tfbs <- meth_var1[queryHits(m1)]
    # mcols(var1_meth_tfbs) <- cbind.data.frame(
    #     mcols(var1_meth_tfbs),
    #     mcols(tfbs_gr[subjectHits(m1)])
    # )

    # m2 <- findOverlaps(meth_var2, tfbs_gr)
    # var2_meth_tfbs <- meth_var2[queryHits(m2)]
    # mcols(var2_meth_tfbs) <- cbind.data.frame(
    #     mcols(var2_meth_tfbs),
    #     mcols(tfbs_gr[subjectHits(m2)])
    # )

    cat("overlap joined methylation data with TFBS file\n")
    m <- findOverlaps(jointed_gr, tfbs_gr)
    joined_meth_tfbs <- jointed_gr[queryHits(m)]
    mcols(joined_meth_tfbs) <- cbind.data.frame(
        mcols(joined_meth_tfbs),
        mcols(tfbs_gr[subjectHits(m)])
    )


    ##################### call DMPs
    cat("call DMPs\n")
    dmp_gr <- calling_DMPs(joined_meth_tfbs, min_cov = 6, fdr = DMP_fdr)
    write.csv(dmp_gr, paste0(out_file_name, ".csv"), row.names = F)

    ##################### edit for the plot
    fam_size_order <- names(sort(table(dmp_gr$TF_family), decreasing = T))
    windowSize_legend_name <- gsub("e\\+0*", "E", format(windowSize, scientific = TRUE, upper.case = TRUE))

    chr_amount <- length(seqnames(ann.gr)@values)

    # heterochromatin positions (TAIR)
    heteroChr <- data.frame(
        Chr = paste0("Chr", c(1:chr_amount)),
        start = c(12500000, 1250000, 11000000, 1666667, 9444444),
        end = c(17500000, 7500000, 16250000, 7000000, 15000000)
    )

    # centromere positions (TAIR)
    cenChr <- data.frame(
        Chr = paste0("Chr", c(1:chr_amount)),
        start = c(14476796, 3462971, 13780083, 3177188, 11207348),
        end = c(15081019, 3650512, 14388500, 3248799, 11555278)
    )

    ##################### plot
    img_device(out_file_name, w = 4.25, h = 4.25)

    par(mar = c(0, 0, 0, 0))

    circos.par(gap.degree = c(rep(5, chr_amount - 1), 40), start.degree = 90, points.overflow.warning = FALSE)
    circos.genomicInitialize(as.data.frame(ann.gr)[, 1:3], sector.names = paste0("Chr ", 1:chr_amount), axis.labels.cex = 0.4, labels.cex = 1.25)
    for (family.i in fam_size_order) {
        cat(context, ">", family.i, "\n")
        message(time_msg(), context, ":\t", family.i)

        suppressMessages({
            density_data <- genomicDensity(dmp_gr[which(dmp_gr$TF_family == family.i)], window.size = windowSize, count_by = "number")
            ylims <- range(density_data$value) * 1.435

            # check if any chromosomes are missing and add it
            missing_chrs <- setdiff(as.character(unique(seqnames(jointed_gr))), unique(density_data$chr))
            if (length(missing_chrs) != 0)
            density_data <- rbind(
                density_data,
                data.frame(chr = missing_chrs, start = 1, end = 1, value = 0)
            )

            circos.genomicTrackPlotRegion(density_data, ylim = range(density_data$value), bg.border = NA, track.height = ifelse(length(fam_size_order) > 7, 0.075, 1), track.margin = c(0, 0), panel.fun = function(region, value, ...) {
                chr.n <- gsub("Chr", "", get.cell.meta.data("sector.index"))

                ### heterocromatin
                circos.rect(heteroChr[chr.n, 2], ylims[1], heteroChr[chr.n, 3], ylims[2], # xleft, ybottom, xright, ytop
                    col = "#fcba0320",
                    border = NA
                )
                ### centromere
                circos.rect(cenChr[chr.n, 2], ylims[1], cenChr[chr.n, 3], ylims[2], # xleft, ybottom, xright, ytop
                    col = "#fcba0360",
                    border = NA
                )

                ### density lines
                colors <- ifelse(value >= 15, "#440154",
                    ifelse(value >= 10, "#31688e",
                        ifelse(value >= 7, "#21918c",
                            ifelse(value >= 5, "#35b779",
                                ifelse(value >= 3, "#90d743",
                                    "#c2c2c2"
                                )
                            )
                        )
                    )
                )
                circos.genomicLines(region, value,
                    col = colors,
                    border = TRUE, lty = 1, lwd = 0.5, type = "h"
                )
            })
            ### y-axis labels
            circos.text("Chr1", x = 0, y = 0.5, labels = paste0(family.i, "  "), facing = "downward", cex = 0.45, adj = c(0.85, -0.15))
        })
    }
    circos.clear()

    legend("topleft",
        legend = c("≥15", "≥10", "≥7", "≥5", "≥3", "<3"),
        fill = c("#440154", "#31688e", "#21918c", "#35b779", "#90d743", "#d9d9d9"),
        title = paste0("TFs count"),
        cex = 0.6,
        bty = "n"
    )

    dev.off()
}
