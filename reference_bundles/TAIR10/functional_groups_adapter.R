################################
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
### TAIR10 legacy functional-group adapter (main function is below)
################
offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms

    child_terms_0 <- as.character(xx[[go_id]])
    child_terms <- child_terms_0

    for (i in 1:length(child_terms_0)) {
        child_terms <- c(child_terms, as.character(xx[[child_terms[i]]]))
    }

    return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
}

################
grep_position <- function(go_ids_vec, df) {
    vec <- NULL
    for (terms_l in go_ids_vec) {
        vec <- c(vec, grep(terms_l, df$GO.biological.process))
    }
    return(unique(vec))
}

################
clean_ASCII <- function(x) {
    x <- gsub("\001", " ", x)
    x <- gsub("\002", " ", x)
    x <- gsub("\036", " ", x)

    # x = gsub("[[:punct:]]", " ", x)
    # x = iconv(x, from = 'UTF-8', to = 'ASCII')
    return(x)
}

################
remove_dup_DMR <- function(y) {
    y <- as.character(unique(unlist(strsplit(y, ","))))
    paste(y, collapse = ",")
}

################################
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
### main functions

DMRs_into_groups <- function(
    treatment,
    ann,
    context = "all",
    DMRs_ann_dir = "./genome_annotation",
    datasets_dir = "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/Arabidopsis/Groups_genes_list",
    gene_sets_path = NULL) {

    if (!is.null(gene_sets_path)) {
        return(DMRs_into_generic_groups(
            treatment, ann, context,
            gene_sets_path = gene_sets_path,
            DMRs_ann_dir = DMRs_ann_dir
        ))
    }

    out_dir <- paste0(DMRs_ann_dir, "/functional_groups/", context)
    dir.create(out_dir, showWarnings = FALSE)

    ### load DMRs files
    if (context != "all") {
        DMR_file <- read.csv(paste0(DMRs_ann_dir, "/", context, "/DMRs_", ann, "_", context, "_genom_annotations.csv"))
        cntx_file <- context
    } else {
        DMR_file <- rbind(
            read.csv(paste0(DMRs_ann_dir, "/CG/DMRs_", ann, "_CG_genom_annotations.csv")),
            read.csv(paste0(DMRs_ann_dir, "/CHG/DMRs_", ann, "_CHG_genom_annotations.csv")),
            read.csv(paste0(DMRs_ann_dir, "/CHH/DMRs_", ann, "_CHH_genom_annotations.csv"))
        )
        cntx_file <- "all_contexts"
    }


    ################################################
    # # # # # # # # # # # # # # # #
    ### load genes-group datasets

    ################
    # from RA lab in github

    ################
    # RdDM pathway
    rddm <- rbind(
        read.table(paste0(datasets_dir, "/rddm_Matzke_et_al.txt"),
            sep = "\t", header = T
        ),
        read.table(paste0(datasets_dir, "/rddm_Cuerda_n_Slotkin.txt"),
            sep = "\t", header = T
        )
    ) %>% distinct()

    names(rddm) <- "gene_id"

    ################
    # Histone Lysine Methyltransferases
    HLM <- rbind(
        read.table(paste0(datasets_dir, "/histone_lysin_MTs_PONTVIANNE_et_al.txt"),
            sep = "\t", header = T
        ),
        read.table(paste0(datasets_dir, "/set_domain_danny_wd_et_al.txt"),
            sep = "\t", header = T
        )
    ) %>%
        distinct()

    names(HLM)[1] <- "gene_id"

    ################
    # Royal Family Proteins
    RF <- read.table(paste0(datasets_dir, "/At_Agenet_Tudor_family_brazil_et_al.txt"),
        sep = "\t", header = T
    )
    names(RF)[1] <- "gene_id"
    RF$gene_id <- gsub("DUF\\d+", "", RF$gene_id)

    ################
    # primary/secondary metabolism (https://doi.org/10.1093/gbe/evv217)
    primary_metabolism_0 <- read.csv(paste0(datasets_dir, "/primary_metabolism_mukherjee_et_al.csv"))

    secondary_metabolism_0 <- read.csv(paste0(datasets_dir, "/secondary_metabolism_mukherjee_et_al.csv"))

    ################
    # Seed dpecific genes
    # https://bar.utoronto.ca/ExpressionAngler/
    gene_list_seed_specific <- read.csv(paste0(datasets_dir, "/seed_specific_TAIR_Dry_seed_n_all_Stages.csv")) %>%
        # filter(r_value >= 0.6) %>%
        dplyr::rename(gene_id = tair_id)

    ################
    # GO term offsprings
    # child_terms_epigenetic <- offspring_fun("GO:0040029") # epigenetic regulation of_gene expression
    child_terms_chromatin_rem <- offspring_fun("GO:0006338") # chromatin remodeling
    child_terms_stress <- offspring_fun("GO:0006950") # response to stress
    child_terms_biotic <- offspring_fun("GO:0009607") # response to biotic stimulus
    child_terms_abiotic <- offspring_fun("GO:0009628") # response to abiotic stimulus

    ################

    ################################################
    # # # # # # # # # # # # # # # #
    ### define groups of genes from DMRs results

    ################
    # GO grep positions of chiled ters
    grep_child_stress <- DMR_file[grep_position(child_terms_stress, DMR_file), ]
    grep_child_biotic <- DMR_file[grep_position(child_terms_biotic, DMR_file), ]
    grep_child_abiotic <- DMR_file[grep_position(child_terms_abiotic, DMR_file), ]
    grep_child_chromatin_rem <- DMR_file[grep_position(child_terms_chromatin_rem, DMR_file), ]

    ################
    # RdDM pathway
    rddm_group <- merge.data.frame(rddm, DMR_file, by = "gene_id") %>%
        arrange(pValue)

    ################
    # Histone Lysine Methyltransferases
    HLM_group <- merge.data.frame(HLM, DMR_file, by = "gene_id") %>%
        arrange(pValue)

    ################
    # Royal Family Proteins
    RF_group <- merge.data.frame(RF, DMR_file, by = "gene_id") %>% dplyr::select(-Agenet.Tudor.Domain, -Other.Domain)

    ################
    # DNA de-Methylases
    DNA_deMTs <- DMR_file[grep("AT4G34060|AT5G04560|AT2G36490|AT3G10010", DMR_file$gene_id), ] %>%
        distinct(gene_id, .keep_all = T) %>%
        arrange(pValue)

    # histone de-Methylases
    histone_deMTs <- DMR_file[grep("LDL|FLD|ELF|IBM|JMJ|REF", DMR_file$Symbol), ] %>%
        distinct(gene_id, .keep_all = T) %>%
        arrange(pValue)

    ################
    # REM Transcription Factor Family (https://www.arabidopsis.org/browse/gene_family/REM)
    # added VDD () and VAL (AT5G60140)
    REM_TFs <- DMR_file[grep("AT4G31610|AT2G24700|AT2G24690|AT2G24680|AT2G24650|AT2G24630|AT4G33280|AT1G26680|AT2G46730|AT1G49480|AT4G31620|AT3G53310|AT3G06220|AT3G46770|AT5G09780|AT4G31630|AT4G31640|AT4G31650|AT4G31660|AT4G31680|AT4G31690|AT5G18000|AT5G60140", DMR_file$gene_id), ] %>%
        distinct(gene_id, .keep_all = T) %>%
        arrange(pValue)

    ################
    # primary/secondary metabolism (https://doi.org/10.1093/gbe/evv217)
    primary_metabolism_v <- distinct(primary_metabolism_0, gene_id) %>%
        merge.data.frame(., DMR_file, by = "gene_id") %>%
        arrange(pValue)

    secondary_metabolism_v <- distinct(secondary_metabolism_0, gene_id) %>%
        merge.data.frame(., DMR_file, by = "gene_id") %>%
        arrange(pValue)

    ################
    # Seed dpecific genes
    # https://bar.utoronto.ca/ExpressionAngler/
    seed_specific_genes <- merge.data.frame(gene_list_seed_specific, DMR_file, by = "gene_id") %>%
        filter(r_value >= 0.7) %>%
        dplyr::select(-r_value)

    ################################################
    # # # # # # # # # # # # # # # #
    ### create final list

    final_0_list <- list(
        DNA_methyltransferase = rbind(
            DMR_file[grep("^2\\.1\\.1\\.37", DMR_file$EC), ], # EC_2.1.1.37
            DMR_file[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(DMR_file$Protein.names)), ], # DNA_C5_MT
            DMR_file[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(DMR_file$short_description)), ], # DNA_C5_MT
            DMR_file[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(DMR_file$Computational_description)), ],
            DMR_file[grep("dna \\(cytosine-5\\)-methyltransferase", tolower(DMR_file$GO.biological.process)), ] # DNA_C5_MT
        ),
        Histone_Lysine_MTs = rbind(
            HLM_group,
            DMR_file[grep("SET domain", gsub("SET-domain", "SET domain", DMR_file$short_description)), ],
            DMR_file[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(DMR_file$Symbol)), ], # SDG
            DMR_file[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(DMR_file$old_symbols)), ], # SDG
            DMR_file[grep("atx[1-9]|atxr[1-9]|SDG[1-100]", tolower(DMR_file$Protein.names)), ], # SDG
            DMR_file[grep("class v-like sam-binding methyltransferase", tolower(DMR_file$Protein.families)), ] # SAM_MT
        ),
        RdDM_pathway = rddm_group,
        Royal_Family_Proteins = RF_group,
        DNA_deMTs = DNA_deMTs,
        histone_deMTs = histone_deMTs,
        REM_TFs = REM_TFs,
        chromatin_remodeling = rbind(
            DMR_file[grep("chromatin remodeling|chromatin remodeler", tolower(DMR_file$GO.biological.process)), ], # chromatin_remodeling_BP
            DMR_file[grep("chromatin remodeling|chromatin remodeler", tolower(DMR_file$Gene_description)), ], # chromatin_remodeling_pro_name
            DMR_file[grep("chromatin remodeling|chromatin remodeler", tolower(DMR_file$Computational_description)), ], # chromatin_remodeling_pro_name
            DMR_file[grep("chromatin remodeling|chromatin remodeler", tolower(DMR_file$Function)), ], # chromatin_remodeling_pro_name
            DMR_file[grep("chromatin remodeling|chromatin remodeler", tolower(DMR_file$short_description)), ] # chromatin remodeling_function
        ),
        other_methylation = rbind(
            DMR_file[grep("class i-like sam-binding methyltransferase", tolower(DMR_file$Protein.families)), ], # SAM_MT
            DMR_file[grep("class iv-like sam-binding methyltransferase", tolower(DMR_file$Protein.families)), ] # SAM_MT
            # DMR_file[grep("^2\\.1\\.1\\.", DMR_file$EC.number),] # EC_2.1.1
        ),
        primary_metabolism = primary_metabolism_v,
        secondary_metabolism = secondary_metabolism_v,
        seed_specific_genes = seed_specific_genes,
        # methionine_biosynthesis = DMR_file[grep("ath00270", DMR_file$KEGG_pathway), ],

        transporters = rbind(
            DMR_file[grep("transporter", tolower(DMR_file$Short_description)), ],
            DMR_file[grep("transporter", tolower(DMR_file$Computational_description)), ],
            DMR_file[grep("transporter", tolower(DMR_file$Function)), ],
            DMR_file[grep("transporter", tolower(DMR_file$note)), ],
            DMR_file[grep("transport", tolower(DMR_file$GO.biological.process)), ],
            DMR_file[grep("transport", tolower(DMR_file$GO.molecular.function)), ]
        ),
        response_to_stress = grep_child_stress,
        response_to_biotic_stress = grep_child_biotic,
        response_to_abiotic_stress = grep_child_abiotic
    )

    ################################################
    # # # # # # # # # # # # # # # #
    ### save csv
    final_list <- lapply(final_0_list, function(x) {
        x %>%
            dplyr::select(gene_id, Symbol, context, regionType, type, pValue, log2FC) %>%
            distinct() %>%
            arrange(pValue)
    })

    final_df <- final_list[order(sapply(final_list, nrow))] %>% bind_rows(., .id = "group")

    write.csv(final_df, paste0(out_dir, "/", cntx_file, "_", ann, "_groups_", treatment, ".csv"), row.names = F)

    ################################################################################

    plot_df <- final_df %>%
        mutate(
            p_adj_for_plot = pmax(pValue, .Machine$double.xmin),
            neglog10p = -log10(p_adj_for_plot)
        )

    plot_df$group <- gsub("_", " ", plot_df$group)
    substr(plot_df$group, 1, 1) <- toupper(substr(plot_df$group, 1, 1))

    ### volano plots
    vol_df <- plot_df %>%
        mutate(group = factor(group, levels = names(sort(table(group), decreasing = TRUE))))

    ggplot_vol_color <- if (context == "all") {
        vol_df %>% ggplot(aes(x = log2FC, y = neglog10p, color = context)) +
            scale_color_manual(values = c(CG = "#847fc5", CHG = "#00BA38", CHH = "#df6d17"))
    } else {
        vol_df %>% ggplot(aes(x = log2FC, y = neglog10p, color = regionType)) +
            scale_color_manual(values = c(gain = "#d96c6c", loss = "#6c96d9"))
    }

    vol_plot <- ggplot_vol_color +
        geom_point(alpha = 0.7, size = 1.2) +
        facet_wrap(~group, scales = "free") +
        theme_bw() +
        guides(color = guide_legend(override.aes = list(size = 4))) +
        labs(x = "log2(fold-change)", y = "-log10(pValue)", color = paste(context, "DMRs\noverlap", ann))

    img_device(paste0(out_dir, "/", cntx_file, "_", ann, "_groups_volcano_", treatment), w = 10, h = 8)
    print(vol_plot)
    dev.off()

    ### bar plots
    plot_df$context <- ifelse(context == "all", "All contexts", plot_df$context)

    return(
        plot_df %>%
            count(group, regionType, context) %>%
            group_by(group) %>%
            mutate(total = sum(n))
    )
}


################################
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
### bar plot functions
groups_barPlots <- function(x) {
    x %>%
        ggplot(aes(x = reorder(group, total), y = n, fill = regionType)) +
        geom_col() +
        geom_text(aes(y = total, label = total), hjust = -0.3, size = 3, check_overlap = TRUE) +
        facet_wrap(~context, ncol = 4, scales = "free") +
        coord_flip() +
        theme_bw() +
        scale_fill_manual(values = c(gain = "#d96c6c", loss = "#6c96d9")) +
        scale_y_continuous(expand = expansion(mult = c(0.01, 0.15))) +
        labs(x = NULL, y = "Number of DMR-annotated genes", fill = "DMR direction")
}
