if (!exists(".getMetaPlotScriptDir", mode = "function")) {
  .getMetaPlotScriptDir <- function() {
    frames <- sys.frames()
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile)) {
        return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
      }
    }
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
}

if (!exists(".loadMetaPlotsCpp", mode = "function")) {
  .loadMetaPlotsCpp <- local({
    loaded <- FALSE
    tried <- FALSE
    function() {
      if (loaded) {
        return(invisible(TRUE))
      }
      if (tried) {
        return(invisible(FALSE))
      }
      tried <<- TRUE

      if (!requireNamespace("Rcpp", quietly = TRUE)) {
        return(invisible(FALSE))
      }

      cppFile <- file.path(.getMetaPlotScriptDir(), "metaPlots_rcpp.cpp")
      if (!file.exists(cppFile)) {
        return(invisible(FALSE))
      }

      ok <- tryCatch({
        Rcpp::sourceCpp(cppFile, env = parent.frame())
        TRUE
      }, error = function(e) {
        message("metaPlots Rcpp acceleration disabled: ", conditionMessage(e))
        FALSE
      })

      loaded <<- isTRUE(ok)
      invisible(loaded)
    }
  })
}

.loadMetaPlotsCpp()

Genes_metaPlot <- function(methylationPool_var1,methylationPool_var2,var1,var2,annotations_file,n.random,minReadsC,n.cores,is_TE=F) {

  if (is_TE) {
    new_path.f = "TEs"
    coding.Genes = annotations_file
  } else {
    new_path.f = "Genes"
    coding.Genes.0 <- annotations_file[which(annotations_file$type == "gene")]
    coding_mask <- coding.Genes.0$gene_model_type == "protein_coding"
    coding.Genes <- if (any(coding_mask, na.rm = TRUE)) coding.Genes.0[coding_mask %in% TRUE] else coding.Genes.0
  }
  
  if (n.random[1] != "all") {
    if (length(n.random) == 1) {
      rndm_genes = sample(1:length(coding.Genes), as.numeric(n.random)) # random genes
      coding.Genes = coding.Genes[rndm_genes]
    } else {
      coding.Genes = coding.Genes[which(coding.Genes$gene_id %in% n.random)] # list of gene IDs
    }
  }
  
  dir.create(new_path.f, showWarnings = F)
  setwd(new_path.f)
  
  cat("\nbin", length(coding.Genes), new_path.f, "body ±2kb in 20bp size and compute average methylation:\n")
  
  # make windowSize ranges with the average value
  genes_metaPlot_fun <- function(methylationData, ann.obj, group_name, n.cores.f = n.cores) {
    methylationData <- methylationData[which(methylationData$readsN >= minReadsC)]
    methylationData$Proportion <- methylationData$readsM / methylationData$readsN

    n_genes <- length(ann.obj)
    if (n_genes == 0) {
      mk_empty <- function(nm) {
        gr <- GRanges(rep(nm, 20), IRanges(1:20, 1:20))
        gr$Proportion <- rep(NA_real_, 20)
        gr
      }
      empty_list <- GRangesList("up.stream" = mk_empty("up.stream"), "gene.body" = mk_empty("gene.body"), "down.stream" = mk_empty("down.stream"))
      return(list(gr_list_CG = empty_list, gr_list_CHG = empty_list, gr_list_CHH = empty_list))
    }

    total_bins <- n_genes * 60L
    bin_chr <- character(total_bins)
    bin_start <- integer(total_bins)
    bin_end <- integer(total_bins)
    bin_gene <- integer(total_bins)
    bin_idx <- integer(total_bins)

    make_equal_bins <- function(minPos, maxPos, n_bins) {
      brks <- seq(minPos, maxPos + 1, length.out = n_bins + 1)
      starts <- as.integer(floor(brks[-(n_bins + 1)]))
      ends <- as.integer(floor(brks[-1] - 1))
      list(starts = starts, ends = ends)
    }

    cat(paste0("\r", group_name, " >> preparing ", n_genes, " ", gsub("s","",new_path.f), " bin maps...\t[0%]     "))
    for (g in seq_len(n_genes)) {
      region <- ann.obj[g][, 0]
      minPos <- as.integer(start(region))
      maxPos <- as.integer(end(region))
      chr_name <- as.character(seqnames(region))
      strand_char <- as.character(strand(region))

      up_starts <- as.integer(seq(minPos - 2000L, minPos - 100L, by = 100L))
      up_ends <- up_starts + 99L
      body_bins <- make_equal_bins(minPos, maxPos, 20L)
      body_starts <- body_bins$starts
      body_ends <- body_bins$ends
      down_starts <- as.integer(seq(maxPos + 1L, maxPos + 1901L, by = 100L))
      down_ends <- down_starts + 99L

      starts60 <- c(up_starts, body_starts, down_starts)
      ends60 <- c(up_ends, body_ends, down_ends)
      if (strand_char == "-") {
        starts60 <- c(rev(down_starts), rev(body_starts), rev(up_starts))
        ends60 <- c(rev(down_ends), rev(body_ends), rev(up_ends))
      }

      idx <- ((g - 1L) * 60L + 1L):(g * 60L)
      bin_chr[idx] <- chr_name
      bin_start[idx] <- starts60
      bin_end[idx] <- ends60
      bin_gene[idx] <- g
      bin_idx[idx] <- 1:60

      if (g %% 100 == 0 || g == n_genes) {
        cat(paste0("\r", group_name, " >> preparing ", n_genes, " ", gsub("s","",new_path.f), " bin maps...\t[", round(g / n_genes * 100, 0), "%]     "))
      }
    }
    cat("\n")

    all_bins <- GRanges(seqnames = bin_chr, ranges = IRanges(bin_start, bin_end))
    mean_bins <- function(mat) {
      if (exists("meta_cpp_row_means_ignore_na", mode = "function")) {
        return(meta_cpp_row_means_ignore_na(mat))
      }
      rowMeans(mat, na.rm = TRUE)
    }

    summarize_context <- function(cntx) {
      cat(paste0(group_name, " >> aggregating context ", cntx, "...\t"))
      ctx_data <- methylationData[methylationData$context == cntx]
      mat <- matrix(NA_real_, nrow = 60, ncol = n_genes)

      if (length(ctx_data) > 0) {
        hits <- findOverlaps(ctx_data, all_bins)
        if (length(hits) > 0) {
          subj <- subjectHits(hits)
          qidx <- queryHits(hits)
          means <- rep(NA_real_, length(all_bins))
          if (exists("meta_cpp_mean_by_group", mode = "function")) {
            means <- meta_cpp_mean_by_group(subj, ctx_data$Proportion[qidx], length(all_bins))
          } else {
            hit_df <- data.frame(subject_idx = subj, proportion = ctx_data$Proportion[qidx])
            bin_means <- aggregate(proportion ~ subject_idx, hit_df, mean, na.rm = TRUE)
            means[bin_means$subject_idx] <- bin_means$proportion
          }
          mat[cbind(bin_idx, bin_gene)] <- means
        }
      }

      avg60 <- mean_bins(mat)
      cat("done\n")
      list(
        up = avg60[1:20],
        body = avg60[21:40],
        down = avg60[41:60]
      )
    }

    ctx_results <- lapply(c("CG", "CHG", "CHH"), summarize_context)
    names(ctx_results) <- c("CG", "CHG", "CHH")

    make_stream_gr <- function(name, vals) {
      gr <- GRanges(rep(name, 20), IRanges(1:20, 1:20))
      gr$Proportion <- vals
      gr
    }

    gr_list_CG <- GRangesList(
      "up.stream" = make_stream_gr("up.stream", ctx_results$CG$up),
      "gene.body" = make_stream_gr("gene.body", ctx_results$CG$body),
      "down.stream" = make_stream_gr("down.stream", ctx_results$CG$down)
    )
    gr_list_CHG <- GRangesList(
      "up.stream" = make_stream_gr("up.stream", ctx_results$CHG$up),
      "gene.body" = make_stream_gr("gene.body", ctx_results$CHG$body),
      "down.stream" = make_stream_gr("down.stream", ctx_results$CHG$down)
    )
    gr_list_CHH <- GRangesList(
      "up.stream" = make_stream_gr("up.stream", ctx_results$CHH$up),
      "gene.body" = make_stream_gr("gene.body", ctx_results$CHH$body),
      "down.stream" = make_stream_gr("down.stream", ctx_results$CHH$down)
    )

    list(
      gr_list_CG = gr_list_CG,
      gr_list_CHG = gr_list_CHG,
      gr_list_CHH = gr_list_CHH
    )
  }
  
  ############################################
  # run main loop
  var1_metaPlot = genes_metaPlot_fun(methylationPool_var1, coding.Genes, var1)
  var2_metaPlot = genes_metaPlot_fun(methylationPool_var2, coding.Genes, var2)

  v1.CG = var1_metaPlot$gr_list_CG
  v1.CHG = var1_metaPlot$gr_list_CHG
  v1.CHH = var1_metaPlot$gr_list_CHH
  
  v2.CG = var2_metaPlot$gr_list_CG
  v2.CHG = var2_metaPlot$gr_list_CHG
  v2.CHH = var2_metaPlot$gr_list_CHH
  
  dir.create("metaPlot_tables", showWarnings = F)
  for (name in c("up.stream","gene.body","down.stream")) {
    try({write.csv(v1.CG[[name]], paste0("metaPlot_tables/",var1,".CG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v1.CHG[[name]], paste0("metaPlot_tables/",var1,".CHG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v1.CHH[[name]], paste0("metaPlot_tables/",var1,".CHH.",name,".csv"), row.names = FALSE)})
    try({write.csv(v2.CG[[name]], paste0("metaPlot_tables/",var2,".CG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v2.CHG[[name]], paste0("metaPlot_tables/",var2,".CHG.",name,".csv"), row.names = FALSE)})
    try({write.csv(v2.CHH[[name]], paste0("metaPlot_tables/",var2,".CHH.",name,".csv"), row.names = FALSE)})
  }
  
  # bind data frames for metaPlot
  v1.CG.bind = data.frame(pos = 1:60, Proportion = c(v1.CG[[1]],v1.CG[[2]],v1.CG[[3]])$Proportion)
  v1.CHG.bind = data.frame(pos = 1:60, Proportion = c(v1.CHG[[1]],v1.CHG[[2]],v1.CHG[[3]])$Proportion)
  v1.CHH.bind = data.frame(pos = 1:60, Proportion = c(v1.CHH[[1]],v1.CHH[[2]],v1.CHH[[3]])$Proportion)
  
  v2.CG.bind = data.frame(pos = 1:60, Proportion = c(v2.CG[[1]],v2.CG[[2]],v2.CG[[3]])$Proportion)
  v2.CHG.bind = data.frame(pos = 1:60, Proportion = c(v2.CHG[[1]],v2.CHG[[2]],v2.CHG[[3]])$Proportion)
  v2.CHH.bind = data.frame(pos = 1:60, Proportion = c(v2.CHH[[1]],v2.CHH[[2]],v2.CHH[[3]])$Proportion)
  
  binsPlot <- function(v1.cntx.stream,v2.cntx.stream,var1,var2,cntx.m) {
    
    v1.cntx.stream$V = "V1"
    v2.cntx.stream$V = "V2"
    v.cntx.stream = rbind(v1.cntx.stream,v2.cntx.stream)
    
    min_value = min(v.cntx.stream$Proportion)
    max_value = max(v.cntx.stream$Proportion)
    q1_value = min_value+((max_value-min_value)/3)
    q2_value = min_value+((max_value-min_value)/3)+((max_value-min_value)/3)
    #middle_value = round(mean(c(min_value, max_value)), 2)
    q_value = min_value+((max_value-min_value)/4)
    
    if (is_TE) {
      legend_labels = c(var1,paste0("\n",var2))
    } else {
      legend_labels = c(paste0("\n",var1),paste0("\n\n",var2))
    }

    main_title = ifelse(is_TE, "TEs", "Gene bodies")
    breaks_and_labels <- list(breaks = c(1.35, 20, 40, 59.65), labels = c("  -2kb", "TSS", "TTS", "+2kb   "))
    
    plot_out = ggplot(data = v.cntx.stream, aes(x = pos, y = Proportion, color = V, group = V)) +
      geom_vline(xintercept = c(20, 40), colour = "gray", linetype = "solid", linewidth = 0.5) +
      geom_line(linewidth = 0.65) +#, aes(linetype = V)) +
      #scale_color_manual(values = c("V1" = "#e37320", "V2" = "#0072B2")) +
      scale_color_manual(values = c("V1" = "gray50", "V2" = "#bf6828")) +
      #scale_linetype_manual(values = c("V1" = "dashed", "V2" = "solid")) +
      theme_classic() +
      labs(title = main_title,
           x = "",
           y = paste0(cntx.m," methylation")) +
      theme(legend.position = "none",
            axis.line.x = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.75),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9)
      ) +
      scale_x_continuous(breaks = breaks_and_labels$breaks, labels = breaks_and_labels$labels, expand = expansion(add = c(0,0))) +
      scale_y_continuous(breaks = c(min_value, q1_value, q2_value, max_value),
                         labels = c(round(min_value, 3), round(q1_value, 3),
                                    round(q2_value, 3), round(max_value, 3))) +
      annotate("text",
               x = 3,
               y = ifelse(is_TE,max_value,q1_value),
               label = legend_labels,
               hjust = 0, vjust = 0.75, size = 3, 
               color = c("gray40","#bf6828"), fontface = "bold")
    
    
    if (is_TE) {
      img_device(paste0("TEs_",cntx.m,"_metaPlot_",var2,"_vs_",var1), w = 1.88, h = 1.94)
    } else {
      img_device(paste0("Genes_",cntx.m,"_metaPlot_",var2,"_vs_",var1), w = 1.88, h = 1.94)
    }
    #par(mar = c(2,2,1,2))
    print(plot_out)
    dev.off()
  }  
  
  binsPlot(v1.CG.bind,v2.CG.bind,var1,var2,"CG")
  binsPlot(v1.CHG.bind,v2.CHG.bind,var1,var2,"CHG")
  binsPlot(v1.CHH.bind,v2.CHH.bind,var1,var2,"CHH")
  
  # message(paste("process average metaPlot to",
  #                length(coding.Genes),
  #                ifelse(!is_TE,"Protein Coding Genes","Transposable Elements\n")))
}
