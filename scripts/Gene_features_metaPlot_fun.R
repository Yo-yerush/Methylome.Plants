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

Genes_features_metaPlot <- function(methylationPool_var1, methylationPool_var2, var1, var2, annotations_file, n.random, minReadsC, binSize, n.cores, promoter_upstream = 2000L) {
  
  new_path.f = "Gene_features"
  dir.create(new_path.f, showWarnings = F)
  setwd(new_path.f)
  
  # Define regions as per your provided code
  has_coding_biotype <- any(annotations_file$gene_model_type == "protein_coding", na.rm = TRUE)
  keep_biotype <- if (has_coding_biotype) annotations_file$gene_model_type == "protein_coding" else rep(TRUE, length(annotations_file))
  keep_biotype[is.na(keep_biotype)] <- FALSE
  Genes = annotations_file[which(annotations_file$type == "gene" & keep_biotype)]
  Promoters = safe_promoters(Genes, promoter_upstream)
  Promoters$type = "promoter"
  
  CDS = annotations_file[which(annotations_file$type == "CDS" & keep_biotype)]
  fiveUTR = annotations_file[which(annotations_file$type == "five_prime_UTR" & keep_biotype)]
  introns = annotations_file[which(annotations_file$type == "intron" & keep_biotype)]
  threeUTR = annotations_file[which(annotations_file$type == "three_prime_UTR" & keep_biotype)]
  
  # genes to analyze
  if (n.random[1] != "all") {
    if (length(n.random) == 1) {
    rndm_genes <- function(x) {return(sample(1:length(x), as.numeric(n.random)))} # random genes
    } else {
      rndm_genes <- function(x) {return(which(x$gene_id %in% n.random))} # list of gene IDs
    }
    
    Promoters = Promoters[rndm_genes(Promoters)]
    fiveUTR = fiveUTR[rndm_genes(fiveUTR)]
    CDS = CDS[rndm_genes(CDS)]
    introns = introns[rndm_genes(introns)]
    threeUTR = threeUTR[rndm_genes(threeUTR)]
  }
  
  regions_list = list("Promoters" = Promoters,
                      "fiveUTR" = fiveUTR,
                      "CDS" = CDS,
                      "introns" = introns,
                      "threeUTR" = threeUTR)
  
  region_names = names(regions_list)
  contexts = c("CG", "CHG", "CHH")
  
  cat(paste0("\nbin ", length(regions_list[[1]]), " protein coding gene fetures in ", binSize, "bp size and compute average methylation:\n"))

  # Function to process methylation data for each region and context
  genes_metaPlot_fun <- function(methylationData, regions_list, group_name, n.cores.f = n.cores) {
    methylationData <- methylationData[which(methylationData$readsN >= minReadsC)]
    methylationData$Proportion <- methylationData$readsM / methylationData$readsN

    methylationData_contexts <- lapply(contexts, function(cntx) methylationData[methylationData$context == cntx])
    names(methylationData_contexts) <- contexts

    mean_bins <- function(mat) {
      if (exists("meta_cpp_row_means_ignore_na", mode = "function")) {
        return(meta_cpp_row_means_ignore_na(mat))
      }
      rowMeans(mat, na.rm = TRUE)
    }

    make_equal_bins <- function(minPos, maxPos, n_bins) {
      brks <- seq(minPos, maxPos + 1, length.out = n_bins + 1)
      starts <- as.integer(floor(brks[-(n_bins + 1)]))
      ends <- as.integer(floor(brks[-1] - 1))
      list(starts = starts, ends = ends)
    }

    result_list <- list()

    for (region_name in region_names) {
      ann.obj <- regions_list[[region_name]]
      n_features <- length(ann.obj)
      region_result <- list()

      if (n_features == 0) {
        for (cntx in contexts) {
          gr_obj <- GRanges(rep(region_name, binSize), IRanges(1:binSize, 1:binSize))
          gr_obj$Proportion <- rep(NA_real_, binSize)
          region_result[[cntx]] <- gr_obj
        }
        result_list[[region_name]] <- region_result
        next
      }

      cat(paste0("\r", group_name, " >> preparing ", region_name, " bins...\t"))
      total_bins <- n_features * as.integer(binSize)
      bin_chr <- character(total_bins)
      bin_start <- integer(total_bins)
      bin_end <- integer(total_bins)
      feature_idx <- integer(total_bins)
      bin_idx <- integer(total_bins)

      for (f in seq_len(n_features)) {
        region <- ann.obj[f][, 0]
        minPos <- as.integer(start(region))
        maxPos <- as.integer(end(region))
        chr_name <- as.character(seqnames(region))
        strand_char <- as.character(strand(region))

        bins <- make_equal_bins(minPos, maxPos, as.integer(binSize))
        starts <- bins$starts
        ends <- bins$ends
        if (strand_char == "-") {
          starts <- rev(starts)
          ends <- rev(ends)
        }

        idx <- ((f - 1L) * as.integer(binSize) + 1L):(f * as.integer(binSize))
        bin_chr[idx] <- chr_name
        bin_start[idx] <- starts
        bin_end[idx] <- ends
        feature_idx[idx] <- f
        bin_idx[idx] <- seq_len(binSize)
      }
      cat("done\n")

      all_bins <- GRanges(seqnames = bin_chr, ranges = IRanges(bin_start, bin_end))

      for (cntx in contexts) {
        cat(paste0(group_name, " >> aggregating ", region_name, " (", cntx, ")...\t"))
        gr_obj <- GRanges(rep(region_name, binSize), IRanges(1:binSize, 1:binSize))
        gr_obj$Proportion <- rep(NaN, binSize)

        ctx_data <- methylationData_contexts[[cntx]]
        if (length(ctx_data) > 0) {
          feature_hits <- findOverlaps(ctx_data, ann.obj[, 0])
          valid_features <- rep(FALSE, n_features)
          if (length(feature_hits) > 0) {
            cnt <- tabulate(subjectHits(feature_hits), nbins = n_features)
            valid_features <- cnt >= 15
          }

          if (any(valid_features)) {
            hits <- findOverlaps(ctx_data, all_bins)
            if (length(hits) > 0) {
              subj <- subjectHits(hits)
              qidx <- queryHits(hits)
              keep <- valid_features[feature_idx[subj]]

              if (any(keep)) {
                means <- rep(NA_real_, length(all_bins))
                subj_keep <- subj[keep]
                qidx_keep <- qidx[keep]

                if (exists("meta_cpp_mean_by_group", mode = "function")) {
                  means <- meta_cpp_mean_by_group(subj_keep, ctx_data$Proportion[qidx_keep], length(all_bins))
                } else {
                  hit_df <- data.frame(subject_idx = subj_keep, proportion = ctx_data$Proportion[qidx_keep])
                  bin_means <- aggregate(proportion ~ subject_idx, hit_df, mean, na.rm = TRUE)
                  means[bin_means$subject_idx] <- bin_means$proportion
                }

                mat <- matrix(NA_real_, nrow = binSize, ncol = n_features)
                mat[cbind(bin_idx, feature_idx)] <- means
                gr_obj$Proportion <- mean_bins(mat)
              }
            }
          }
        }

        region_result[[cntx]] <- gr_obj
        cat("done\n")
      }

      result_list[[region_name]] <- region_result
    }

    result_list
  }
  
  ############################################
  # run main loop
  var1_metaPlot = genes_metaPlot_fun(methylationPool_var1, regions_list, var1, n.cores.f = n.cores)
  var2_metaPlot = genes_metaPlot_fun(methylationPool_var2, regions_list, var2, n.cores.f = n.cores)
  
  # Save the data
  dir.create("metaPlot_tables", showWarnings = F)
  for (region_name in region_names) {
    for (cntx in c("CG", "CHG", "CHH")) {
      try({write.csv(as.data.frame(var1_metaPlot[[region_name]][[cntx]]), paste0("metaPlot_tables/", var1, ".", cntx, ".", region_name, ".features.csv"), row.names = FALSE)})
      try({write.csv(as.data.frame(var2_metaPlot[[region_name]][[cntx]]), paste0("metaPlot_tables/", var2, ".", cntx, ".", region_name, ".features.csv"), row.names = FALSE)})
    }
  }
  
  # Plotting the data
  # bind data frames for metaPlot
  pos.end = binSize * length(region_names)
  
  #binsPlot <- function(v1.cntx.stream,v2.cntx.stream,var1,var2,cntx.m) {
  
  for (cntx in contexts) {
    
    # bind data frames by feature
    var1_proportions = c()
    var2_proportions = c()
    for (region_name in region_names) {
      var1_proportions = c(var1_proportions, var1_metaPlot[[region_name]][[cntx]]$Proportion)
      var2_proportions = c(var2_proportions, var2_metaPlot[[region_name]][[cntx]]$Proportion)
    }
    v.cntx = rbind(data.frame(pos = 1:pos.end, Proportion = var1_proportions, V = "V1"),
                   data.frame(pos = 1:pos.end, Proportion = var2_proportions, V = "V2"))
    
    # plot configuration
    min_value = min(v.cntx$Proportion)
    max_value = max(v.cntx$Proportion)
    q1_value = min_value+((max_value-min_value)/3)
    q2_value = min_value+((max_value-min_value)/3)+((max_value-min_value)/3)
    #middle_value = round(mean(c(min_value, max_value)), 2)
    q_value = min_value+((max_value-min_value)/4)
    
    legend_labels = c(paste0(" ",var1),paste0(" ",var2))
    
    region_names_plot = c("Promoter","5'UTR","CDS","Intron","3'UTR")
    main_title = "Protein Coding Genes"
    br.max = binSize * length(region_names)
    breaks_x = seq(0,br.max,by=binSize)
    breaks_vline = breaks_x[-c(1, length(breaks_x))]
    breaks_labels = breaks_x[1:length(region_names)] + binSize/2
    #breaks_and_labels <- list(breaks = seq(0,100,by=binSize), labels = region_names)
    
    plot_out = ggplot(data = v.cntx, aes(x = pos, y = Proportion, color = V, group = V)) +
      geom_vline(xintercept = breaks_vline, colour = "gray", linetype = "solid", linewidth = 0.6) +
      geom_line(linewidth = 0.65) +#, aes(linetype = V)) +
      scale_color_manual(values = c("V1" = "gray50", "V2" = "#bf6828")) +
      theme_classic() +
      labs(title = main_title, x = "", y = paste0(cntx," Methylation")) +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8, face = "bold")
      ) +
      scale_x_continuous(breaks = breaks_labels, labels = region_names_plot, minor_breaks = breaks_x, expand = expansion(add = c(0, 0))) +
      scale_y_continuous(breaks = c(min_value, q1_value, q2_value, max_value),
                         labels = c(round(min_value, 3), round(q1_value, 3),
                                    round(q2_value, 3), round(max_value, 3))) +
      annotate("text",
               x = max(breaks_vline),
               y = c(max_value*0.95, max_value*0.85),
               label = legend_labels,
               hjust = 0, vjust = 0.75, size = 3, 
               color = c("gray40","#bf6828"), fontface = "bold")
    
    
    img_device(paste0("Genes_features_",cntx,"_metaPlot_",var2,"_vs_",var1), w = 3.5, h = 2)
    #par(mar = c(2,2,1,2))
    print(plot_out)
    dev.off()
  } 
  
  # message(paste("Processed average metaPlot for", length(n.random), "Protein Coding Genes\n"))
}
