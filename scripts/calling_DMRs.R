calling_DMRs <- function(GRreplicates_joints, meth_var1, meth_var2,
                         var1, var2, var1_path, var2_path, comparison_name,
                         context, minProportionDiff, binSize, pValueThreshold,
                         minCytosinesCount, minReadsPerCytosine, call_ncores, is_Replicates,
                         analysis_name = "DMRs", save_csv = TRUE, chr_starts = NULL) {
  condition <- c(
    rep(var1, length(var1_path)),
    rep(var2, length(var2_path))
  )

  call_ncores_chr <- ifelse(call_ncores >= 5, 5, 1)
  call_ncores_dmr <- ifelse(call_ncores >= 10, call_ncores/5, 1)

  if (context == "CG") {
    minProportionDifference_var <- minProportionDiff[1]
  } else if (context == "CHG") {
    minProportionDifference_var <- minProportionDiff[2]
  } else if (context == "CHH") {
    minProportionDifference_var <- minProportionDiff[3]
  }
  # message(paste0(time_msg("\t"), "min difference in ",context," methylation proportion: ", minProportionDifference_var))


  # get the range for each chromosome (needed to run in parallel)
  methData_split <- split(GRreplicates_joints, seqnames(GRreplicates_joints))
  methData_split <- methData_split[order(names(methData_split))]
  start_min <- if(is.null(chr_starts)) min(start(methData_split)) else chr_starts
  end_max <- max(end(methData_split))

  if (analysis_name == "DMRs") {
    chromosome_ranges <- GRanges(seqnames = names(methData_split), IRanges(start = start_min, end = end_max))

  } else if (analysis_name == "DMVs") {
    # 200bp steps for DMVs analysis
    chromosome_ranges <- GRanges(
      seqnames = rep(names(methData_split), each = 5),
      IRanges(
        start = unlist(lapply(start_min, `+`, seq(0, 800, by = 200))),
        end = unlist(lapply(end_max - 1000, `+`, seq(0, 800, by = 200)))
      )
    )
  } else {

    stop(pate0("'analysis_name' vector must contain 'DMRs' or 'DMVs'"))
  }

  ######### run DMRcaller functions
  DMRs_gr <- GRanges()

  DMRs_list <- parallel::mclapply(
    seq_along(chromosome_ranges),
    function(i_chr) {
      tryCatch(
        {
          if (is_Replicates) {
            runReplicates <- function(cores) {
              invisible(
                capture.output(
                  x <- computeDMRsReplicates(
                    GRreplicates_joints,
                    condition = condition,
                    regions = chromosome_ranges[i_chr],
                    context = context,
                    method = "bins",
                    binSize = binSize,
                    test = "betareg",
                    pseudocountM = 1,
                    pseudocountN = 2,
                    pValueThreshold = pValueThreshold,
                    minCytosinesCount = minCytosinesCount,
                    minProportionDifference = minProportionDifference_var,
                    minGap = ifelse(analysis_name == "DMVs", 200, 0),
                    minSize = 1,
                    minReadsPerCytosine = minReadsPerCytosine,
                    cores = cores
                  )
                )
              )
              return(x)
            }
            DMRs_gr.loop <- tryCatch(
              {
                runReplicates(call_ncores_dmr)
              },
              error = function(e) {
                if (call_ncores_dmr > 10) {
                  cat("Error encountered. Retrying with cores = 10\n")
                  runReplicates(10)
                } else if (call_ncores_dmr > 1) {
                  cat("Error encountered. Retrying with cores = 1\n")
                  runReplicates(1)
                } else {
                  stop(e)
                }
              }
            )
          } else {
            # single sample in one or more of the treatments
            runPooled <- function(cores) {
              invisible(
                capture.output(
                  x <- computeDMRs(
                    meth_var1,
                    meth_var2,
                    regions = chromosome_ranges[i_chr],
                    context = context,
                    method = "bins",
                    binSize = binSize,
                    test = "fisher",
                    pValueThreshold = pValueThreshold,
                    minCytosinesCount = minCytosinesCount,
                    minProportionDifference = minProportionDifference_var,
                    minGap = ifelse(analysis_name == "DMVs", 200, 0),
                    minSize = 1,
                    minReadsPerCytosine = minReadsPerCytosine,
                    cores = cores
                  )
                )
              )
              return(x)
            }
            DMRs_gr.loop <- tryCatch(
              {
                runPooled(call_ncores_dmr)
              },
              error = function(e) {
                if (call_ncores_dmr > 10) {
                  cat("Error encountered. Retrying with cores = 10\n")
                  runPooled(10)
                } else if (call_ncores_dmr > 1) {
                  cat("Error encountered. Retrying with cores = 1\n")
                  runPooled(1)
                } else {
                  stop(e)
                }
              }
            )
          }
          return(DMRs_gr.loop)
        },
        error = function(cond) {
          message("\t* fail to calculate ", analysis_name, " in chromosome ", seqnames(chromosome_ranges)[i_chr])
          return(NULL)
        }
      )
    },
    mc.cores = call_ncores_chr
  )

  # Combine results
  DMRs_gr <- do.call(c, DMRs_list)

  if (length(DMRs_gr) != 0) {
    mcols(DMRs_gr)[, paste0("proportionsR", 1:length(condition))] <- NULL # remove 'NA' cols

    # normelized log2FC (to not get INF values)
    DMRs_gr$log2FC <- log2((DMRs_gr$proportion2 + 1e-5) / (DMRs_gr$proportion1 + 1e-5))

    DMRs_gr <- as.data.frame(DMRs_gr) %>%
      dplyr::relocate(pValue, log2FC, context, .after = strand) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T) %>%
      sort()
  }

  if (save_csv) {
    write.csv(DMRs_gr,
      paste0(analysis_name, "_", context, "_", comparison_name, ".csv"),
      row.names = F
    )
  }

  return(DMRs_gr)
}
