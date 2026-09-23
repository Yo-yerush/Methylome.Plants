calling_DMRs <- function(GRreplicates_joints, meth_var1, meth_var2,
                         var1, var2, var1_path, var2_path, comparison_name,
                         context, minProportionDiff, binSize, pValueThreshold,
                         minCytosinesCount, minReadsPerCytosine, call_ncores, is_Replicates,
                         analysis_name = "DMRs", save_csv = TRUE, chr_starts = NULL) {
  condition <- c(
    rep(var1, length(var1_path)),
    rep(var2, length(var2_path))
  )

  call_ncores <- max(1L, as.integer(floor(call_ncores)))

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
  call_ncores_chr <- min(5L, length(methData_split), call_ncores)
  call_ncores_dmr <- max(1L, as.integer(floor(call_ncores / call_ncores_chr)))
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

  run_with_retries <- function(run_call, chromosome_name) {
    attempts <- unique(c(call_ncores_dmr, if (call_ncores_dmr > 10L) 10L, 1L))
    for (attempt_index in seq_along(attempts)) {
      cores <- attempts[[attempt_index]]
      result <- tryCatch(
        list(value = run_call(cores), error = NULL),
        error = function(e) list(value = NULL, error = e)
      )
      if (is.null(result$error)) return(result$value)

      if (attempt_index < length(attempts)) {
        message(
          "DMR call failed on ", chromosome_name, " with cores = ", cores,
          ": ", conditionMessage(result$error),
          "; retrying with cores = ", attempts[[attempt_index + 1L]]
        )
      } else {
        stop(result$error)
      }
    }
  }

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
            DMRs_gr.loop <- run_with_retries(
              runReplicates,
              as.character(seqnames(chromosome_ranges)[i_chr])
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
            DMRs_gr.loop <- run_with_retries(
              runPooled,
              as.character(seqnames(chromosome_ranges)[i_chr])
            )
          }
          return(DMRs_gr.loop)
        },
        error = function(cond) {
          message(
            "\t* fail to calculate ", analysis_name, " in chromosome ",
            seqnames(chromosome_ranges)[i_chr], ": ", conditionMessage(cond)
          )
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

#################################################################

calling_DMRs_queue <- function(GRreplicates_joints, meth_var1, meth_var2,
                               var1, var2, var1_path, var2_path, comparison_name,
                               contexts, minProportionDiff, binSize, pValueThreshold,
                               minCytosinesCount, minReadsPerCytosine, call_ncores,
                               is_Replicates) {
  call_ncores <- max(1L, as.integer(floor(call_ncores)))
  contexts <- as.character(contexts)
  condition <- c(rep(var1, length(var1_path)), rep(var2, length(var2_path)))
  min_difference_by_context <- stats::setNames(
    as.numeric(minProportionDiff),
    c("CG", "CHG", "CHH")[seq_along(minProportionDiff)]
  )
  if (!all(contexts %in% names(min_difference_by_context))) {
    stop("minProportionDiff must provide thresholds for every requested context")
  }

  chromosomes <- as.character(seqlevelsInUse(GRreplicates_joints))
  chromosomes <- chromosomes[chromosomes %in% unique(as.character(seqnames(GRreplicates_joints)))]
  chromosome_ranges <- range(GRreplicates_joints, ignore.strand = TRUE)
  chromosome_ranges <- chromosome_ranges[match(chromosomes, as.character(seqnames(chromosome_ranges)))]

  tasks <- expand.grid(
    context = contexts,
    chromosome = chromosomes,
    stringsAsFactors = FALSE
  )
  tasks$key <- paste(tasks$context, tasks$chromosome, sep = "::")

  split_for_tasks <- function(methylation_data) {
    group <- factor(
      paste(as.character(methylation_data$context), as.character(seqnames(methylation_data)), sep = "::"),
      levels = tasks$key
    )
    split(methylation_data, group, drop = FALSE)
  }

  message(time_msg(), "pre-splitting methylation data into context/chromosome tasks")
  replicate_tasks <- split_for_tasks(GRreplicates_joints)
  tasks$sites <- lengths(replicate_tasks)
  tasks <- tasks[tasks$sites > 0L, , drop = FALSE]
  replicate_tasks <- replicate_tasks[tasks$key]

  if (!is_Replicates) {
    meth_var1_tasks <- split_for_tasks(meth_var1)[tasks$key]
    meth_var2_tasks <- split_for_tasks(meth_var2)[tasks$key]
  }

  task_results <- vector("list", nrow(tasks))
  pending <- seq_len(nrow(tasks))

  # Long tasks start first; mc.preschedule=FALSE lets every free worker take the
  # next chromosome/context rather than waiting for a fixed context batch.
  pending <- pending[order(tasks$sites[pending], decreasing = TRUE)]
  task_count <- length(pending)
  workers <- min(call_ncores, task_count)

  if (length(pending)) {
    message(
      time_msg(),
      "DMR worker queue: ", task_count,
      " tasks, ", workers, " workers, 1 DMRcaller core per task"
    )

    computed <- parallel::mclapply(
      pending,
      function(i) {
        context <- tasks$context[[i]]
        chromosome <- tasks$chromosome[[i]]
        min_difference <- min_difference_by_context[[context]]
        region <- chromosome_ranges[match(chromosome, as.character(seqnames(chromosome_ranges)))]

        result <- tryCatch({
          if (is_Replicates) {
            invisible(capture.output(
              value <- computeDMRsReplicates(
                replicate_tasks[[tasks$key[[i]]]],
                condition = condition,
                regions = region,
                context = context,
                method = "bins",
                binSize = binSize,
                test = "betareg",
                pseudocountM = 1,
                pseudocountN = 2,
                pValueThreshold = pValueThreshold,
                minCytosinesCount = minCytosinesCount,
                minProportionDifference = min_difference,
                minGap = 0,
                minSize = 1,
                minReadsPerCytosine = minReadsPerCytosine,
                cores = 1L
              )
            ))
          } else {
            invisible(capture.output(
              value <- computeDMRs(
                meth_var1_tasks[[tasks$key[[i]]]],
                meth_var2_tasks[[tasks$key[[i]]]],
                regions = region,
                context = context,
                method = "bins",
                binSize = binSize,
                test = "fisher",
                pValueThreshold = pValueThreshold,
                minCytosinesCount = minCytosinesCount,
                minProportionDifference = min_difference,
                minGap = 0,
                minSize = 1,
                minReadsPerCytosine = minReadsPerCytosine,
                cores = 1L
              )
            ))
          }

          list(ok = TRUE, index = i, value = value)
        }, error = function(error) {
          list(ok = FALSE, index = i, error = conditionMessage(error))
        })
        result
      },
      mc.cores = workers,
      mc.preschedule = FALSE
    )

    computed <- lapply(seq_along(computed), function(result_index) {
      result <- computed[[result_index]]
      if (is.list(result) && !is.null(result$ok)) return(result)
      list(
        ok = FALSE,
        index = pending[[result_index]],
        error = paste0("worker terminated without a result: ", paste(as.character(result), collapse = " "))
      )
    })

    failures <- vapply(computed, function(x) !isTRUE(x$ok), logical(1))
    for (result in computed[!failures]) task_results[[result$index]] <- result$value
    if (any(failures)) {
      failed_messages <- vapply(
        computed[failures],
        function(x) paste0(tasks$context[[x$index]], "/", tasks$chromosome[[x$index]], ": ", x$error),
        character(1)
      )
      stop(
        "DMR queue failed for ", sum(failures), " task(s). ",
        paste(failed_messages, collapse = "; ")
      )
    }
  }

  DMRs_results <- lapply(contexts, function(context) {
    context_results <- task_results[tasks$context == context]
    context_results <- context_results[vapply(context_results, function(x) !is.null(x), logical(1))]
    DMRs_gr <- if (length(context_results)) do.call(c, context_results) else GRanges()

    if (length(DMRs_gr)) {
      proportion_columns <- grep("^proportionsR", names(mcols(DMRs_gr)), value = TRUE)
      if (length(proportion_columns)) mcols(DMRs_gr)[, proportion_columns] <- NULL
      DMRs_gr$log2FC <- log2((DMRs_gr$proportion2 + 1e-5) / (DMRs_gr$proportion1 + 1e-5))
      DMRs_gr <- as.data.frame(DMRs_gr) %>%
        dplyr::relocate(pValue, log2FC, context, .after = strand) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
        sort()
    }

    write.csv(
      DMRs_gr,
      paste0("DMRs_", context, "_", comparison_name, ".csv"),
      row.names = FALSE
    )
    DMRs_gr
  })
  names(DMRs_results) <- contexts
  DMRs_results
}
