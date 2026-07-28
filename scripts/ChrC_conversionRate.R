conversionRate <- function(methylationData, var,
                           chloroplast_seqlevels = character(),
                           alias_map = character()) {
  if (length(chloroplast_seqlevels) == 0L) {
    message(time_msg(), var, ": skipped (no chloroplast sequence configured)")
    return(data.frame(sample = character(), conversion_rate = numeric()))
  }

  source_names <- as.character(GenomicRanges::seqnames(methylationData))
  canonical_names <- source_names
  mapped <- source_names %in% names(alias_map)
  canonical_names[mapped] <- unname(alias_map[source_names[mapped]])
  methylationData <- methylationData[canonical_names %in% chloroplast_seqlevels]

  if (length(methylationData) == 0L) {
    message(time_msg(), var, ": skipped (configured chloroplast sequence is absent)")
    return(data.frame(sample = character(), conversion_rate = numeric()))
  }

  methylation_df <- as.data.frame(S4Vectors::mcols(methylationData))
  readsM_pos <- grep("^readsM", names(methylation_df))
  rows <- lapply(readsM_pos, function(i) {
    reads_m <- methylation_df[[i]]
    reads_n <- methylation_df[[i + 1L]]
    ratio <- reads_m / reads_n
    rate <- (1 - mean(ratio[is.finite(ratio)], na.rm = TRUE)) * 100
    message(time_msg(), var, ":	", round(rate, 2), "%")
    data.frame(sample = var, conversion_rate = round(rate, 2))
  })
  do.call(rbind, rows)
}
