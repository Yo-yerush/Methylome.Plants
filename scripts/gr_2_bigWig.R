gr_2_bigWig <- function(gr, out_file_name = "gr_2_bigWig.bw",
                        out_format = "bigWig", reference_bundle = NULL) {
  if (is.null(gr) || length(gr) == 0L) return(invisible(FALSE))

  gr_out <- gr
  S4Vectors::mcols(gr_out) <- NULL
  score <- if ("log2FC" %in% names(S4Vectors::mcols(gr))) {
    gr$log2FC
  } else if ("score" %in% names(S4Vectors::mcols(gr))) {
    gr$score
  } else {
    rep(0, length(gr))
  }
  gr_out$score <- round(as.numeric(score), 2)

  configured_sizes <- bundle_chromosome_sizes(reference_bundle)
  present <- as.character(GenomeInfoDb::seqlevels(gr_out))
  if (!is.null(configured_sizes)) {
    missing <- setdiff(present, names(configured_sizes))
    if (length(missing)) {
      stop("No chromosome size configured for bigWig sequence levels: ",
           paste(missing, collapse = ", "))
    }
    GenomeInfoDb::seqlengths(gr_out) <- configured_sizes[present]
  } else {
    existing <- GenomeInfoDb::seqlengths(gr_out)
    if (any(is.na(existing))) {
      inferred <- vapply(
        split(GenomicRanges::end(gr_out), GenomicRanges::seqnames(gr_out)),
        max,
        numeric(1)
      )
      existing[names(inferred)] <- inferred
      GenomeInfoDb::seqlengths(gr_out) <- existing
      warning("bigWig sequence lengths were inferred from observed data; provide chromosome_sizes for exact bounds.")
    }
  }

  rtracklayer::export(gr_out, out_file_name, format = out_format)
  invisible(TRUE)
}
