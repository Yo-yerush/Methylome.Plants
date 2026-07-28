harmonize_seqlevels <- function(gr_obj, reference_bundle = NULL, keep_non_primary = NULL) {
  if (is.null(gr_obj) || length(gr_obj) == 0) return(gr_obj)

  aliases <- bundle_alias_map(reference_bundle)
  current <- as.character(GenomeInfoDb::seqlevels(gr_obj))
  renamed <- current
  has_alias <- current %in% names(aliases)
  renamed[has_alias] <- unname(aliases[current[has_alias]])

  if (anyDuplicated(renamed)) {
    conflicting <- unique(renamed[duplicated(renamed)])
    stop(
      "Sequence alias map collapses multiple sequence levels onto: ",
      paste(conflicting, collapse = ", "),
      ". Normalize the source files before running."
    )
  }
  if (!identical(current, renamed)) {
    gr_obj <- GenomeInfoDb::renameSeqlevels(gr_obj, stats::setNames(renamed, current))
  }

  metadata_aliases <- c(
    biotype = "gene_model_type",
    gene_biotype = "gene_model_type",
    locus_tag = "gene_id"
  )
  for (source in names(metadata_aliases)) {
    target <- unname(metadata_aliases[[source]])
    if (source %in% names(S4Vectors::mcols(gr_obj)) &&
        !target %in% names(S4Vectors::mcols(gr_obj))) {
      S4Vectors::mcols(gr_obj)[[target]] <- S4Vectors::mcols(gr_obj)[[source]]
    }
  }

  primary <- as.character(bundle_get(reference_bundle, "genome.primary_seqlevels", character()))
  if (is.null(keep_non_primary)) {
    keep_non_primary <- isTRUE(bundle_get(reference_bundle, "genome.include_non_primary", TRUE))
  }
  if (length(primary) > 0 && !keep_non_primary) {
    present_primary <- intersect(primary, as.character(GenomeInfoDb::seqlevels(gr_obj)))
    if (length(present_primary) == 0L) {
      stop('None of the configured primary_seqlevels occur in this genomic input.')
    }
    dropped <- setdiff(as.character(GenomeInfoDb::seqlevels(gr_obj)), present_primary)
    if (length(dropped)) {
      message('Reference bundle excludes non-primary sequence levels: ', paste(dropped, collapse = ', '))
    }
    gr_obj <- GenomeInfoDb::keepSeqlevels(gr_obj, present_primary, pruning.mode = "coarse")
  }

  sizes <- bundle_chromosome_sizes(reference_bundle)
  if (!is.null(sizes)) {
    present <- intersect(names(sizes), as.character(GenomeInfoDb::seqlevels(gr_obj)))
    current_lengths <- GenomeInfoDb::seqlengths(gr_obj)
    current_lengths[present] <- sizes[present]
    GenomeInfoDb::seqlengths(gr_obj) <- current_lengths
  }

  desired_order <- c(
    primary,
    as.character(bundle_get(reference_bundle, "genome.chloroplast_seqlevels", character())),
    as.character(bundle_get(reference_bundle, "genome.mitochondrial_seqlevels", character()))
  )
  desired_order <- unique(c(
    intersect(desired_order, as.character(GenomeInfoDb::seqlevels(gr_obj))),
    setdiff(as.character(GenomeInfoDb::seqlevels(gr_obj)), desired_order)
  ))
  GenomeInfoDb::seqlevels(gr_obj) <- desired_order
  GenomicRanges::sort(gr_obj, by = ~seqnames + start)
}

# Backward-compatible name. Without a bundle this preserves source sequence
# names instead of guessing an organism and deleting non-numeric sequences.
trimm_and_rename <- function(gr_obj, reference_bundle = NULL) {
  harmonize_seqlevels(gr_obj, reference_bundle = reference_bundle)
}
