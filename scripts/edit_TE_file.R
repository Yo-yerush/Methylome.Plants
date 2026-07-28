# Backward-compatible adapter for already tabular TE annotations. New code
# should use read_te_annotation(), which also supports BED and GFF.
edit_TE_file <- function(TE_df,
                         seqnames_field = "seqnames",
                         start_field = "start",
                         end_field = "end",
                         strand_field = "strand",
                         id_field = "gene_id") {
  required <- c(seqnames_field, start_field, end_field)
  if (!all(required %in% names(TE_df))) {
    stop(
      "TE table must provide explicit coordinates. Missing: ",
      paste(setdiff(required, names(TE_df)), collapse = ", ")
    )
  }
  if (id_field %in% names(TE_df) && id_field != "gene_id") {
    names(TE_df)[names(TE_df) == id_field] <- "gene_id"
  }
  if (!"gene_id" %in% names(TE_df)) {
    TE_df$gene_id <- paste0("TE_", seq_len(nrow(TE_df)))
  }
  if (!strand_field %in% names(TE_df)) {
    TE_df[[strand_field]] <- '*'
  }
  GenomicRanges::makeGRangesFromDataFrame(
    TE_df,
    seqnames.field = seqnames_field,
    start.field = start_field,
    end.field = end_field,
    strand.field = strand_field,
    keep.extra.columns = TRUE
  )
}
