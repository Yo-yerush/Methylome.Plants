read_gene_sets <- function(path) {
  if (is.null(path) || !file.exists(path)) stop("Gene-set file does not exist: ", path)
  if (grepl("\\.gmt$", path, ignore.case = TRUE)) {
    lines <- readLines(path, warn = FALSE)
    rows <- lapply(lines, function(line) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 3L) return(NULL)
      data.frame(group = fields[1], gene_id = fields[-c(1, 2)])
    })
    return(do.call(rbind, rows))
  }
  table <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  if (!all(c("group", "gene_id") %in% names(table))) {
    stop("Gene-set TSV/CSV must contain 'group' and 'gene_id' columns.")
  }
  unique(table[, c("group", "gene_id")])
}

DMRs_into_generic_groups <- function(treatment, ann, context,
                                     gene_sets_path,
                                     DMRs_ann_dir = "./genome_annotation") {
  out_dir <- file.path(DMRs_ann_dir, "functional_groups", context)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  contexts <- if (context == "all") c("CG", "CHG", "CHH") else context
  dmr_files <- file.path(
    DMRs_ann_dir,
    contexts,
    paste0("DMRs_", ann, "_", contexts, "_genom_annotations.csv")
  )
  dmr <- do.call(rbind, lapply(seq_along(dmr_files), function(i) {
    value <- utils::read.csv(dmr_files[i], check.names = FALSE)
    if (!"context" %in% names(value)) value$context <- contexts[i]
    value
  }))

  gene_sets <- read_gene_sets(gene_sets_path)
  final_df <- merge(gene_sets, dmr, by = "gene_id")
  final_df <- final_df[order(final_df$group, final_df$pValue), , drop = FALSE]
  output_context <- if (context == "all") "all_contexts" else context
  utils::write.csv(
    final_df,
    file.path(out_dir, paste0(output_context, "_", ann, "_groups_", treatment, ".csv")),
    row.names = FALSE
  )

  final_df %>%
    dplyr::count(group, regionType, context) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(total = sum(n))
}
