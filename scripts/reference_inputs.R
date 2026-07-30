.annotation_path_format <- function(path, configured = NULL) {
  if (!is.null(configured)) return(tolower(configured))
  clean <- sub("\\.gz$", "", tolower(path))
  sub("^.*\\.", "", clean)
}

.rename_metadata_field <- function(gr, source, target) {
  if (!is.null(source) && source %in% names(S4Vectors::mcols(gr)) && source != target) {
    if (!target %in% names(S4Vectors::mcols(gr))) {
      S4Vectors::mcols(gr)[[target]] <- S4Vectors::mcols(gr)[[source]]
    }
  }
  gr
}

.metadata_first_value <- function(value) {
  if (inherits(value, 'List')) {
    return(vapply(as.list(value), function(x) {
      if (length(x)) as.character(x[[1]]) else NA_character_
    }, character(1)))
  }
  value <- as.character(value)
  sub(',.*$', '', value)
}

.resolve_gff_gene_ids <- function(gr) {
  metadata <- S4Vectors::mcols(gr)
  if (!all(c('ID', 'Parent', 'type') %in% names(metadata))) return(rep(NA_character_, length(gr)))

  ids <- .metadata_first_value(metadata[['ID']])
  parents <- .metadata_first_value(metadata[['Parent']])
  feature_types <- tolower(as.character(metadata[['type']]))
  gene_rows <- grepl('gene$', feature_types)

  id_to_gene <- stats::setNames(rep(NA_character_, length(ids)), ids)
  valid_roots <- gene_rows & !is.na(ids) & nzchar(ids)
  id_to_gene[ids[valid_roots]] <- ids[valid_roots]

  # Resolve transcript and child-feature ancestry one level at a time.
  for (iteration in seq_len(20L)) {
    unresolved <- !is.na(ids) & nzchar(ids) & is.na(id_to_gene[ids]) &
      !is.na(parents) & nzchar(parents)
    if (!any(unresolved)) break
    parent_gene <- unname(id_to_gene[parents[unresolved]])
    resolved <- !is.na(parent_gene)
    if (!any(resolved)) break
    child_ids <- ids[unresolved][resolved]
    id_to_gene[child_ids] <- parent_gene[resolved]
  }

  result <- unname(id_to_gene[ids])
  parent_gene <- unname(id_to_gene[parents])
  use_parent <- !gene_rows & !is.na(parent_gene)
  result[use_parent] <- parent_gene[use_parent]
  result
}

read_gene_annotation <- function(path = NULL, bundle = NULL) {
  path <- path %||% bundle_get(bundle, "annotation.genes")
  if (is.null(path) || !nzchar(path)) return(GenomicRanges::GRanges())
  if (!file.exists(path)) stop("Gene annotation does not exist: ", path)

  format <- .annotation_path_format(path, bundle_get(bundle, "annotation.format"))
  if (format == "csv") {
    annotation <- utils::read.csv(path, check.names = FALSE)
    gr <- GenomicRanges::makeGRangesFromDataFrame(annotation, keep.extra.columns = TRUE)
  } else if (format %in% c("gff", "gff3", "gtf")) {
    gr <- rtracklayer::import(path)
  } else {
    stop("Unsupported gene annotation format: ", format)
  }

  gr <- .rename_metadata_field(
    gr,
    bundle_get(bundle, "annotation.gene_id_field", "gene_id"),
    "gene_id"
  )
  gr <- .rename_metadata_field(
    gr,
    bundle_get(bundle, "annotation.biotype_field", "gene_model_type"),
    "gene_model_type"
  )
  gr <- .rename_metadata_field(
    gr,
    bundle_get(bundle, "annotation.feature_type_field", "type"),
    "type"
  )
  gr <- .rename_metadata_field(
    gr,
    bundle_get(bundle, "annotation.transcript_id_field", "transcript_id"),
    "transcript_id"
  )

  if (format %in% c('gff', 'gff3')) {
    resolved_gene_ids <- .resolve_gff_gene_ids(gr)
    if (!'gene_id' %in% names(S4Vectors::mcols(gr))) gr$gene_id <- NA_character_
    resolved <- !is.na(resolved_gene_ids) & nzchar(resolved_gene_ids)
    gr$gene_id[resolved] <- resolved_gene_ids[resolved]
  }

  if (!"type" %in% names(S4Vectors::mcols(gr))) {
    stop("Gene annotation has no feature type column after field mapping.")
  }
  if (!"gene_id" %in% names(S4Vectors::mcols(gr))) {
    candidates <- intersect(c("ID", "Name", "locus_tag", "Parent"), names(S4Vectors::mcols(gr)))
    if (length(candidates) == 0) {
      stop("Gene annotation has no gene identifier column after field mapping.")
    }
    gr$gene_id <- as.character(S4Vectors::mcols(gr)[[candidates[1]]])
  }
  if (!"gene_model_type" %in% names(S4Vectors::mcols(gr))) {
    gr$gene_model_type <- NA_character_
  }
  coding_values <- as.character(bundle_get(bundle, 'annotation.protein_coding_values', 'protein_coding'))
  is_coding <- as.character(gr$gene_model_type) %in% coding_values
  gr$gene_model_type[is_coding] <- 'protein_coding'

  feature_map <- bundle_get(bundle, "annotation.feature_type_map")
  if (!is.null(feature_map)) {
    replacements <- unlist(feature_map, use.names = TRUE)
    matched <- replacements[as.character(gr$type)]
    gr$type[!is.na(matched)] <- matched[!is.na(matched)]
  }
  gr
}

read_te_annotation <- function(path = NULL, bundle = NULL) {
  path <- path %||% bundle_get(bundle, "annotation.transposable_elements")
  if (is.null(path) || !nzchar(path)) return(GenomicRanges::GRanges())
  if (!file.exists(path)) stop("TE annotation does not exist: ", path)

  format <- .annotation_path_format(path, bundle_get(bundle, "annotation.te_format"))
  if (format == "tair_legacy") {
    te <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c(
      "Transposon_Name", "orientation_is_5prime",
      "Transposon_min_Start", "Transposon_max_End"
    )
    if (!all(required %in% names(te))) {
      stop("TAIR legacy TE file is missing required columns: ", path)
    }
    chromosome <- sub("^AT([0-9]+)TE.*$", "Chr\\1", te$Transposon_Name)
    invalid <- chromosome == te$Transposon_Name
    if (any(invalid)) stop("Could not derive chromosome for some TAIR TE identifiers.")
    te$seqnames <- chromosome
    te$start <- as.integer(te$Transposon_min_Start)
    te$end <- as.integer(te$Transposon_max_End)
    te$strand <- ifelse(tolower(te$orientation_is_5prime) == "true", "+", "-")
    te$gene_id <- te$Transposon_Name
    return(GenomicRanges::makeGRangesFromDataFrame(te, keep.extra.columns = TRUE))
  }

  if (format %in% c("gff", "gff3", "gtf", "bed")) {
    gr <- rtracklayer::import(path)
  } else if (format %in% c("csv", "tsv", "txt")) {
    sep <- if (format == "csv") "," else "\t"
    te <- utils::read.table(
      path,
      header = TRUE,
      sep = sep,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    gr <- GenomicRanges::makeGRangesFromDataFrame(te, keep.extra.columns = TRUE)
  } else {
    stop("Unsupported TE annotation format: ", format)
  }

  field_map <- list(
    gene_id = bundle_get(bundle, "annotation.te_id_field", "gene_id"),
    Transposon_Family = bundle_get(bundle, "annotation.te_family_field", "Transposon_Family"),
    Transposon_Super_Family = bundle_get(bundle, "annotation.te_superfamily_field", "Transposon_Super_Family")
  )
  for (target in names(field_map)) {
    gr <- .rename_metadata_field(gr, field_map[[target]], target)
  }
  if (!"gene_id" %in% names(S4Vectors::mcols(gr))) {
    candidates <- intersect(c("ID", "Name", "name", "repeat_id", "te_id"), names(S4Vectors::mcols(gr)))
    if (length(candidates)) gr$gene_id <- as.character(S4Vectors::mcols(gr)[[candidates[1]]])
  }
  if (!"gene_id" %in% names(S4Vectors::mcols(gr))) {
    gr$gene_id <- paste0("TE_", seq_along(gr))
  }
  if (!"Transposon_Family" %in% names(S4Vectors::mcols(gr))) {
    gr$Transposon_Family <- NA_character_
  }
  if (!"Transposon_Super_Family" %in% names(S4Vectors::mcols(gr))) {
    gr$Transposon_Super_Family <- 'unclassified'
  }
  gr
}

read_gene_descriptions <- function(path = NULL, bundle = NULL) {
  path <- path %||% bundle_get(bundle, "annotation.descriptions")
  if (is.null(path) || !nzchar(path)) return(data.frame(gene_id = character()))
  if (!file.exists(path)) stop("Gene description table does not exist: ", path)
  sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  descriptions <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  id_field <- bundle_get(bundle, "annotation.description_gene_id_field", names(descriptions)[1])
  if (!id_field %in% names(descriptions)) {
    stop("Description gene ID field is absent: ", id_field)
  }
  names(descriptions)[names(descriptions) == id_field] <- "gene_id"
  descriptions
}

validate_reference_inputs <- function(methylation, annotation = NULL, tes = NULL, bundle = NULL) {
  objects <- list(methylation = methylation, annotation = annotation, transposable_elements = tes)
  objects <- objects[vapply(objects, function(x) !is.null(x) && length(x) > 0, logical(1))]
  if (length(objects) < 2) return(invisible(TRUE))

  meth_levels <- as.character(GenomeInfoDb::seqlevels(methylation))
  for (name in names(objects)[-1]) {
    other <- as.character(GenomeInfoDb::seqlevels(objects[[name]]))
    shared <- intersect(meth_levels, other)
    if (length(shared) == 0) {
      stop("No sequence names are shared between methylation and ", name, ". Check the alias map and assembly.")
    }
  }

  sizes <- bundle_chromosome_sizes(bundle)
  if (!is.null(sizes)) {
    for (name in names(objects)) {
      gr <- objects[[name]]
      known <- as.character(GenomicRanges::seqnames(gr)) %in% names(sizes)
      out_of_bounds <- known & GenomicRanges::end(gr) >
        sizes[as.character(GenomicRanges::seqnames(gr))]
      if (any(out_of_bounds, na.rm = TRUE)) {
        stop(name, " contains coordinates outside configured chromosome sizes.")
      }
    }
  }
  invisible(TRUE)
}
