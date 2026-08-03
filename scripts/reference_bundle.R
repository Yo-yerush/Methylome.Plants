`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

bundle_get <- function(bundle, path, default = NULL) {
  if (is.null(bundle)) return(default)
  keys <- strsplit(path, "\\.", fixed = FALSE)[[1]]
  value <- bundle
  for (key in keys) {
    if (!is.list(value) || is.null(value[[key]])) return(default)
    value <- value[[key]]
  }
  value
}

.bundle_resolve_path <- function(value, bundle_dir) {
  if (is.null(value) || length(value) == 0 || !nzchar(value)) return(NULL)
  if (grepl("^(https?|ftp)://", value, ignore.case = TRUE)) return(value)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", value)) return(normalizePath(value, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(bundle_dir, value), winslash = "/", mustWork = FALSE)
}

read_reference_bundle <- function(path) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) stop("Reference bundle does not exist: ", path)
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' R package is required to read a reference bundle.")
  }

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  document <- yaml::read_yaml(path)
  is_resolved_document <- is.list(document$reference_bundle)
  bundle <- if (is_resolved_document) document$reference_bundle else document
  if (!identical(as.integer(bundle$schema_version), 1L)) {
    stop("Unsupported reference-bundle schema_version in ", path)
  }
  if (is.null(bundle$species$id) || is.null(bundle$species$assembly)) {
    stop("Reference bundle must define species.id and species.assembly.")
  }

  source_bundle_dir <- if (is_resolved_document &&
      !is.null(bundle$bundle_dir) && dir.exists(bundle$bundle_dir)) {
    normalizePath(bundle$bundle_dir, winslash = "/", mustWork = TRUE)
  } else {
    dirname(path)
  }
  bundle$bundle_path <- path
  bundle$bundle_dir <- source_bundle_dir
  path_fields <- list(
    c("genome", "fasta"),
    c("genome", "chromosome_sizes"),
    c("genome", "seqname_aliases"),
    c("annotation", "genes"),
    c("annotation", "descriptions"),
    c("annotation", "transposable_elements"),
    c("regions", "centromeres"),
    c("regions", "heterochromatin"),
    c("functional", "tfbs"),
    c("functional", "stable_gbm"),
    c("functional", "dynamic_gbm"),
    c("functional", "gene_sets"),
    c("functional", "gene_to_go"),
    c("functional", "kegg_id_map"),
    c("functional", "gene_sets_adapter")
  )
  for (keys in path_fields) {
    if (!is.null(bundle[[keys[1]]][[keys[2]]])) {
      resolved <- .bundle_resolve_path(
        bundle[[keys[1]]][[keys[2]]],
        bundle$bundle_dir
      )
      if (!grepl('^(https?|ftp)://', resolved, ignore.case = TRUE) && !file.exists(resolved)) {
        stop('Configured reference resource does not exist (', paste(keys, collapse = '.'), '): ', resolved)
      }
      bundle[[keys[1]]][[keys[2]]] <- resolved
    }
  }
  class(bundle) <- c("methylome_reference_bundle", class(bundle))
  sizes <- bundle_chromosome_sizes(bundle)
  primary <- as.character(bundle_get(bundle, 'genome.primary_seqlevels', character()))
  if (length(primary) && !is.null(sizes)) {
    missing_sizes <- setdiff(primary, names(sizes))
    if (length(missing_sizes)) {
      stop('No chromosome size is defined for primary sequence levels: ',
           paste(missing_sizes, collapse = ', '))
    }
  }
  bundle
}

read_bundle_table <- function(path, header = TRUE) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  data.table::fread(path, data.table = FALSE, header = header, showProgress = FALSE)
}

bundle_chromosome_sizes <- function(bundle) {
  path <- bundle_get(bundle, "genome.chromosome_sizes")
  sizes <- read_bundle_table(path)
  if (is.null(sizes)) return(NULL)
  if (!all(c("seqname", "length") %in% names(sizes))) {
    stop("chromosome_sizes must contain 'seqname' and 'length' columns: ", path)
  }
  values <- as.numeric(sizes$length)
  if (any(!is.finite(values) | values <= 0) || anyDuplicated(sizes$seqname)) {
    stop("chromosome_sizes contains duplicate sequence names or invalid lengths: ", path)
  }
  names(values) <- as.character(sizes$seqname)
  values
}

bundle_alias_map <- function(bundle) {
  path <- bundle_get(bundle, "genome.seqname_aliases")
  aliases <- read_bundle_table(path)
  if (is.null(aliases)) return(character())
  if (!all(c("alias", "canonical") %in% names(aliases))) {
    stop("seqname_aliases must contain 'alias' and 'canonical' columns: ", path)
  }
  if (anyDuplicated(aliases$alias)) {
    stop("seqname_aliases contains duplicate aliases: ", path)
  }
  stats::setNames(as.character(aliases$canonical), as.character(aliases$alias))
}

read_bundle_regions <- function(bundle, key) {
  path <- bundle_get(bundle, paste0("regions.", key))
  if (is.null(path) || !file.exists(path)) return(GenomicRanges::GRanges())
  regions <- data.table::fread(
    path,
    header = FALSE,
    data.table = FALSE,
    showProgress = FALSE,
    col.names = c("seqnames", "start0", "end")
  )
  GenomicRanges::makeGRangesFromDataFrame(
    transform(regions, start = start0 + 1L),
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = FALSE
  )
}

bundle_capabilities <- function(bundle) {
  path_exists <- function(path) {
    !is.null(path) && (grepl("^(https?|ftp)://", path) || file.exists(path))
  }
  c(
    gene_annotation = path_exists(bundle_get(bundle, "annotation.genes")),
    descriptions = path_exists(bundle_get(bundle, "annotation.descriptions")),
    transposable_elements = path_exists(bundle_get(bundle, "annotation.transposable_elements")),
    centromeres = path_exists(bundle_get(bundle, "regions.centromeres")),
    heterochromatin = path_exists(bundle_get(bundle, "regions.heterochromatin")),
    chloroplast = length(bundle_get(bundle, "genome.chloroplast_seqlevels", character())) > 0,
    tfbs = path_exists(bundle_get(bundle, "functional.tfbs")),
    go = !is.null(bundle_get(bundle, "functional.orgdb_package")) ||
      path_exists(bundle_get(bundle, "functional.gene_to_go")),
    kegg = !is.null(bundle_get(bundle, "functional.kegg_organism")),
    functional_groups = path_exists(bundle_get(bundle, "functional.gene_sets")) ||
      path_exists(bundle_get(bundle, "functional.gene_sets_adapter"))
  )
}

reference_bundle_summary <- function(bundle) {
  if (is.null(bundle)) return("custom inputs (no reference bundle)")
  paste0(
    bundle_get(bundle, "species.display_name", bundle_get(bundle, "species.id")),
    " / ",
    bundle_get(bundle, "species.assembly")
  )
}
