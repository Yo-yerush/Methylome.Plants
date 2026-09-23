#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg)) normalizePath(file_arg[[1]], winslash = "/", mustWork = TRUE) else NA_character_
}

.wizard_script <- script_path()
.wizard_repo <- if (!is.na(.wizard_script)) {
  normalizePath(file.path(dirname(.wizard_script), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(.wizard_repo, "scripts", "reference_bundle.R"))
source(file.path(.wizard_repo, "scripts", "reference_inputs.R"))

parse_cli <- function(args) {
  if (!length(args)) return(list(command = "help", options = list()))
  result <- list(command = args[[1]], options = list())
  i <- 2L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    key <- gsub("-", "_", sub("^--", "", token))
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      result$options[[key]] <- TRUE
      i <- i + 1L
    } else {
      result$options[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  result
}

optional_value <- function(x) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(trimws(x)) ||
      tolower(trimws(x)) %in% c("none", "null", "off")) return(NULL)
  trimws(x)
}

existing_path <- function(x, label, required = FALSE) {
  x <- optional_value(x)
  if (is.null(x)) {
    if (required) stop(label, " is required.")
    return(NULL)
  }
  if (!file.exists(x)) stop(label, " does not exist: ", x)
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

package_directories <- function() {
  unique(unlist(lapply(.libPaths(), function(lib) {
    if (!dir.exists(lib)) return(character())
    list.files(lib, full.names = FALSE, recursive = FALSE)
  }), use.names = FALSE))
}

get_orgdb <- function(package) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("AnnotationDbi is required for OrgDb discovery.")
  }
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("OrgDb package is not installed: ", package)
  }
  getExportedValue(package, package)
}

orgdb_chromosome_lengths <- function(package) {
  if (is.null(package)) return(NULL)
  exports <- tryCatch(getNamespaceExports(package), error = function(e) character())
  candidates <- grep("CHRLENGTHS$", exports, value = TRUE)
  if (!length(candidates)) {
    conventional <- paste0(sub("\\.db$", "", package), "CHRLENGTHS")
    candidates <- conventional
  }
  for (candidate in candidates) {
    value <- tryCatch(getExportedValue(package, candidate), error = function(e) NULL)
    if (!is.null(value) && length(value) && !is.null(names(value))) {
      value <- suppressWarnings(as.numeric(value))
      names(value) <- names(tryCatch(
        getExportedValue(package, candidate),
        error = function(e) numeric()
      ))
      valid <- is.finite(value) & value > 0 & !is.na(names(value)) & nzchar(names(value))
      if (any(valid)) return(value[valid])
    }
  }
  NULL
}

orgdb_information <- function(package) {
  db <- get_orgdb(package)
  lengths <- orgdb_chromosome_lengths(package)
  list(
    package = package,
    species = tryCatch(AnnotationDbi::species(db), error = function(e) package),
    taxonomy_id = tryCatch(as.character(AnnotationDbi::taxonomyId(db)), error = function(e) NA_character_),
    keytypes = tryCatch(AnnotationDbi::keytypes(db), error = function(e) character()),
    has_chromosome_lengths = !is.null(lengths)
  )
}

available_orgdb_packages <- function(cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) {
    cached <- tryCatch(utils::read.delim(cache, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(cached) && all(c('package', 'title') %in% names(cached))) return(cached)
  }
  if (!requireNamespace('BiocManager', quietly = TRUE)) return(NULL)
  value <- tryCatch({
    old_timeout <- getOption('timeout')
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = max(30, old_timeout))
    repositories <- BiocManager::repositories()
    matrix <- utils::available.packages(
      contriburl = utils::contrib.url(repositories, type = 'source'),
      filters = list()
    )
    keep <- grepl('^org[.].*[.]db$', matrix[, 'Package'], ignore.case = TRUE)
    data.frame(
      package = matrix[keep, 'Package'],
      title = matrix[keep, 'Title'],
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
  if (!is.null(value) && !is.null(cache)) {
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(value, cache, sep = '\t', quote = FALSE, row.names = FALSE)
  }
  value
}

list_orgdb <- function(query = NULL, include_available = FALSE, cache = NULL) {
  packages <- grep("^org\\..+\\.db$", package_directories(), value = TRUE, ignore.case = TRUE)
  packages <- sort(unique(packages))
  rows <- lapply(packages, function(package) {
    tryCatch(orgdb_information(package), error = function(e) NULL)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!is.null(query) && nzchar(query)) {
    keep <- vapply(rows, function(x) {
      grepl(tolower(query), tolower(paste(x$package, x$species, x$taxonomy_id)), fixed = TRUE)
    }, logical(1))
    rows <- rows[keep]
  }
  for (row in rows) {
    description <- paste0(
      row$species,
      if (isTRUE(row$has_chromosome_lengths)) " [chromosome lengths available]" else "",
      if (!is.na(row$taxonomy_id)) paste0(" [taxid ", row$taxonomy_id, "]") else ""
    )
    cat(row$package, description, 'installed', sep = "\t")
    cat("\n")
  }
  if (isTRUE(include_available)) {
    available <- available_orgdb_packages(cache)
    if (!is.null(available)) {
      available <- available[!available$package %in% packages, , drop = FALSE]
      if (!is.null(query) && nzchar(query)) {
        keep <- grepl(tolower(query), tolower(paste(available$package, available$title)), fixed = TRUE)
        available <- available[keep, , drop = FALSE]
      }
      for (i in seq_len(nrow(available))) {
        cat(
          available$package[[i]],
          paste0(available$title[[i]], ' [installation required]'),
          'available',
          sep = '\t'
        )
        cat('\n')
      }
    }
  }
  invisible(rows)
}

read_kegg_cache <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  value <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("code", "name")
  if (!all(required %in% names(value))) return(NULL)
  value
}

fetch_kegg_organisms <- function(cache = NULL, refresh = FALSE) {
  cached <- read_kegg_cache(cache)
  cache_is_fresh <- !is.null(cached) &&
    as.numeric(difftime(Sys.time(), file.info(cache)$mtime, units = "days")) <= 30
  if (!refresh && cache_is_fresh) return(cached)

  fetched <- tryCatch({
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
      stop("KEGGREST is not installed.")
    }
    organisms <- KEGGREST::keggList("organism")
    data.frame(
      code = as.character(organisms[, 2]),
      name = as.character(organisms[, 3]),
      taxonomy = as.character(organisms[, 4]),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    if (!is.null(cached)) return(cached)
    stop("Could not retrieve the KEGG organism list: ", conditionMessage(e))
  })

  if (!is.null(cache) && !identical(fetched, cached)) {
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(fetched, cache, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  fetched
}

list_kegg <- function(query = NULL, cache = NULL, refresh = FALSE, limit = 200L) {
  organisms <- fetch_kegg_organisms(cache, refresh)
  if (!is.null(query) && nzchar(query)) {
    keep <- grepl(
      tolower(query),
      tolower(paste(organisms$code, organisms$name, organisms$taxonomy)),
      fixed = TRUE
    )
    organisms <- organisms[keep, , drop = FALSE]
  }
  organisms <- utils::head(organisms, as.integer(limit))
  for (i in seq_len(nrow(organisms))) {
    cat(organisms$code[[i]], organisms$name[[i]], sep = "\t")
    cat("\n")
  }
  invisible(organisms)
}

annotation_format <- function(path, configured = NULL) {
  if (!is.null(configured)) return(tolower(configured))
  clean <- sub("\\.gz$", "", tolower(path))
  sub("^.*\\.", "", clean)
}

first_present <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit)) hit[[1]] else NULL
}

inspect_annotation <- function(path, configured_format = NULL) {
  format <- annotation_format(path, configured_format)
  if (format %in% c("gtf", "gff", "gff3")) {
    raw <- rtracklayer::import(path)
  } else if (format == "csv") {
    table <- utils::read.csv(path, check.names = FALSE)
    raw <- GenomicRanges::makeGRangesFromDataFrame(table, keep.extra.columns = TRUE)
  } else {
    stop("The reference wizard supports GTF, GFF3, GFF, and CSV gene annotations.")
  }

  fields <- names(S4Vectors::mcols(raw))
  gene_id_field <- first_present(c("gene_id", "gene", "locus_tag", "ID"), fields)
  transcript_id_field <- first_present(c("transcript_id", "transcript"), fields)
  biotype_field <- first_present(
    c("gene_biotype", "gene_type", "gene_model_type", "biotype"),
    fields
  )
  feature_type_field <- first_present(c("type", "feature", "feature_type"), fields)
  if (is.null(feature_type_field)) stop("Could not identify a feature-type field in the annotation.")

  mini_bundle <- list(annotation = list(
    genes = path,
    format = format,
    gene_id_field = gene_id_field %||% "gene_id",
    transcript_id_field = transcript_id_field,
    biotype_field = biotype_field,
    feature_type_field = feature_type_field,
    protein_coding_values = c("protein_coding", "protein_coding_gene")
  ))
  normalized <- read_gene_annotation(bundle = mini_bundle)
  feature_type <- tolower(as.character(normalized$type))
  gene_rows <- feature_type %in% c("gene", "pseudogene", "transposable_element_gene")
  gene_ids <- unique(as.character(normalized$gene_id[gene_rows]))
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  if (!length(gene_ids)) {
    gene_ids <- unique(as.character(normalized$gene_id))
    gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  }

  biotype_values <- if (!is.null(biotype_field)) {
    unique(as.character(S4Vectors::mcols(raw)[[biotype_field]]))
  } else {
    character()
  }
  coding_values <- biotype_values[grepl("protein.*coding|coding.*protein", biotype_values, ignore.case = TRUE)]
  if (!length(coding_values)) coding_values <- "protein_coding"

  seqnames <- unique(as.character(GenomicRanges::seqnames(normalized)))
  list(
    format = format,
    raw = raw,
    normalized = normalized,
    gene_ids = gene_ids,
    seqnames = seqnames,
    gene_id_field = gene_id_field %||% "gene_id",
    transcript_id_field = transcript_id_field,
    biotype_field = biotype_field,
    feature_type_field = feature_type_field,
    protein_coding_values = coding_values
  )
}

read_size_table <- function(path) {
  clean <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (all(c("seqname", "length") %in% names(clean))) {
    result <- clean[, c("seqname", "length")]
  } else if (ncol(clean) >= 2L) {
    result <- clean[, 1:2]
    names(result) <- c("seqname", "length")
  } else {
    stop("Chromosome sizes must contain seqname and length columns.")
  }
  result$seqname <- as.character(result$seqname)
  result$length <- suppressWarnings(as.numeric(result$length))
  if (any(!is.finite(result$length) | result$length <= 0) || anyDuplicated(result$seqname)) {
    stop("Chromosome sizes contain invalid lengths or duplicate names.")
  }
  result
}

fasta_sizes <- function(path) {
  if (grepl("\\.fai$", path, ignore.case = TRUE)) {
    fai <- utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(fai) < 2L) stop("FASTA index must contain at least two columns.")
    return(data.frame(seqname = as.character(fai[[1]]), length = as.numeric(fai[[2]])))
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings is required to obtain chromosome lengths from FASTA.")
  }
  index <- Biostrings::fasta.index(path, seqtype = "DNA")
  data.frame(
    seqname = sub("\\s.*$", "", as.character(index$desc)),
    length = as.numeric(index$seqlength),
    stringsAsFactors = FALSE
  )
}

normal_seq_token <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("^chromosome[_ .-]*", "", x)
  x <- sub("^chr", "", x)
  numeric_name <- grepl("^[0-9]+$", x)
  x[numeric_name] <- as.character(as.integer(x[numeric_name]))
  x
}

map_annotation_to_reference <- function(annotation_names, reference_names) {
  mapped <- rep(NA_character_, length(annotation_names))
  exact <- match(annotation_names, reference_names)
  mapped[!is.na(exact)] <- reference_names[exact[!is.na(exact)]]

  unresolved <- which(is.na(mapped))
  reference_tokens <- normal_seq_token(reference_names)
  for (i in unresolved) {
    hits <- which(reference_tokens == normal_seq_token(annotation_names[[i]]))
    if (length(hits) == 1L) mapped[[i]] <- reference_names[[hits]]
  }
  stats::setNames(mapped, annotation_names)
}

merge_aliases <- function(user_alias_path, inferred_aliases, output_dir, file_prefix = '') {
  aliases <- data.frame(alias = character(), canonical = character())
  if (!is.null(user_alias_path)) {
    aliases <- utils::read.delim(user_alias_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("alias", "canonical") %in% names(aliases))) {
      stop("Sequence alias file must contain alias and canonical columns.")
    }
    aliases <- aliases[, c("alias", "canonical")]
  }
  if (length(inferred_aliases)) {
    inferred <- data.frame(
      alias = names(inferred_aliases),
      canonical = unname(inferred_aliases),
      stringsAsFactors = FALSE
    )
    aliases <- rbind(aliases, inferred)
  }
  aliases <- unique(aliases[nzchar(aliases$alias) & nzchar(aliases$canonical), , drop = FALSE])
  if (!nrow(aliases)) return(NULL)
  if (anyDuplicated(aliases$alias)) stop("Sequence alias mappings contain duplicate aliases.")
  path <- file.path(output_dir, paste0(file_prefix, "seqname_aliases.tsv"))
  utils::write.table(aliases, path, sep = "\t", quote = FALSE, row.names = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

select_orgdb_keytype <- function(package, gene_ids, requested = "auto") {
  db <- get_orgdb(package)
  available <- AnnotationDbi::keytypes(db)
  if (!identical(tolower(requested), "auto")) {
    if (!requested %in% available) {
      stop("OrgDb key type is unavailable: ", requested)
    }
    keys <- AnnotationDbi::keys(db, keytype = requested)
    overlap <- sum(gene_ids %in% keys)
    return(list(keytype = requested, overlap = overlap, total = length(gene_ids)))
  }

  scores <- vapply(available, function(keytype) {
    keys <- tryCatch(AnnotationDbi::keys(db, keytype = keytype), error = function(e) character())
    sum(gene_ids %in% keys)
  }, numeric(1))
  best <- if (length(scores)) names(which.max(scores)) else NA_character_
  list(
    keytype = if (!is.na(best)) best else NULL,
    overlap = if (!is.na(best)) unname(max(scores)) else 0,
    total = length(gene_ids)
  )
}

inspect_te <- function(path) {
  if (is.null(path)) return(list())
  format <- annotation_format(path)
  if (format %in% c("gff", "gff3", "gtf", "bed")) {
    raw <- rtracklayer::import(path)
    fields <- names(S4Vectors::mcols(raw))
  } else {
    separator <- if (format == "csv") "," else "\t"
    raw <- utils::read.table(
      path, header = TRUE, sep = separator, quote = "\"", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, nrows = 100
    )
    fields <- names(raw)
  }
  if ("Transposon_Name" %in% fields) {
    return(list(format = "tair_legacy"))
  }
  list(
    format = format,
    id_field = first_present(c("gene_id", "ID", "Name", "name", "repeat_id", "te_id"), fields) %||% "gene_id",
    family_field = first_present(c("Transposon_Family", "family", "te_family", "Family"), fields),
    superfamily_field = first_present(
      c("Transposon_Super_Family", "superfamily", "te_superfamily", "Superfamily"),
      fields
    )
  )
}

description_id_field <- function(path) {
  if (is.null(path)) return(NULL)
  separator <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
  value <- utils::read.table(
    path, header = TRUE, sep = separator, quote = "\"", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, nrows = 5
  )
  first_present(c("gene_id", "gene", "ID", "locus_tag"), names(value)) %||% names(value)[[1]]
}

validate_mapping_table <- function(path, required_columns, label, gene_ids = NULL) {
  if (is.null(path)) return(NULL)
  separator <- if (grepl('[.]csv([.]gz)?$', path, ignore.case = TRUE)) ',' else '\t'
  value <- utils::read.table(
    path, header = TRUE, sep = separator, quote = '', comment.char = '',
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!all(required_columns %in% names(value))) {
    stop(label, ' must contain columns: ', paste(required_columns, collapse = ', '))
  }
  overlap <- if (!is.null(gene_ids) && 'gene_id' %in% required_columns) {
    sum(unique(as.character(value$gene_id)) %in% gene_ids)
  } else {
    NA_integer_
  }
  list(rows = nrow(value), overlap = overlap)
}

validate_bed_file <- function(path, label) {
  if (is.null(path)) return(invisible(TRUE))
  value <- utils::read.table(
    path, header = FALSE, sep = '\t', quote = '', comment.char = '',
    stringsAsFactors = FALSE, nrows = 100
  )
  if (ncol(value) < 3L) stop(label, ' must contain at least three BED columns.')
  starts <- suppressWarnings(as.numeric(value[[2]]))
  ends <- suppressWarnings(as.numeric(value[[3]]))
  if (any(!is.finite(starts) | !is.finite(ends) | starts < 0 | ends <= starts)) {
    stop(label, ' contains invalid BED coordinates.')
  }
  invisible(TRUE)
}

absolute_optional <- function(value, label) existing_path(value, label, required = FALSE)

generate_bundle <- function(options) {
  organism <- optional_value(options$organism)
  assembly <- optional_value(options$assembly)
  annotation_path <- existing_path(options$annotation, "Gene annotation", required = TRUE)
  output_path <- optional_value(options$output)
  if (is.null(organism)) stop("--organism is required.")
  if (is.null(assembly)) stop("--assembly is required.")
  if (is.null(output_path)) stop("--output is required.")

  output_path <- normalizePath(output_path, winslash = "/", mustWork = FALSE)
  output_dir <- dirname(output_path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  bundle_stem <- tools::file_path_sans_ext(basename(output_path))
  auxiliary_prefix <- paste0(bundle_stem, '.')

  annotation <- inspect_annotation(annotation_path, optional_value(options$annotation_format))
  requested_primary <- optional_value(options$primary_seqlevels)
  annotation_primary <- if (is.null(requested_primary)) {
    annotation$seqnames
  } else {
    trimws(strsplit(requested_primary, ",", fixed = TRUE)[[1]])
  }
  missing_primary <- setdiff(annotation_primary, annotation$seqnames)
  if (length(missing_primary)) {
    stop("Configured primary sequence names are absent from the annotation: ",
         paste(missing_primary, collapse = ", "))
  }

  orgdb_package <- optional_value(options$orgdb_package)
  orgdb_keytype <- NULL
  orgdb_match <- NULL
  warnings <- character()
  if (!is.null(orgdb_package)) {
    match <- select_orgdb_keytype(
      orgdb_package,
      annotation$gene_ids,
      optional_value(options$orgdb_keytype) %||% "auto"
    )
    orgdb_keytype <- match$keytype
    orgdb_match <- match
    if (is.null(orgdb_keytype) || match$overlap == 0L) {
      warnings <- c(
        warnings,
        paste0("No GTF gene identifiers matched the selected OrgDb package ", orgdb_package, ".")
      )
    }
  }

  chromosome_sizes_path <- existing_path(
    options$chromosome_sizes, "Chromosome-size table", required = FALSE
  )
  fasta_path <- existing_path(options$fasta, "Reference FASTA/FAI", required = FALSE)
  user_alias_path <- existing_path(options$seqname_aliases, "Sequence alias table", required = FALSE)

  size_table <- NULL
  length_source <- "not supplied"
  strict_size_source <- FALSE
  if (!is.null(chromosome_sizes_path)) {
    size_table <- read_size_table(chromosome_sizes_path)
    length_source <- "chromosome-size table"
    strict_size_source <- TRUE
  } else if (!is.null(fasta_path)) {
    size_table <- fasta_sizes(fasta_path)
    length_source <- if (grepl("\\.fai$", fasta_path, ignore.case = TRUE)) "FASTA index" else "FASTA"
    strict_size_source <- TRUE
  } else if (!is.null(orgdb_package) && !isTRUE(options$disable_orgdb_lengths)) {
    org_lengths <- orgdb_chromosome_lengths(orgdb_package)
    if (!is.null(org_lengths)) {
      size_table <- data.frame(
        seqname = names(org_lengths),
        length = as.numeric(org_lengths),
        stringsAsFactors = FALSE
      )
      length_source <- paste0(orgdb_package, " CHRLENGTHS")
    }
  }

  canonical_primary <- annotation_primary
  sequence_mapping <- stats::setNames(annotation_primary, annotation_primary)
  inferred_aliases <- character()
  generated_sizes_path <- NULL
  if (!is.null(size_table)) {
    mapping <- map_annotation_to_reference(annotation_primary, size_table$seqname)
    unmatched <- names(mapping)[is.na(mapping)]
    if (length(unmatched)) {
      message <- paste0(
        "Chromosome lengths could not be matched to annotation sequences: ",
        paste(unmatched, collapse = ", ")
      )
      if (strict_size_source) stop(message)
      warnings <- c(warnings, message, "OrgDb chromosome lengths were not used.")
      size_table <- NULL
      length_source <- "not supplied"
    } else {
      sequence_mapping <- mapping
      canonical_primary <- unname(mapping)
      changed <- names(mapping) != unname(mapping)
      inferred_aliases <- unname(mapping[changed])
      names(inferred_aliases) <- names(mapping)[changed]
      selected_sizes <- size_table[match(canonical_primary, size_table$seqname), , drop = FALSE]
      names(selected_sizes) <- c("seqname", "length")
      annotation_seqnames <- as.character(GenomicRanges::seqnames(annotation$normalized))
      annotation_canonical <- unname(mapping[annotation_seqnames])
      size_lookup <- stats::setNames(selected_sizes$length, selected_sizes$seqname)
      out_of_bounds <- !is.na(annotation_canonical) &
        GenomicRanges::end(annotation$normalized) > size_lookup[annotation_canonical]
      if (any(out_of_bounds, na.rm = TRUE)) {
        stop(
          "The annotation contains coordinates beyond the selected chromosome lengths. ",
          "Check that the GTF/GFF and genome-length source use the same assembly."
        )
      }
      generated_sizes_path <- file.path(
        output_dir, paste0(auxiliary_prefix, "chromosome_sizes.tsv")
      )
      utils::write.table(
        selected_sizes, generated_sizes_path, sep = "\t",
        quote = FALSE, row.names = FALSE
      )
      generated_sizes_path <- normalizePath(
        generated_sizes_path, winslash = "/", mustWork = TRUE
      )
    }
  } else {
    warnings <- c(
      warnings,
      "Exact chromosome lengths are unavailable; length-dependent exports will use observed ranges."
    )
  }

  aliases_path <- merge_aliases(
    user_alias_path, inferred_aliases, output_dir, auxiliary_prefix
  )
  descriptions_path <- absolute_optional(options$descriptions, "Gene-description table")
  te_path <- absolute_optional(options$te, "TE annotation")
  te_info <- inspect_te(te_path)
  centromeres_path <- absolute_optional(options$centromeres, "Centromere BED")
  heterochromatin_path <- absolute_optional(options$heterochromatin, "Heterochromatin BED")
  stable_gbm_path <- absolute_optional(options$stable_gbm, "Stable gbM list")
  dynamic_gbm_path <- absolute_optional(options$dynamic_gbm, "Dynamic gbM list")
  tfbs_path <- absolute_optional(options$tfbs, "TFBS BED")
  gene_to_go_path <- absolute_optional(options$gene_to_go, "Gene-to-GO table")
  gene_sets_path <- absolute_optional(options$gene_sets, "Functional gene-set file")
  kegg_id_map_path <- absolute_optional(options$kegg_id_map, "KEGG ID mapping table")

  gene_to_go_info <- validate_mapping_table(
    gene_to_go_path, c('gene_id', 'go_id'), 'Gene-to-GO table', annotation$gene_ids
  )
  kegg_map_info <- validate_mapping_table(
    kegg_id_map_path, c('gene_id', 'kegg_id'), 'KEGG ID mapping table', annotation$gene_ids
  )
  if (!is.null(gene_to_go_info) && gene_to_go_info$overlap == 0L) {
    warnings <- c(warnings, 'No GTF gene identifiers matched the custom gene-to-GO table.')
  }
  if (!is.null(kegg_map_info) && kegg_map_info$overlap == 0L) {
    warnings <- c(warnings, 'No GTF gene identifiers matched the custom KEGG ID table.')
  }
  if (!is.null(optional_value(options$kegg_organism)) && is.null(kegg_id_map_path)) {
    warnings <- c(
      warnings,
      'KEGG was configured without an ID map; verify that GTF gene IDs are the identifiers returned by KEGG.'
    )
  }
  validate_bed_file(centromeres_path, 'Centromere BED')
  validate_bed_file(heterochromatin_path, 'Heterochromatin BED')
  validate_bed_file(tfbs_path, 'TFBS BED')

  validate_gene_list <- function(path, label) {
    if (is.null(path)) return(invisible(NULL))
    ids <- trimws(readLines(path, warn = FALSE))
    ids <- unique(ids[nzchar(ids)])
    if (!length(intersect(ids, annotation$gene_ids))) {
      warnings <<- c(warnings, paste0('No GTF gene identifiers matched the ', label, '.'))
    }
    invisible(ids)
  }
  validate_gene_list(stable_gbm_path, 'stable gbM list')
  validate_gene_list(dynamic_gbm_path, 'dynamic gbM list')

  chloroplast <- optional_value(options$chloroplast_seqlevels)
  chloroplast <- if (is.null(chloroplast)) character() else trimws(strsplit(chloroplast, ",", fixed = TRUE)[[1]])
  mitochondria <- optional_value(options$mitochondrial_seqlevels)
  mitochondria <- if (is.null(mitochondria)) character() else trimws(strsplit(mitochondria, ",", fixed = TRUE)[[1]])
  canonicalize_declared <- function(values) {
    mapped <- unname(sequence_mapping[values])
    mapped[is.na(mapped)] <- values[is.na(mapped)]
    unique(mapped)
  }
  chloroplast <- canonicalize_declared(chloroplast)
  mitochondria <- canonicalize_declared(mitochondria)

  species_id <- optional_value(options$species_id)
  if (is.null(species_id)) {
    species_id <- tolower(gsub("[^A-Za-z0-9]+", "_", organism))
    species_id <- gsub("^_|_$", "", species_id)
  }

  bundle <- list(
    schema_version = 1L,
    species = list(
      id = species_id,
      display_name = organism,
      assembly = assembly
    ),
    genome = list(
      fasta = if (!is.null(fasta_path) && !grepl("\\.fai$", fasta_path, ignore.case = TRUE)) fasta_path else NULL,
      chromosome_sizes = generated_sizes_path,
      seqname_aliases = aliases_path,
      primary_seqlevels = unique(canonical_primary),
      chloroplast_seqlevels = chloroplast,
      mitochondrial_seqlevels = mitochondria,
      include_non_primary = FALSE
    ),
    annotation = list(
      genes = annotation_path,
      format = annotation$format,
      gene_id_field = annotation$gene_id_field,
      transcript_id_field = annotation$transcript_id_field,
      biotype_field = annotation$biotype_field,
      feature_type_field = annotation$feature_type_field,
      protein_coding_values = as.list(annotation$protein_coding_values),
      promoter_upstream = as.integer(options$promoter_upstream %||% 2000L),
      descriptions = descriptions_path,
      description_gene_id_field = description_id_field(descriptions_path),
      transposable_elements = te_path,
      te_format = te_info$format %||% NULL,
      te_id_field = te_info$id_field %||% NULL,
      te_family_field = te_info$family_field %||% NULL,
      te_superfamily_field = te_info$superfamily_field %||% NULL
    ),
    regions = list(
      centromeres = centromeres_path,
      heterochromatin = heterochromatin_path
    ),
    functional = list(
      orgdb_package = orgdb_package,
      orgdb_keytype = orgdb_keytype,
      gene_to_go = gene_to_go_path,
      kegg_organism = optional_value(options$kegg_organism),
      kegg_id_map = kegg_id_map_path,
      tfbs = tfbs_path,
      gene_sets = gene_sets_path,
      stable_gbm = stable_gbm_path,
      dynamic_gbm = dynamic_gbm_path
    ),
    generated = list(
      created_by = "Methylome.Plants reference wizard",
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      chromosome_length_source = length_source,
      warnings = as.list(unique(warnings))
    )
  )

  if (!requireNamespace("yaml", quietly = TRUE)) stop("The yaml package is required.")
  yaml::write_yaml(bundle, output_path)
  resolved <- read_reference_bundle(output_path)

  cat("Bundle", normalizePath(output_path, winslash = "/", mustWork = TRUE), sep = "\t")
  cat("\n")
  cat("Organism", reference_bundle_summary(resolved), sep = "\t")
  cat("\n")
  cat("Annotation", paste(length(annotation$normalized), "features;", length(annotation$gene_ids), "gene IDs"), sep = "\t")
  cat("\n")
  cat("Sequences", paste(length(canonical_primary), "primary sequence levels"), sep = "\t")
  cat("\n")
  cat("Chromosome lengths", length_source, sep = "\t")
  cat("\n")
  if (!is.null(orgdb_package)) {
    match_percent <- if (!is.null(orgdb_match) && orgdb_match$total > 0L) {
      round(100 * orgdb_match$overlap / orgdb_match$total, 1)
    } else {
      0
    }
    cat(
      "GO / OrgDb",
      paste0(orgdb_package, " [", orgdb_keytype %||% "no key type", "; ", match_percent, "% gene-ID match]"),
      sep = "\t"
    )
    cat("\n")
  } else if (!is.null(gene_to_go_path)) {
    cat("GO", "custom gene-to-GO table", sep = "\t")
    cat("\n")
  } else {
    cat("GO", "not configured", sep = "\t")
    cat("\n")
  }
  cat("KEGG", optional_value(options$kegg_organism) %||% "not configured", sep = "\t")
  cat("\n")
  configured_resources <- c(
    if (!is.null(te_path)) "TE",
    if (!is.null(stable_gbm_path) || !is.null(dynamic_gbm_path)) "gbM",
    if (!is.null(centromeres_path)) "centromeres",
    if (!is.null(heterochromatin_path)) "heterochromatin",
    if (!is.null(tfbs_path)) "TFBS",
    if (!is.null(gene_sets_path)) "gene sets",
    if (!is.null(descriptions_path)) "descriptions"
  )
  cat(
    "Optional resources",
    if (length(configured_resources)) paste(configured_resources, collapse = ", ") else "none",
    sep = "\t"
  )
  cat("\n")
  available_analyses <- c(
    "core methylation", "DMRs", "gene annotation",
    if (!is.null(te_path)) "TE analysis",
    if (!is.null(te_path) && !is.null(centromeres_path)) "TE-centromere distance",
    if (!is.null(stable_gbm_path) || !is.null(dynamic_gbm_path)) "gbM QC",
    if (!is.null(tfbs_path)) "TFBS",
    if (!is.null(gene_sets_path)) "functional groups",
    if (!is.null(orgdb_package) || !is.null(gene_to_go_path)) "GO",
    if (!is.null(optional_value(options$kegg_organism))) "KEGG"
  )
  cat("Analysis support", paste(available_analyses, collapse = ", "), sep = "\t")
  cat("\n")
  for (warning in unique(warnings)) {
    cat("Warning", warning, sep = "\t")
    cat("\n")
  }
  invisible(output_path)
}

print_help <- function() {
  cat(
    "Methylome.Plants reference wizard backend\n\n",
    "Commands:\n",
    "  list-orgdb [--query text] [--include-available] [--cache path]\n",
    "  list-kegg [--query text] [--cache path] [--refresh]\n",
    "  generate --organism name --assembly name --annotation file --output bundle.yaml [options]\n",
    sep = ""
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- parse_cli(args)
  command <- parsed$command
  options <- parsed$options
  if (command == "list-orgdb") {
    list_orgdb(
      optional_value(options$query),
      isTRUE(options$include_available),
      optional_value(options$cache)
    )
  } else if (command == "list-kegg") {
    list_kegg(
      optional_value(options$query),
      optional_value(options$cache),
      isTRUE(options$refresh),
      as.integer(options$limit %||% 200L)
    )
  } else if (command == "generate") {
    generate_bundle(options)
  } else {
    print_help()
  }
}

if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(e) {
      message("Reference wizard error: ", conditionMessage(e))
      quit(status = 1L)
    }
  )
}
