safe_promoters <- function(genes, upstream = 2000L) {
  if (length(genes) == 0L) return(GenomicRanges::GRanges())
  value <- suppressWarnings(GenomicRanges::promoters(
    genes, upstream = as.integer(upstream), downstream = 0L, use.names = TRUE
  ))
  GenomicRanges::start(value) <- pmax(1L, GenomicRanges::start(value))
  seq_lengths <- GenomeInfoDb::seqlengths(value)
  bounds <- seq_lengths[as.character(GenomicRanges::seqnames(value))]
  bounded <- !is.na(bounds)
  GenomicRanges::end(value)[bounded] <- pmin(GenomicRanges::end(value)[bounded], bounds[bounded])
  value
}

.txdb_from_annotation <- function(annotation) {
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    return(txdbmaker::makeTxDbFromGRanges(annotation))
  }

  if ("makeTxDbFromGRanges" %in% getNamespaceExports("GenomicFeatures")) {
    return(GenomicFeatures::makeTxDbFromGRanges(annotation))
  }

  stop(
    "Deriving transcript features requires the Bioconductor package 'txdbmaker'.",
    call. = FALSE
  )
}

.first_metadata_value <- function(x) {
  if (inherits(x, "List")) {
    return(vapply(as.list(x), function(value) {
      if (length(value)) as.character(value[[1]]) else NA_character_
    }, character(1)))
  }
  as.character(x)
}

.derive_txdb_feature <- function(txdb, extractor, feature_type) {
  grouped <- extractor(txdb, use.names = FALSE)
  feature_counts <- lengths(grouped)
  if (!sum(feature_counts)) return(GenomicRanges::GRanges())

  transcript_keys <- rep(names(grouped), feature_counts)
  if (length(transcript_keys) != sum(feature_counts)) {
    stop("TxDb feature groups do not contain transcript identifiers.", call. = FALSE)
  }

  transcript_ranges <- GenomicFeatures::transcripts(
    txdb,
    columns = c("tx_id", "tx_name", "gene_id")
  )
  transcript_metadata <- S4Vectors::mcols(transcript_ranges)
  tx_id <- as.character(transcript_metadata$tx_id)
  tx_name <- as.character(transcript_metadata$tx_name)
  gene_id <- .first_metadata_value(transcript_metadata$gene_id)
  matched <- match(transcript_keys, tx_id)

  derived <- unlist(grouped, use.names = FALSE)
  derived$transcript_id <- tx_name[matched]
  missing_tx_name <- is.na(derived$transcript_id) | !nzchar(derived$transcript_id)
  derived$transcript_id[missing_tx_name] <- transcript_keys[missing_tx_name]
  derived$gene_id <- gene_id[matched]
  derived$type <- feature_type
  derived
}

prepare_gene_features <- function(annotation, promoter_upstream = 2000L) {
  if (is.null(annotation)) annotation <- GenomicRanges::GRanges()

  feature_type <- if ("type" %in% names(S4Vectors::mcols(annotation))) {
    tolower(as.character(annotation$type))
  } else {
    rep("", length(annotation))
  }
  pick <- function(values) annotation[feature_type %in% tolower(values)]

  features <- list(
    Genes = pick("gene"),
    CDS = pick("CDS"),
    Introns = pick(c("intron", "intronic")),
    fiveUTRs = pick(c("five_prime_UTR", "five_prime_utr", "5UTR", "5_prime_UTR")),
    threeUTRs = pick(c("three_prime_UTR", "three_prime_utr", "3UTR", "3_prime_UTR"))
  )
  features$Promoters <- safe_promoters(features$Genes, promoter_upstream)
  features$Promoters$type <- "promoter"

  missing_features <- names(features)[
    names(features) %in% c("Introns", "fiveUTRs", "threeUTRs") &
      vapply(features, length, integer(1)) == 0L
  ]

  if (length(missing_features) && length(annotation)) {
    txdb <- tryCatch(
      suppressWarnings(.txdb_from_annotation(annotation)),
      error = function(e) {
        warning(
          "Could not derive missing transcript features: ", conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )

    if (!is.null(txdb)) {
      derivations <- list(
        Introns = list(GenomicFeatures::intronsByTranscript, "intron"),
        fiveUTRs = list(GenomicFeatures::fiveUTRsByTranscript, "five_prime_UTR"),
        threeUTRs = list(GenomicFeatures::threeUTRsByTranscript, "three_prime_UTR")
      )

      for (feature_name in missing_features) {
        specification <- derivations[[feature_name]]
        features[[feature_name]] <- tryCatch(
          .derive_txdb_feature(txdb, specification[[1]], specification[[2]]),
          error = function(e) {
            warning(
              "Could not derive ", feature_name, ": ", conditionMessage(e),
              call. = FALSE
            )
            GenomicRanges::GRanges()
          }
        )
      }
    }
  }

  gene_metadata <- S4Vectors::mcols(features$Genes)
  if (all(c("gene_id", "gene_model_type") %in% names(gene_metadata))) {
    gene_ids <- as.character(gene_metadata$gene_id)
    gene_biotypes <- as.character(gene_metadata$gene_model_type)
    for (feature_name in c("CDS", "Introns", "fiveUTRs", "threeUTRs")) {
      feature <- features[[feature_name]]
      metadata <- S4Vectors::mcols(feature)
      if (length(feature) && "gene_id" %in% names(metadata) &&
          !"gene_model_type" %in% names(metadata)) {
        feature$gene_model_type <- gene_biotypes[match(as.character(feature$gene_id), gene_ids)]
        features[[feature_name]] <- feature
      }
    }
  }

  features[c("Genes", "Promoters", "CDS", "Introns", "fiveUTRs", "threeUTRs")]
}

genome_ann <- function(ann.gr, TE_file.f = GenomicRanges::GRanges(),
                       promoter_upstream = 2000L) {
  if (is.null(ann.gr)) ann.gr <- GenomicRanges::GRanges()
  if (is.null(TE_file.f)) TE_file.f <- GenomicRanges::GRanges()

  feature <- if ("type" %in% names(S4Vectors::mcols(ann.gr))) {
    tolower(as.character(ann.gr$type))
  } else {
    rep("", length(ann.gr))
  }

  pick <- function(values) ann.gr[feature %in% tolower(values)]
  annotation_vec <- prepare_gene_features(ann.gr, promoter_upstream)
  if (length(TE_file.f)) annotation_vec$Transposable_Elements <- TE_file.f

  teg <- pick("transposable_element_gene")
  pseudogene <- pick("pseudogene")
  contains_TEG_n_pseudogene <<- length(teg) > 0L || length(pseudogene) > 0L
  if (length(teg)) annotation_vec$TEG <- teg
  if (length(pseudogene)) annotation_vec$pseudogene <- pseudogene

  annotation_vec[vapply(annotation_vec, length, integer(1)) > 0L]
}
