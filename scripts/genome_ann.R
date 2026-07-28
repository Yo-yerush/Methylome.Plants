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
  genes <- pick("gene")
  cds <- pick("CDS")
  introns <- pick(c("intron", "intronic"))
  five_utrs <- pick(c("five_prime_UTR", "five_prime_utr", "5UTR", "5_prime_UTR"))
  three_utrs <- pick(c("three_prime_UTR", "three_prime_utr", "3UTR", "3_prime_UTR"))

  promoters_gr <- safe_promoters(genes, promoter_upstream)
  promoters_gr$type <- "promoter"

  # Derive absent transcript features when the input carries enough
  # gene/transcript structure for GenomicFeatures.
  needs_derivation <- length(introns) == 0L ||
    length(five_utrs) == 0L ||
    length(three_utrs) == 0L
  if (needs_derivation && length(ann.gr)) {
    txdb <- tryCatch(
      GenomicFeatures::makeTxDbFromGRanges(ann.gr),
      error = function(e) NULL
    )
    if (!is.null(txdb)) {
      if (length(introns) == 0L) {
        introns <- tryCatch(
          unlist(GenomicFeatures::intronsByTranscript(txdb), use.names = FALSE),
          error = function(e) GenomicRanges::GRanges()
        )
        introns$type <- "intron"
      }
      if (length(five_utrs) == 0L) {
        five_utrs <- tryCatch(
          unlist(GenomicFeatures::fiveUTRsByTranscript(txdb), use.names = FALSE),
          error = function(e) GenomicRanges::GRanges()
        )
        five_utrs$type <- "five_prime_UTR"
      }
      if (length(three_utrs) == 0L) {
        three_utrs <- tryCatch(
          unlist(GenomicFeatures::threeUTRsByTranscript(txdb), use.names = FALSE),
          error = function(e) GenomicRanges::GRanges()
        )
        three_utrs$type <- "three_prime_UTR"
      }
    }
  }

  annotation_vec <- list(
    Genes = genes,
    Promoters = promoters_gr,
    CDS = cds,
    Introns = introns,
    fiveUTRs = five_utrs,
    threeUTRs = three_utrs
  )
  if (length(TE_file.f)) annotation_vec$Transposable_Elements <- TE_file.f

  teg <- pick("transposable_element_gene")
  pseudogene <- pick("pseudogene")
  contains_TEG_n_pseudogene <<- length(teg) > 0L || length(pseudogene) > 0L
  if (length(teg)) annotation_vec$TEG <- teg
  if (length(pseudogene)) annotation_vec$pseudogene <- pseudogene

  annotation_vec[vapply(annotation_vec, length, integer(1)) > 0L]
}
