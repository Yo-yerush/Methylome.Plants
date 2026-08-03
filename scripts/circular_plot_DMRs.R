DMRs_circular_plot <- function(ann.file, TE_4_dens, comparison_name,
                               up_col = "#FF0000", down_col = "#304ed1",
                               include_genes = TRUE, reference_bundle = NULL) {
  cntx_file <- function(context) {
    dmrs_file <- read.csv(paste0("DMRs_", context, "_", comparison_name, ".csv"))
    dmrs_file[, c("seqnames", "start", "end", "log2FC")]
  }

  dmrs <- list(CG = cntx_file("CG"), CHG = cntx_file("CHG"), CHH = cntx_file("CHH"))
  configured_sizes <- bundle_chromosome_sizes(reference_bundle)
  configured_primary <- as.character(bundle_get(reference_bundle, "genome.primary_seqlevels", character()))

  if (!is.null(configured_sizes)) {
    chromosomes <- if (length(configured_primary)) {
      intersect(configured_primary, names(configured_sizes))
    } else {
      names(configured_sizes)
    }
    chromosome_sizes <- configured_sizes[chromosomes]
  } else {
    chromosomes <- unique(c(
      as.character(seqlevels(ann.file)),
      unlist(lapply(dmrs, function(x) as.character(x$seqnames)))
    ))
    annotation_lengths <- seqlengths(ann.file)[chromosomes]
    observed_lengths <- vapply(chromosomes, function(chr) {
      ends <- c(
        end(ann.file[as.character(seqnames(ann.file)) == chr]),
        unlist(lapply(dmrs, function(x) x$end[x$seqnames == chr]))
      )
      if (length(ends)) max(ends, na.rm = TRUE) else NA_real_
    }, numeric(1))
    chromosome_sizes <- annotation_lengths
    chromosome_sizes[!is.finite(chromosome_sizes)] <- observed_lengths[!is.finite(chromosome_sizes)]
  }

  valid_chromosomes <- is.finite(chromosome_sizes) & chromosome_sizes > 0
  chromosomes <- chromosomes[valid_chromosomes]
  chromosome_sizes <- chromosome_sizes[valid_chromosomes]
  if (!length(chromosomes)) stop("No chromosome lengths are available for the DMR density plot.")

  # Scale from 3 inches for 5 chromosomes to 4.5 inches for 12,
  # while keeping unusually small or large genomes within practical limits.
  plot_size <- 3 + (length(chromosomes) - 5L) * (1.5 / 7)
  plot_size <- max(3, min(6, plot_size))

  chromosome_frame <- data.frame(
    seqnames = chromosomes,
    start = 0,
    end = as.numeric(chromosome_sizes),
    stringsAsFactors = FALSE
  )
  dmrs <- lapply(dmrs, function(x) {
    x <- x[x$seqnames %in% chromosomes, , drop = FALSE]
    limits <- chromosome_sizes[match(x$seqnames, chromosomes)]
    x[x$start >= 0 & x$end <= limits, , drop = FALSE]
  })

  genes_type <- ann.file[ann.file$type == "gene"]
  annotation_tracks <- list()
  annotation_colours <- character()
  annotation_labels <- character()
  if (include_genes && length(genes_type)) {
    annotation_tracks$Genes <- as.data.frame(genes_type)[, 1:3]
    annotation_colours <- c(annotation_colours, "#9c9c9c")
    annotation_labels <- c(annotation_labels, "Genes")
  }
  if (length(TE_4_dens)) {
    annotation_tracks$TEs <- as.data.frame(TE_4_dens)[, 1:3]
    annotation_colours <- c(annotation_colours, "#fcba0360")
    annotation_labels <- c(annotation_labels, "TEs")
  }

  output_name <- paste0("DMRs_Density_", comparison_name)
  output_type <- as.character(formals(img_device)$img_type)
  output_file <- paste0(output_name, ".", output_type)
  original_device <- dev.cur()
  completed <- FALSE
  img_device(output_name, w = plot_size, h = plot_size)
  on.exit({
    try(circos.clear(), silent = TRUE)
    if (dev.cur() != original_device) try(dev.off(), silent = TRUE)
    if (!completed && file.exists(output_file)) unlink(output_file)
  }, add = TRUE)

  par(mar = c(0, 0, 0, 0))
  circos.clear()
  circos.par(
    gap.degree = c(rep(1, length(chromosomes) - 1L), 28),
    start.degree = 90,
    points.overflow.warning = FALSE
  )
  circos.genomicInitialize(
    chromosome_frame,
    sector.names = chromosomes,
    axis.labels.cex = 0.45,
    labels.cex = 0.8
  )

  for (context in names(dmrs)) {
    context_dmrs <- dmrs[[context]]
    direction_tracks <- list(
      gain = context_dmrs[context_dmrs$log2FC > 0, 1:3, drop = FALSE],
      loss = context_dmrs[context_dmrs$log2FC < 0, 1:3, drop = FALSE]
    )
    direction_colours <- c(gain = paste0(up_col, 80), loss = paste0(down_col, 80))
    populated <- vapply(direction_tracks, nrow, integer(1)) > 0L
    direction_tracks <- direction_tracks[populated]
    direction_colours <- direction_colours[populated]
    if (!length(direction_tracks)) next

    circos.genomicDensity(
      unname(direction_tracks),
      bg.col = "#fafcff",
      bg.border = NA,
      count_by = "number",
      col = unname(direction_colours),
      border = TRUE,
      track.height = 0.145,
      track.margin = c(0, 0)
    )
    circos.text(
      chromosomes[1], x = 0, y = 1, labels = context,
      facing = "downward", cex = 0.65, adj = c(1, -0.4)
    )
  }

  if (length(annotation_tracks)) {
    circos.genomicDensity(
      unname(annotation_tracks),
      bg.col = "#fafcff",
      bg.border = NA,
      count_by = "number",
      col = annotation_colours,
      border = TRUE,
      track.height = 0.145,
      track.margin = c(0, 0)
    )
  }

  circos.clear()
  legend(
    "topleft",
    legend = c("Hyper-DMRs", "Hypo-DMRs", "Overlay", annotation_labels),
    fill = c(paste0(up_col, 95), paste0(down_col, 95), "#8208b695", annotation_colours),
    cex = 0.55,
    bty = "n"
  )

  dev.off()
  completed <- TRUE
  invisible(output_file)
}
