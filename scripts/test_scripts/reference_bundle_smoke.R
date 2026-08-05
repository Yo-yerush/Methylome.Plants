args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(yaml)
})

source(file.path(repo_dir, "scripts", "reference_bundle.R"))
source(file.path(repo_dir, "scripts", "reference_inputs.R"))
source(file.path(repo_dir, "scripts", "trimm_and_rename_seq.R"))
source(file.path(repo_dir, "scripts", "genome_ann.R"))

bundle_path <- file.path(script_dir, "fixtures", "generic_bundle.yaml")
bundle <- read_reference_bundle(bundle_path)
expect_error <- function(expr) {
  failed <- FALSE
  tryCatch(force(expr), error = function(e) failed <<- TRUE)
  stopifnot(failed)
}

stopifnot(reference_bundle_summary(bundle) == "Synthetic plant / test_v1")
stopifnot(identical(
  names(bundle_chromosome_sizes(bundle)),
  c("linkage_A", "linkage_B", "plastome")
))

methylation <- GRanges(
  seqnames = c("scaffold_01", "scaffold_02", "cpDNA"),
  ranges = IRanges(c(100, 200, 50), width = 1),
  context = c("CG", "CHG", "CHH"),
  readsM = 1L,
  readsN = 2L
)
methylation <- harmonize_seqlevels(methylation, bundle)
stopifnot(identical(seqlevels(methylation), c("linkage_A", "linkage_B")))

unconfigured <- GRanges("unplaced_contig", IRanges(1, 10))
stopifnot(identical(seqlevels(harmonize_seqlevels(unconfigured)), "unplaced_contig"))

genes <- harmonize_seqlevels(read_gene_annotation(bundle = bundle), bundle)
tes <- harmonize_seqlevels(read_te_annotation(bundle = bundle), bundle)
descriptions <- read_gene_descriptions(bundle = bundle)
stopifnot(all(c("type", "gene_id", "gene_model_type") %in% names(mcols(genes))))
stopifnot(all(genes$gene_model_type == "protein_coding"))
stopifnot(all(c("gene_id", "Transposon_Family", "Transposon_Super_Family") %in% names(mcols(tes))))
stopifnot(identical(descriptions$gene_id, c("SYN001", "SYN002")))
stopifnot(validate_reference_inputs(methylation, genes, tes, bundle))

gff_like <- GRanges(
  seqnames = rep("linkage_A", 3),
  ranges = IRanges(c(100, 100, 120), c(500, 500, 300)),
  type = c("gene", "mRNA", "CDS"),
  ID = c("GENE_A", "TX_A", "CDS_A"),
  Parent = c(NA, "GENE_A", "TX_A")
)
stopifnot(identical(.resolve_gff_gene_ids(gff_like), rep("GENE_A", 3)))

gff_bundle <- list(annotation = list(
  format = "gff3",
  genes = file.path(script_dir, "fixtures", "genes.gff3"),
  gene_id_field = "gene_id",
  biotype_field = "gene_biotype",
  feature_type_field = "type"
))
gff_genes <- read_gene_annotation(bundle = gff_bundle)
stopifnot(identical(as.character(gff_genes$gene_id), rep("GENE_A", length(gff_genes))))
stopifnot(all(c("ID", "Parent") %in% names(mcols(gff_genes))))

gff_features <- prepare_gene_features(gff_genes, promoter_upstream = 100L)
stopifnot(
  length(gff_features$Introns) == 1L,
  start(gff_features$Introns) == 301L,
  end(gff_features$Introns) == 399L,
  identical(as.character(gff_features$Introns$transcript_id), "TX_A"),
  identical(as.character(gff_features$Introns$gene_id), "GENE_A")
)

gff_descriptions <- read_gene_descriptions(bundle = gff_bundle, annotation = gff_genes)
stopifnot(
  identical(gff_descriptions$gene_id, "GENE_A"),
  is.na(gff_descriptions$Symbol),
  identical(gff_descriptions$Short_description, "Synthetic GFF description")
)

alias_collision <- GRanges(
  seqnames = c("scaffold_01", "linkage_A"),
  ranges = IRanges(c(1, 20), width = 1)
)
expect_error(harmonize_seqlevels(alias_collision, bundle))

out_of_bounds <- methylation
seqlengths(out_of_bounds) <- rep(NA_integer_, length(seqlevels(out_of_bounds)))
end(out_of_bounds)[1] <- bundle_chromosome_sizes(bundle)[["linkage_A"]] + 1L
expect_error(validate_reference_inputs(out_of_bounds, genes, tes, bundle))

annotations <- genome_ann(genes, tes)
stopifnot(all(c("Genes", "Promoters", "CDS", "Introns", "Transposable_Elements") %in% names(annotations)))

caps <- bundle_capabilities(bundle)
stopifnot(caps[["gene_annotation"]], caps[["transposable_elements"]],
          caps[["centromeres"]], caps[["go"]], caps[["functional_groups"]])

message("reference_bundle_smoke: OK")
