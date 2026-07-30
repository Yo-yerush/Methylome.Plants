args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else getwd()
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(yaml)
})

source(file.path(repo_dir, "scripts", "reference_wizard", "reference_wizard.R"))
source(file.path(repo_dir, "scripts", "trimm_and_rename_seq.R"))

fixture_dir <- file.path(script_dir, "fixtures")
output_dir <- tempfile("reference_wizard_")
dir.create(output_dir)
bundle_path <- file.path(output_dir, "synthetic_test.yaml")

summary <- capture.output(generate_bundle(list(
  organism = "Synthetic plant",
  assembly = "test_v1",
  annotation = file.path(fixture_dir, "genes.gff3"),
  chromosome_sizes = file.path(fixture_dir, "chromosome_sizes.tsv"),
  te = file.path(fixture_dir, "tes.tsv"),
  descriptions = file.path(fixture_dir, "descriptions.tsv"),
  centromeres = file.path(fixture_dir, "centromeres.bed"),
  gene_to_go = file.path(fixture_dir, "gene_to_go.tsv"),
  gene_sets = file.path(fixture_dir, "gene_sets.tsv"),
  kegg_organism = "syn",
  output = bundle_path
)))

bundle <- read_reference_bundle(bundle_path)
stopifnot(reference_bundle_summary(bundle) == "Synthetic plant / test_v1")
stopifnot(bundle_get(bundle, "annotation.format") == "gff3")
stopifnot(bundle_get(bundle, "annotation.gene_id_field") == "ID")
stopifnot(bundle_get(bundle, "functional.kegg_organism") == "syn")
stopifnot(file.exists(bundle_get(bundle, "genome.chromosome_sizes")))
stopifnot(grepl("synthetic_test[.]chromosome_sizes[.]tsv$",
               bundle_get(bundle, "genome.chromosome_sizes")))
stopifnot(any(grepl("^Bundle\t", summary)))
generated_genes <- harmonize_seqlevels(read_gene_annotation(bundle = bundle), bundle)
generated_tes <- harmonize_seqlevels(read_te_annotation(bundle = bundle), bundle)
stopifnot(length(generated_genes) == 4L, length(generated_tes) == 1L)
stopifnot(identical(as.character(seqlevels(generated_tes)), "linkage_A"))
stopifnot(all(generated_genes$gene_id == "GENE_A"))
generated_capabilities <- bundle_capabilities(bundle)
stopifnot(generated_capabilities[["gene_annotation"]])
stopifnot(generated_capabilities[["transposable_elements"]])
stopifnot(generated_capabilities[["go"]])
stopifnot(generated_capabilities[["kegg"]])

mapping <- map_annotation_to_reference(
  c("Chr01", "Chr02", "scaffold_A"),
  c("1", "2", "scaffold_A")
)
stopifnot(identical(unname(mapping), c("1", "2", "scaffold_A")))

kegg_cache <- file.path(output_dir, "kegg_organisms.tsv")
write.table(
  data.frame(
    code = c("syn", "oth"),
    name = c("Synthetic plant", "Other plant"),
    taxonomy = c("Eukaryotes;Plants", "Eukaryotes;Plants")
  ),
  kegg_cache, sep = "\t", quote = FALSE, row.names = FALSE
)
kegg_rows <- capture.output(list_kegg("Synthetic", cache = kegg_cache))
stopifnot(identical(kegg_rows, "syn\tSynthetic plant"))
invisible(capture.output(list_orgdb("Synthetic")))

no_size_path <- file.path(output_dir, "no_sizes.yaml")
no_size_summary <- capture.output(generate_bundle(list(
  organism = "Synthetic plant",
  assembly = "test_v1",
  annotation = file.path(fixture_dir, "genes.gff3"),
  output = no_size_path
)))
no_size_bundle <- read_reference_bundle(no_size_path)
stopifnot(is.null(bundle_get(no_size_bundle, "genome.chromosome_sizes")))
stopifnot(any(grepl("Exact chromosome lengths are unavailable", no_size_summary)))

fai_path <- file.path(output_dir, "from_fai.yaml")
fai_summary <- capture.output(generate_bundle(list(
  organism = "Synthetic plant",
  assembly = "test_v1",
  annotation = file.path(fixture_dir, "genes.gff3"),
  fasta = file.path(fixture_dir, "genome.fa.fai"),
  output = fai_path
)))
fai_bundle <- read_reference_bundle(fai_path)
stopifnot(bundle_chromosome_sizes(fai_bundle)[["linkage_A"]] == 10000)
stopifnot(bundle_get(fai_bundle, "generated.chromosome_length_source") == "FASTA index")

bad_sizes <- file.path(output_dir, "bad_sizes.tsv")
write.table(
  data.frame(seqname = "linkage_A", length = 200),
  bad_sizes, sep = "\t", quote = FALSE, row.names = FALSE
)
failed <- FALSE
tryCatch(
  generate_bundle(list(
    organism = "Synthetic plant",
    assembly = "wrong_assembly",
    annotation = file.path(fixture_dir, "genes.gff3"),
    chromosome_sizes = bad_sizes,
    output = file.path(output_dir, "bad.yaml")
  )),
  error = function(e) failed <<- grepl("same assembly", conditionMessage(e))
)
stopifnot(failed)

message("reference_wizard_smoke: OK")
