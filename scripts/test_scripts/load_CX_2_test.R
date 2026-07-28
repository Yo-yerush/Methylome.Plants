### load libreries
lib_packages <- c(
  "dplyr", "tidyr", "ggplot2", "DMRcaller", "rtracklayer", "lattice",
  "PeakSegDisk", "topGO", "KEGGREST", "Rgraphviz", "org.At.tair.db",
  "GenomicFeatures", "geomtextpath", "plyranges", "parallel",
  "RColorBrewer", "circlize", "cowplot", "knitr", "data.table"
)
for (n.pkg in seq(lib_packages)) {
  tryCatch(
    {
      suppressWarnings(suppressMessages(library(lib_packages[n.pkg], character.only = TRUE)))
      perc_val <- (n.pkg / length(lib_packages)) * 100
      cat(paste0("\rloading libraries [", round(perc_val, 1), "%] "))
    },
    error = function(e) {
      cat(paste0("\nError loading ", lib_packages[n.pkg], ": ", e$message, "\n"))
      message(paste0("\n* Error loading ", lib_packages[n.pkg], "package\n"))
    }
  )
}

Methylome.Plants_path <- "/home/yoyerush/yo/methylome_pipeline/Methylome.Plants_290126/Methylome.Plants/"
CX_files_dir_path <- "/home/yoyerush/yo/methylome_pipeline/other_mutants/stroud_et_al_2013/bismark_results/"

# configurations
var1 <- "wt"
var2 <- "suvh8"
var1_path <- c(
  paste0(CX_files_dir_path, "wt_2_bismark_se.CX_report.txt.gz"),
  paste0(CX_files_dir_path, "wt_3_bismark_se.CX_report.txt.gz")
)
var2_path <- paste0(CX_files_dir_path, "suvh8_bismark_se.CX_report.txt.gz")
annotation_file <- paste0(Methylome.Plants_path, "annotation_files/At_custom_annotations.csv.gz")
description_file <- paste0(Methylome.Plants_path, "annotation_files/At_custom_description_file.csv.gz")
TEs_file <- "annotation_files/TAIR10_Transposable_Elements.txt"
minProportionDiff = c(0.4, 0.2, 0.1) # CG, CHG, CHH
binSize = 100
minCytosinesCount = 4
minReadsPerCytosine = 6
pValueThreshold = 0.05
methyl_files_type = "CX_report"
img_type = "png"
n.cores = 8
analyze_DMRs = TRUE
run_PCA_plot = TRUE
run_total_meth_plot = TRUE
run_CX_Chrplot = TRUE
run_TEs_distance_n_size = TRUE
total_meth_annotation = TRUE
run_TF_motifs = TRUE
run_functional_groups = TRUE
run_GO_analysis = FALSE
run_KEGG_pathways = FALSE
analyze_strand_asymmetry_DMRs = TRUE
analyze_DMVs = TRUE
analyze_dH = FALSE
run_TE_metaPlots = TRUE
run_GeneBody_metaPlots = TRUE
run_GeneFeatures_metaPlots = TRUE
gene_features_binSize = 10
metaPlot.random.genes = 10000
comparison_name <- paste0(var2, "_vs_", var1)
exp_path <- paste0(Methylome.Plants_path, "results/", comparison_name)

###########################################################################

time_msg <<- function(suffix = " test:\t") paste0(format(Sys.time(), "[%H:%M]"), suffix)

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.Plants/main/scripts/trimm_and_rename_seq.R")

##### Read annotation and description files #####

annotation.gr <- read.csv(annotation_file) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  trimm_and_rename()

# TAIR10 Transposable Elements file

source("https://raw.githubusercontent.com/Yo-yerush/Methylome.Plants/main/scripts/edit_TE_file.R")
TEs_file.df <- read.csv(TEs_file, sep = "\t")
TE_gr <- edit_TE_file(TEs_file.df)


# upload description file

description_df <- read.csv(description_file)
names(description_df)[1] <- "gene_id"


###########################################################################

is_single <- (length(var1_path) == 1 & length(var2_path) == 1) # both genotypes includes 1 sample
is_Replicates <- (length(var1_path) > 1 & length(var2_path) > 1) # both genotypes includes >1 samples

var_args <- list(
  list(path = var1_path, name = var1),
  list(path = var2_path, name = var2)
)

##### load methylation data ('CX_report' file) #####
message("load CX methylation data...")
source("https://raw.githubusercontent.com/Yo-yerush/Methylome.Plants/main/scripts/load_replicates.R")
n.cores <- 10
methyl_files_type <- "CX_report"

# load 'CX_reports'
n.cores.load <- ifelse(n.cores > 1, 2, 1)
load_vars <- mclapply(var_args, function(x) load_replicates(x$path, n.cores, x$name, F, methyl_files_type), mc.cores = n.cores.load)

# trimm seqs objects (rename if not 'TAIR10' Chr seqnames)
meth_var1 <- trimm_and_rename(load_vars[[1]]$methylationData_pool)
meth_var2 <- trimm_and_rename(load_vars[[2]]$methylationData_pool)
meth_var1_replicates <- trimm_and_rename(load_vars[[1]]$methylationDataReplicates)
meth_var2_replicates <- trimm_and_rename(load_vars[[2]]$methylationDataReplicates)

methylationDataReplicates_joints <- joinReplicates(
  meth_var1,
  meth_var2
)
