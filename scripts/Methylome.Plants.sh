#!/bin/bash

samples_file=""

# Get the directory where the Bash script is located
Methylome_At_path=$(pwd)
# Methylome_At_path=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
# Methylome_At_path="$Methylome_At_path/../"

# Default files path
annotation_file=""
description_file=""
TEs_file=""
reference_bundle="$Methylome_At_path/reference_bundles/arabidopsis_thaliana_TAIR10.yaml"
reference_bundle_explicit=false

# Default values for optional arguments
minProportionDiff_CG=0.4
minProportionDiff_CHG=0.2
minProportionDiff_CHH=0.1
binSize=100
minCytosinesCount=4
minReadsPerCytosine=6
pValueThreshold=0.05
methyl_files_type=CX_report
img_type=png
n_cores=8

DMR_analysis=TRUE
QC_plots=TRUE
pca=FALSE
total_methylation=FALSE
CX_ChrPlot=FALSE
TEs_dis_n_size=FALSE
total_meth_ann=FALSE
TF_motifs=FALSE
func_groups=FALSE
GO_analysis=FALSE
KEGG_pathways=FALSE
strand_DMRs=FALSE
DMV_analysis=FALSE
dH_scripts=FALSE
TEs_mp=FALSE
Genes_mp=FALSE
Gene_features_mp=FALSE
bin_size_features=10
metaPlot_random_genes=10000

# Function to display help text
usage() {
  echo ""
  echo "Usage: $0 [samples_file] [options]"
  echo ""
  echo "Required argument:"
  echo "  --samples_file                Path to samples file [required]"
  echo ""
  echo "Optional arguments:"
  echo "  --minReadsPerCytosine         Minimum reads per cytosine [default: $minReadsPerCytosine]"
  echo "  --n_cores                     Number of cores [default: $n_cores]"
  echo "  --image_type                  Output images format [default: '$img_type']"
  echo "  --file_type                   Post-alignment file type - 'CX_report', 'bedMethyl' and 'CGmap' [default: '$methyl_files_type' OR determine automatically]"
  echo "  --annotation_file             Override the bundle gene annotation file"
  echo "  --description_file            Override the bundle description file"
  echo "  --TEs_file                    Override the bundle Transposable Elements file"
  echo "  --reference_bundle            Species/assembly reference bundle YAML [default: bundled TAIR10 profile]"
  echo "  --pipeline_path               Path to Methylome.Plants [default: $Methylome_At_path]"
  echo "  --Methylome_At_path           Deprecated alias for --pipeline_path"
  echo ""
  echo "DMRs analysis arguments:"
  echo "  --minProportionDiff_CG        Minimum proportion difference for CG [default: $minProportionDiff_CG]"
  echo "  --minProportionDiff_CHG       Minimum proportion difference for CHG [default: $minProportionDiff_CHG]"
  echo "  --minProportionDiff_CHH       Minimum proportion difference for CHH [default: $minProportionDiff_CHH]"
  echo "  --binSize                     DMRs bin size [default: $binSize]"
  echo "  --minCytosinesCount           Minimum cytosines count [default: $minCytosinesCount]"
  echo "  --pValueThreshold             P-value (padj) threshold [default: $pValueThreshold]"
  echo ""
  echo "  --QC_off                      Disable sample-level QC plots (distributions, scatter, correlation, variance)"
  echo "  --pca                         Perform PCA for total methylation levels"
  echo "  --total_methylation           Total methylation bar-plot"
  echo "  --CX_ChrPlot                  total methylation chromosome plot"
  echo "  --TEs_distance_n_size         Analyze TEs total methylation by size and distance from centromere"
  echo "  --total_meth_ann              Total methylation per genic annotations"
  echo "  --TF_motifs                   Transcription factors motif analysis"
  echo "  --func_groups                 Functional groups genes overlap DMRs"
  echo "  --GO_analysis                 Perform GO analysis over DMRs"
  echo "  --KEGG_pathways               Perform KEGG pathways analysis over DMRs"
  echo "  --all_analyses                Enable all analyses (sets all analysis flags above [pca --> KEGG_pathways])"
  echo ""
  echo "  --DMRs_off                    Disable the main DMRs analysis workflow"
  echo "  --strand_DMRs                 Analyze strand-specific (+/-) DMRs"
  echo "  --DMVs                        Analyze differentially methylated vallies (1kbp)"
  echo "  --dH                          instead of DMRs worflow (calculated by ratios [p]), analyze SurpDMRs:"
  echo "                                delta-H = -(p * log2(p) + (1 - p) * log2(1 - p))"
  echo ""
  echo "MetaPlots analysis arguments:"
  echo "  --MP_TEs                      Analyze of TEs metaPlots"
  echo "  --MP_Genes                    Analyze of Genes-body metaPlots"
  echo "  --MP_Gene_features            Analyze Gene Features metaPlots"
  echo "  --all_metaplots               Enable all metaPlots analyses"
  echo "  --MP_features_bin_size        Bin-size (set only for 'Gene_features' analysis!) [default: $bin_size_features]"
  echo "  --metaPlot_random             Number of random genes/TEs for metaPlots. 'all' for all the coding-genes and TEs [default: $metaPlot_random_genes]"
  echo ""
  exit 1
}

# Assign first positional argument to samples_file if it does not start with --
if [[ -n "$1" && "$1" != --* ]]; then
    samples_file="$1"
    shift
fi

# Check if the first argument is 'help' or a help flag
#if [[ "$1" == "-h" || "$1" == "--help" ]]; then
#  usage
#  exit 1
#fi

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --samples_file) samples_file="$2"; shift ;;
    --minProportionDiff_CG) minProportionDiff_CG="$2"; shift ;;
    --minProportionDiff_CHG) minProportionDiff_CHG="$2"; shift ;;
    --minProportionDiff_CHH) minProportionDiff_CHH="$2"; shift ;;
    --binSize) binSize="$2"; shift ;;
    --minCytosinesCount) minCytosinesCount="$2"; shift ;;
    --minReadsPerCytosine) minReadsPerCytosine="$2"; shift ;;
    --pValueThreshold) pValueThreshold="$2"; shift ;;
    --file_type) methyl_files_type="$2"; shift ;;
    --image_type) img_type="$2"; shift ;;
    --n_cores) n_cores="$2"; shift ;;
    --annotation_file) annotation_file="$2"; shift ;;
    --description_file) description_file="$2"; shift ;;
    --TEs_file) TEs_file="$2"; shift ;;
    --reference_bundle) reference_bundle="$2"; reference_bundle_explicit=true; shift ;;
    --pipeline_path | --Methylome_At_path) Methylome_At_path="$2"; shift ;;
    --DMRs_off) DMR_analysis=FALSE ;;
    --QC_off) QC_plots=FALSE ;;
    --strand_DMRs) strand_DMRs=TRUE ;;
    --DMVs) DMV_analysis=TRUE ;;
    --dH) dH_scripts=TRUE ;;
    --pca) pca=TRUE ;;
    --total_methylation) total_methylation=TRUE ;;
    --CX_ChrPlot) CX_ChrPlot=TRUE ;;
    --TEs_distance_n_size) TEs_dis_n_size=TRUE ;;
    --total_meth_ann) total_meth_ann=TRUE ;;
    --TF_motifs) TF_motifs=TRUE ;;
    --func_groups) func_groups=TRUE ;;
    --GO_analysis) GO_analysis=TRUE ;;
    --KEGG_pathways) KEGG_pathways=TRUE ;;
    --all_analyses)
      pca=TRUE
      total_methylation=TRUE
      CX_ChrPlot=TRUE
      TEs_dis_n_size=TRUE
      total_meth_ann=TRUE
      TF_motifs=TRUE
      func_groups=TRUE
      GO_analysis=TRUE
      KEGG_pathways=TRUE
      ;;
    --MP_TEs) TEs_mp=TRUE ;;
    --MP_Genes) Genes_mp=TRUE ;;
    --MP_Gene_features) Gene_features_mp=TRUE ;;
    --all_metaplots)
      TEs_mp=TRUE
      Genes_mp=TRUE
      Gene_features_mp=TRUE
      ;;
    --MP_features_bin_size) bin_size_features="$2"; shift ;;
    --metaPlot_random) metaPlot_random_genes="$2"; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
  shift
done

if [[ "$reference_bundle_explicit" == "false" ]]; then
  reference_bundle="$Methylome_At_path/reference_bundles/arabidopsis_thaliana_TAIR10.yaml"
fi

# Check for required argument
if [ -z "$samples_file" ]; then
  echo "Error: --samples_file is a required argument."
  usage
fi

if [ ! -f "$reference_bundle" ]; then
  echo "Error: Reference bundle '$reference_bundle' does not exist."
  exit 1
fi

for resource_file in "$annotation_file" "$description_file" "$TEs_file"; do
  if [[ -n "$resource_file" && ! -f "$resource_file" ]]; then
    echo "Error: Override resource '$resource_file' does not exist."
    exit 1
  fi
done

# Preserve user-supplied relative paths before changing into the pipeline directory.
samples_file=$(readlink -f "$samples_file")
reference_bundle=$(readlink -f "$reference_bundle")
[ -n "$annotation_file" ] && annotation_file=$(readlink -f "$annotation_file")
[ -n "$description_file" ] && description_file=$(readlink -f "$description_file")
[ -n "$TEs_file" ] && TEs_file=$(readlink -f "$TEs_file")

# ensure files unix line endings [sed -i 's/\r$//' "$samples_file"]
dos2unix "$samples_file" 2>/dev/null
[ -n "$annotation_file" ] && dos2unix "$annotation_file" 2>/dev/null
[ -n "$description_file" ] && dos2unix "$description_file" 2>/dev/null
[ -n "$TEs_file" ] && dos2unix "$TEs_file" 2>/dev/null

# Change to the Methylome_At_path directory
cd "$Methylome_At_path" || {
    echo "Error: Cannot change directory to $Methylome_At_path"
    exit 1
}

# Extract samples unique values from the first column
is_txt=$(awk -F'\t' 'NR==1 {print NF}' "$samples_file")
is_csv=$(awk -F',' 'NR==1 {print NF}' "$samples_file")
if [ "$is_txt" -eq 2 ]; then
  samples_var=$(cut -f1 "$samples_file" | uniq)
  samples_path=$(cut -f2 "$samples_file" | head -n 1)
elif [ "$is_csv" -eq 2 ]; then
  samples_var=$(cut -d',' -f1 "$samples_file" | uniq)
  samples_path=$(cut -d',' -f2 "$samples_file" | head -n 1)
else
  echo "Error: Cannot read the samples file. Please provide a tab or comma separated file with two columns and no headers."
  exit 1
fi

# # Determine file type based on samples_path content
# if [[ "$methyl_files_type" != "CX_report" ]]; then 
#   if [[ $(basename "$samples_path" | awk -F. '{print $NF}') == "bed" ]]; then
#     methyl_files_type="bedMethyl"
#   elif [[ $(basename "$samples_path" .gz | awk -F. '{print $NF}') == "CGmap" ]]; then
#     methyl_files_type="CGmap"
#   else
#     echo "Error: Cannot read the samples file. Please provide a tab or comma separated file with two columns and no headers."
#     exit 1
#   fi
# fi

# Use head and sed to get the first and second lines
control_s=$(echo "$samples_var" | head -n1)
treatment_s=$(echo "$samples_var" | head -n2 | tail -n1)
condition_s="${treatment_s}_vs_${control_s}"

# Output the configuration
echo ""
echo "-------------------------------------"
echo ""
echo "$treatment_s VS. $control_s"
echo ""
echo "DMRs Min Proportion Diff CG: $minProportionDiff_CG"
echo "DMRs Min Proportion Diff CHG: $minProportionDiff_CHG"
echo "DMRs Min Proportion Diff CHH: $minProportionDiff_CHH"
echo "DMRs Bin size: $binSize"
echo "DMRs Min Cytosines Count: $minCytosinesCount"
echo "DMRs P-value Threshold: $pValueThreshold"
echo ""
echo "Min Reads Per Cytosine: $minReadsPerCytosine"
echo "Number of Cores: $n_cores"
echo "Post-alignment file type: $methyl_files_type"
echo "Output images format: $img_type"
echo ""
echo "QC plots: $QC_plots"
echo "PCA: $pca"
echo "Total methylation bar-plot: $total_methylation"
echo "Total methylation chromosome plot: $CX_ChrPlot"
echo "Analyze TEs total methylation by size and distance from centromere: $TEs_dis_n_size"
echo "Total methylation per genic annotations: $total_meth_ann"
echo "Transcription factors motif analysis: $TF_motifs"
echo "Functional groups genes overlap DMRs: $func_groups"
echo "GO Analysis: $GO_analysis"
echo "KEGG Pathways: $KEGG_pathways"
echo ""
echo "Gene-body metaPlots: $Genes_mp"
echo "TEs metaPlots: $TEs_mp"
echo "Gene features metaPlots: $Gene_features_mp"
echo "Bin size for gene features: $bin_size_features"
echo "Number of random genes/TEs for metaPlots: $metaPlot_random_genes"
echo ""
echo "Samples file: $samples_file"
echo "Annotation file: $annotation_file"
echo "Description file: $description_file"
echo "Transposable Elements file: $TEs_file"
echo "Reference bundle: $reference_bundle"
echo "Methylome.Plants directory path: $Methylome_At_path"
echo ""
echo "Analyze DMRs workflow: $DMR_analysis"
echo "Analyze strand-specific DMRs: $strand_DMRs"
echo "Analyze differentially methylated vallies (1kbp): $DMV_analysis"
echo "Analyze delta-H scripts: $dH_scripts"
echo ""
echo "-------------------------------------"
echo ""

# create results directory
mkdir -p results
mkdir -p results/${condition_s}

# Generate log file with a timestamp
log_file="results/${condition_s}/${condition_s}_$(date +"%d-%m-%y").log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "**  $treatment_s VS. $control_s" >> "$log_file"
echo "" >> "$log_file"

# Call the R script with the arguments
Rscript ./scripts/Methylome.Plants_run.R \
"$samples_file" \
"$Methylome_At_path" \
"$annotation_file" \
"$description_file" \
"$TEs_file" \
"$minProportionDiff_CG" \
"$minProportionDiff_CHG" \
"$minProportionDiff_CHH" \
"$binSize" \
"$minCytosinesCount" \
"$minReadsPerCytosine" \
"$pValueThreshold" \
"$methyl_files_type" \
"$img_type" \
"$n_cores" \
"$DMR_analysis" \
"$pca" \
"$QC_plots" \
"$total_methylation" \
"$CX_ChrPlot" \
"$TEs_dis_n_size" \
"$total_meth_ann" \
"$TF_motifs" \
"$func_groups" \
"$GO_analysis" \
"$KEGG_pathways" \
"$strand_DMRs" \
"$DMV_analysis" \
"$dH_scripts" \
"$TEs_mp" \
"$Genes_mp" \
"$Gene_features_mp" \
"$bin_size_features" \
"$metaPlot_random_genes" \
"$reference_bundle" \
    2>> "$log_file"

# Output again the configurations, now to the 'log' file
echo "" >> "$log_file"
echo "" >> "$log_file"
echo "" >> "$log_file"
echo "configuration:" >> "$log_file"
echo "--------------" >> "$log_file"
echo "DMRs Min Proportion Diff CG: $minProportionDiff_CG" >> "$log_file"
echo "DMRs Min Proportion Diff CHG: $minProportionDiff_CHG" >> "$log_file"
echo "DMRs Min Proportion Diff CHH: $minProportionDiff_CHH" >> "$log_file"
echo "DMRs Bin size: $binSize" >> "$log_file"
echo "DMRs Min Cytosines Count: $minCytosinesCount" >> "$log_file"
echo "DMRs P-value Threshold: $pValueThreshold" >> "$log_file"
echo "Number of Cores: $n_cores" >> "$log_file"
echo "Min Reads Per Cytosine: $minReadsPerCytosine" >> "$log_file"
echo "Post-alignment file type: $methyl_files_type" >> "$log_file"
echo "Output images format: $img_type" >> "$log_file"
echo "QC plots: $QC_plots" >> "$log_file"
echo "PCA: $pca" >> "$log_file"
echo "Total methylation bar-plot: $total_methylation" >> "$log_file"
echo "Total methylation chromosome plot: $CX_ChrPlot" >> "$log_file"
echo "Analyze TEs total methylation by size and distance from centromere: $TEs_dis_n_size" >> "$log_file"
echo "Total methylation per genic annotations: $total_meth_ann" >> "$log_file"
echo "Transcription factors motif analysis: $TF_motifs" >> "$log_file"
echo "Functional groups genes overlap DMRs: $func_groups" >> "$log_file"
echo "GO Analysis: $GO_analysis" >> "$log_file"
echo "KEGG Pathways: $KEGG_pathways" >> "$log_file"
echo "Analyze DMRs workflow: $DMR_analysis" >> "$log_file"
echo "Analyze strand-specific DMRs: $strand_DMRs" >> "$log_file"
echo "Analyze differentially methylated vallies (1kbp): $DMV_analysis" >> "$log_file"
echo "Analyze delta-H scripts: $dH_scripts" >> "$log_file"
echo "Gene-body metaPlots: $Genes_mp" >> "$log_file"
echo "TEs metaPlots: $TEs_mp" >> "$log_file"
echo "Gene features metaPlots: $Gene_features_mp" >> "$log_file"
echo "Bin size for gene features: $bin_size_features" >> "$log_file"
echo "Number of random genes/TEs for metaPlots: $metaPlot_random_genes" >> "$log_file"
echo "Samples file: $samples_file" >> "$log_file"
echo "Annotation file: $annotation_file" >> "$log_file"
echo "Description file: $description_file" >> "$log_file"
echo "Transposable Elements file: $TEs_file" >> "$log_file"
echo "Reference bundle: $reference_bundle" >> "$log_file"
echo "Methylome_At_path: $Methylome_At_path" >> "$log_file"
