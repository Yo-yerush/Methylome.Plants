#!/usr/bin/env bash

# Get the directory where the Bash script is located
Methylome_At_path=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
cd "$Methylome_At_path"

#####################################
# Initialize and Activate Conda Env
#####################################
if [ "$CONDA_DEFAULT_ENV" != "Methylome.Plants_env" ]; then
  eval "$(conda shell.bash hook)"
  conda activate Methylome.Plants_env
  
  if [ "$CONDA_DEFAULT_ENV" != "Methylome.Plants_env" ]; then
    echo "Error: Failed to activate the 'Methylome.Plants_env' Conda environment."
    echo "Please activate it manually using 'conda activate Methylome.Plants_env' and rerun the script."
    exit 1
  fi
fi

###############
# CONFIGURATION
###############
# Default parameters for run_bismark.sh:
SCRIPT_BIS_DEFAULT_genome="TAIR10"
SCRIPT_BIS_DEFAULT_ncores="8"

# Default parameters for Methylome.Plants.sh:
SCRIPT1_DEFAULT_minProportionDiff_CG="0.4"
SCRIPT1_DEFAULT_minProportionDiff_CHG="0.2"
SCRIPT1_DEFAULT_minProportionDiff_CHH="0.1"
SCRIPT1_DEFAULT_binSize="100"
SCRIPT1_DEFAULT_minCytosinesCount="4"
SCRIPT1_DEFAULT_minReadsPerCytosine="6"
SCRIPT1_DEFAULT_pValueThreshold="0.05"
SCRIPT1_DEFAULT_n_cores="8"
SCRIPT1_DEFAULT_pca="yes"
SCRIPT1_DEFAULT_QC_plots="yes"
SCRIPT1_DEFAULT_total_methylation="yes"
SCRIPT1_DEFAULT_CX_ChrPlot="yes"
SCRIPT1_DEFAULT_TEs_distance_n_size="yes"
SCRIPT1_DEFAULT_total_meth_ann="yes"
SCRIPT1_DEFAULT_TF_motifs="yes"
SCRIPT1_DEFAULT_func_groups="yes"
SCRIPT1_DEFAULT_GO_analysis="no"
SCRIPT1_DEFAULT_KEGG_pathways="no"
SCRIPT1_DEFAULT_file_type="CX_report"
SCRIPT1_DEFAULT_img_type="png"
SCRIPT1_DEFAULT_annotation_file=""
SCRIPT1_DEFAULT_description_file=""
SCRIPT1_DEFAULT_TEs_file=""
SCRIPT1_DEFAULT_reference_bundle="reference_bundles/arabidopsis_thaliana_TAIR10.yaml"
SCRIPT1_DEFAULT_disable_DMRs="no"
SCRIPT1_DEFAULT_strand_DMRs="no"
SCRIPT1_DEFAULT_DMVs="no"
SCRIPT1_DEFAULT_delta_H="no"
SCRIPT1_DEFAULT_TEs_metaplots="yes"
SCRIPT1_DEFAULT_Genes_metaplots="yes"
SCRIPT1_DEFAULT_Gene_features_metaplots="no"
SCRIPT1_DEFAULT_bin_size_features="10"
SCRIPT1_DEFAULT_metaPlot_random_genes="10000"

# Paths to the scripts we want to run (adjust if needed)
SCRIPT_BIS_PATH="./scripts/run_bismark.sh"
SCRIPT1_PATH="./scripts/Methylome.Plants.sh"

##################
# WHIPTAIL DIALOGS
##################

# Prompt user: which scripts do you want to run?
CHOICE=$(whiptail --title "Choose scripts to run" \
  --checklist "Select which pipeline(s) to run. Use SPACE to toggle selection, ENTER to confirm, ESC to cancle." \
  18 70 3 \
  "Bismark" "Run genome alignment with Bismark" OFF \
  "Methylome.Plants" "Run main methylome pipeline'" ON \
  3>&1 1>&2 2>&3)

# If user hits Cancel or ESC, exit
if [ $? -ne 0 ]; then
  echo "No scripts selected. Exiting."
  exit 1
fi

# Convert whiptail’s checklist output into an array
SELECTED_SCRIPTS=()
for item in $CHOICE; do
  # Remove surrounding quotes
  item=$(echo $item | sed 's/"//g')
  SELECTED_SCRIPTS+=("$item")
done

#############################
# Prompt for the samples_file
#############################
SAMPLES_FILE=$(whiptail --title "samples_file" --inputbox \
  "Enter the path to the samples table file:" \
  10 70 \
  "" \
  3>&1 1>&2 2>&3)
SAMPLES_FILE="${SAMPLES_FILE//\"/}"
if [ $? -ne 0 ] || [ -z "$SAMPLES_FILE" ]; then
  echo "samples_file is required. Exiting."
  exit 1
fi

########################################
# Function to edit parameters interactively
########################################

# A generic function to create a parameter menu for bismrk scripts
edit_script_bis_parameters() {
  # Parameters are expected to be set before calling this function
  while true; do
    OPTION=$(whiptail --title "'Bismark alignment' Parameters" --menu "Select a parameter to change or proceed with current settings." 25 78 16 \
      "Proceed." "Use current parameters" \
      "$SCRIPT_BIS_DEFAULT_genome" "Reference FASTA file" \
      "$SCRIPT_BIS_DEFAULT_ncores" "Number of cores" \
      3>&1 1>&2 2>&3)

    # Check if user cancelled
    [ $? -ne 0 ] && return 1

    if [ "$OPTION" = "Proceed." ]; then
      break
    elif [ "$OPTION" = "$SCRIPT_BIS_DEFAULT_genome" ]; then
      SCRIPT_BIS_DEFAULT_genome=$(whiptail --inputbox "Path to reference genome file" 10 70 "$SCRIPT_BIS_DEFAULT_genome" 3>&1 1>&2 2>&3 || echo "$SCRIPT_BIS_DEFAULT_genome")
    elif [ "$OPTION" = "$SCRIPT_BIS_DEFAULT_ncores" ]; then
      SCRIPT_BIS_DEFAULT_ncores=$(whiptail --inputbox "Number of cores" 10 70 "$SCRIPT_BIS_DEFAULT_ncores" 3>&1 1>&2 2>&3 || echo "$SCRIPT_BIS_DEFAULT_ncores")
    fi
  done
}

# A generic function to create a parameter menu for script1
edit_script1_parameters() {
    # nice aligned display helpers (monospace)
    fmt() { printf '%-10s  %-62s' "$1" "$2"; }

  # Parameters are expected to be set before calling this function
  while true; do
        OPTION=$(whiptail --title "'Methylome.Plants' Parameters" \
            --menu "Select a parameter to change or proceed with current settings." 40 50 33 \
            "Proceed."                "$(fmt '' 'Use current parameters')" \
            "Set off"       "$(fmt '' 'Turn OFF all analyses')" \
            "Min diff (CG)"           "$(fmt "$SCRIPT1_minProportionDiff_CG" 'Min methylation proportion difference to call CG DMRs')" \
            "Min diff (CHG)"          "$(fmt "$SCRIPT1_minProportionDiff_CHG" 'Min methylation proportion difference to call CHG DMRs')" \
            "Min diff (CHH)"          "$(fmt "$SCRIPT1_minProportionDiff_CHH" 'Min methylation proportion difference to call CHH DMRs')" \
            "DMR bin size (bp)"       "$(fmt "$SCRIPT1_binSize" 'Bin-size window used for DMR calling')" \
            "Min cytosines / bin"     "$(fmt "$SCRIPT1_minCytosinesCount" 'Minimum cytosine count required in a bin to test it')" \
            "Min reads / cytosine"    "$(fmt "$SCRIPT1_minReadsPerCytosine" 'Minimum read coverage per cytosine to include it')" \
            "Adj. p-value cutoff"     "$(fmt "$SCRIPT1_pValueThreshold" 'Adjusted p-value threshold for calling significant DMRs')" \
            "CPU cores"               "$(fmt "$SCRIPT1_n_cores" 'Number of cores for parallel steps')" \
            "Input file format"       "$(fmt "$SCRIPT1_file_type" 'Methylation input format (CX_report/bedMethyl/CGmap)')" \
            "Figure format"           "$(fmt "$SCRIPT1_img_type" 'Output image format for plots (pdf/svg/png/tiff/...)')" \
            "PCA"                     "$(fmt "$SCRIPT1_pca" 'Perform PCA analysis')" \
            "QC plots"                "$(fmt "$SCRIPT1_QC_plots" 'Sample-level QC plots (distributions, scatter, correlation, variance)')" \
            "Total methylation bar-plots"  "$(fmt "$SCRIPT1_total_methylation" 'Produce total methylation bar-plots')" \
            "CX Chromosome Plot"      "$(fmt "$SCRIPT1_CX_ChrPlot" 'Generate chromosome-wide CX plots')" \
            "TEs distance and size"   "$(fmt "$SCRIPT1_TEs_distance_n_size" 'Analyze TEs total methylation by size and distance from centromere')" \
            "Total methylation annotations" "$(fmt "$SCRIPT1_total_meth_ann" 'Total methylation per genic annotations')" \
            "TF motifs"               "$(fmt "$SCRIPT1_TF_motifs" 'Analyze transcription factor motifs')" \
            "Functional groups"       "$(fmt "$SCRIPT1_func_groups" 'Functional groups genes overlap DMRs')" \
            "GO enrichment"           "$(fmt "$SCRIPT1_GO_analysis" 'GO enrichment on DMR-associated gene-bodies/promoters')" \
            "KEGG enrichment"         "$(fmt "$SCRIPT1_KEGG_pathways" 'KEGG pathway enrichment on DMR-associated gene-bodies/promoters')" \
            "TE meta-plots"           "$(fmt "$SCRIPT1_TEs_metaplots" 'Metaplots over transposable elements (TE bodies/flanks)')" \
            "Gene-body meta-plots"    "$(fmt "$SCRIPT1_Genes_metaplots" 'Metaplots across gene bodies (TSS→TES/flank)')" \
            "Gene-feature meta-plots" "$(fmt "$SCRIPT1_Gene_features_metaplots" 'Metaplots over gene features (promoter/CDS/intron/UTR)')" \
            "Feature bin size"        "$(fmt "$SCRIPT1_bin_size_features" 'Bins-size for each gene feature region in feature meta-plots')" \
            "Random genes"            "$(fmt "$SCRIPT1_metaPlot_random_genes" 'Number of genes to sample for meta-plots (or all)')" \
            "Disable DMRs analysis"   "$(fmt "$SCRIPT1_disable_DMRs" 'Disable the main DMRs analysis workflow')" \
            "Strand-specific DMRs"    "$(fmt "$SCRIPT1_strand_DMRs" 'Analyze strand-specific (+/-) DMRs')" \
            "DMVs analysis"           "$(fmt "$SCRIPT1_DMVs" 'Analyze differentially methylated vallies (1kbp)')" \
            "dH analysis"             "$(fmt "$SCRIPT1_delta_H" 'instead of DMRs worflow (calculated by ratios [p]), analyze SurpDMRs')" \
            "Annotation file (gtf/gff/csv)"  "$(fmt "$SCRIPT1_annotation_file" '')" \
            "Gene descriptions (txt/csv)"    "$(fmt "$SCRIPT1_description_file" '')" \
            "TE annotation file (txt)"       "$(fmt "$SCRIPT1_TEs_file" '')" \
            "Reference bundle (yaml)"        "$(fmt "$SCRIPT1_reference_bundle" '')" \
            3>&1 1>&2 2>&3)

    # Check if user cancelled
    # Check if user cancelled
    [ $? -ne 0 ] && return 1

    case "$OPTION" in
      "Proceed.")
        break
        ;;

      "Set off")
        if whiptail --yesno "Turn OFF all optional analyses?" 18 70; then
          SCRIPT1_QC_plots="no"
          SCRIPT1_pca="no"
          SCRIPT1_total_methylation="no"
          SCRIPT1_CX_ChrPlot="no"
          SCRIPT1_TEs_distance_n_size="no"
          SCRIPT1_total_meth_ann="no"
          SCRIPT1_TF_motifs="no"
          SCRIPT1_func_groups="no"
          SCRIPT1_GO_analysis="no"
          SCRIPT1_KEGG_pathways="no"
          SCRIPT1_disable_DMRs="yes"
          SCRIPT1_strand_DMRs="no"
          SCRIPT1_DMVs="no"
          SCRIPT1_delta_H="no"
          SCRIPT1_TEs_metaplots="no"
          SCRIPT1_Genes_metaplots="no"
          SCRIPT1_Gene_features_metaplots="no"
        fi
        ;;

      "Min diff (CG)")
        SCRIPT1_minProportionDiff_CG=$(
          whiptail --inputbox "Minimum methylation proportion difference for CG DMRs" 10 80 \
            "$SCRIPT1_minProportionDiff_CG" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minProportionDiff_CG="$SCRIPT1_minProportionDiff_CG"
        ;;

      "Min diff (CHG)")
        SCRIPT1_minProportionDiff_CHG=$(
          whiptail --inputbox "Minimum methylation proportion difference for CHG DMRs" 10 80 \
            "$SCRIPT1_minProportionDiff_CHG" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minProportionDiff_CHG="$SCRIPT1_minProportionDiff_CHG"
        ;;

      "Min diff (CHH)")
        SCRIPT1_minProportionDiff_CHH=$(
          whiptail --inputbox "Minimum methylation proportion difference for CHH DMRs" 10 80 \
            "$SCRIPT1_minProportionDiff_CHH" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minProportionDiff_CHH="$SCRIPT1_minProportionDiff_CHH"
        ;;

      "DMR bin size (bp)")
        SCRIPT1_binSize=$(
          whiptail --inputbox "Bin-size (bp) for DMR calling" 10 80 \
            "$SCRIPT1_binSize" 3>&1 1>&2 2>&3
        ) || SCRIPT1_binSize="$SCRIPT1_binSize"
        ;;

      "Min cytosines / bin")
        SCRIPT1_minCytosinesCount=$(
          whiptail --inputbox "Minimum cytosines count per bin" 10 80 \
            "$SCRIPT1_minCytosinesCount" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minCytosinesCount="$SCRIPT1_minCytosinesCount"
        ;;

      "Min reads / cytosine")
        SCRIPT1_minReadsPerCytosine=$(
          whiptail --inputbox "Minimum reads per cytosine" 10 80 \
            "$SCRIPT1_minReadsPerCytosine" 3>&1 1>&2 2>&3
        ) || SCRIPT1_minReadsPerCytosine="$SCRIPT1_minReadsPerCytosine"
        ;;

      "Adj. p-value cutoff")
        SCRIPT1_pValueThreshold=$(
          whiptail --inputbox "Adjusted p-value threshold for DMRs" 10 80 \
            "$SCRIPT1_pValueThreshold" 3>&1 1>&2 2>&3
        ) || SCRIPT1_pValueThreshold="$SCRIPT1_pValueThreshold"
        ;;

      "CPU cores")
        SCRIPT1_n_cores=$(
          whiptail --inputbox "Number of cores" 10 80 \
            "$SCRIPT1_n_cores" 3>&1 1>&2 2>&3
        ) || SCRIPT1_n_cores="$SCRIPT1_n_cores"
        ;;

      "Input file format")
        SCRIPT1_file_type=$(
          whiptail --radiolist "Select methylation file type:" 15 80 3 \
            "CX_report" "Bismark CX_report (.txt)"   $([ "$SCRIPT1_file_type" = "CX_report" ] && echo ON || echo OFF) \
            "bedMethyl" "BED methylation (.bed)"     $([ "$SCRIPT1_file_type" = "bedMethyl" ] && echo ON || echo OFF) \
            "CGmap"     "CGmap format (.CGmap)"      $([ "$SCRIPT1_file_type" = "CGmap" ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_file_type="$SCRIPT1_file_type"
        ;;

      "Figure format")
        SCRIPT1_img_type=$(
          whiptail --radiolist "Select image format:" 18 80 6 \
            "pdf"  "Vector; publication-friendly"                     $([ "$SCRIPT1_img_type" = "pdf"  ] && echo ON || echo OFF) \
            "svg"  "Vector; editable"                                 $([ "$SCRIPT1_img_type" = "svg"  ] && echo ON || echo OFF) \
            "png"  "Raster (lossless)"                                $([ "$SCRIPT1_img_type" = "png"  ] && echo ON || echo OFF) \
            "tiff" "Raster (lossless - LZW); journal standard"        $([ "$SCRIPT1_img_type" = "tiff" ] && echo ON || echo OFF) \
            "jpeg" "Raster (lossy); small file"                       $([ "$SCRIPT1_img_type" = "jpeg" ] && echo ON || echo OFF) \
            "bmp"  "Raster (uncompressed); avoid"                     $([ "$SCRIPT1_img_type" = "bmp"  ] && echo ON || echo OFF) \
            3>&1 1>&2 2>&3
        ) || SCRIPT1_img_type="$SCRIPT1_img_type"
        ;;

      "GO enrichment")
        if whiptail --yesno "Perform GO enrichment analysis?" 10 60; then
          SCRIPT1_GO_analysis="yes"
        else
          SCRIPT1_GO_analysis="no"
        fi
        ;;

      "KEGG enrichment")
        if whiptail --yesno "Perform KEGG pathway enrichment analysis?" 10 60; then
          SCRIPT1_KEGG_pathways="yes"
        else
          SCRIPT1_KEGG_pathways="no"
        fi
        ;;

      "TE meta-plots")
        if whiptail --yesno "Generate metaplots for Transposable Elements (TEs)?" 10 60; then
          SCRIPT1_TEs_metaplots="yes"
        else
          SCRIPT1_TEs_metaplots="no"
        fi
        ;;

      "Gene-body meta-plots")
        if whiptail --yesno "Generate metaplots across gene bodies?" 10 60; then
          SCRIPT1_Genes_metaplots="yes"
        else
          SCRIPT1_Genes_metaplots="no"
        fi
        ;;

      "Gene-feature meta-plots")
        if whiptail --yesno "Generate metaplots over gene features (promoter/CDS/UTRs/introns)?" 10 60; then
          SCRIPT1_Gene_features_metaplots="yes"
        else
          SCRIPT1_Gene_features_metaplots="no"
        fi
        ;;

      "Feature bin size")
        SCRIPT1_bin_size_features=$(
          whiptail --inputbox "Bin size for feature metaplots" 10 80 \
            "$SCRIPT1_bin_size_features" 3>&1 1>&2 2>&3
        ) || SCRIPT1_bin_size_features="$SCRIPT1_bin_size_features"
        ;;

      "Random genes")
        SCRIPT1_metaPlot_random_genes=$(
          whiptail --inputbox "Number of random genes for metaplots (or 'all')" 10 80 \
            "$SCRIPT1_metaPlot_random_genes" 3>&1 1>&2 2>&3
        ) || SCRIPT1_metaPlot_random_genes="$SCRIPT1_metaPlot_random_genes"
        ;;

      "PCA")
        if whiptail --yesno "Perform PCA analysis?" 10 60; then
          SCRIPT1_pca="yes"
        else
          SCRIPT1_pca="no"
        fi
        ;;

      "QC plots")
        if whiptail --yesno "Generate sample-level QC plots (distributions, scatter, correlation, variance)?" 10 60; then
          SCRIPT1_QC_plots="yes"
        else
          SCRIPT1_QC_plots="no"
        fi
        ;;

      "Total methylation bar-plots")
        if whiptail --yesno "Produce total methylation bar-plots?" 10 60; then
          SCRIPT1_total_methylation="yes"
        else
          SCRIPT1_total_methylation="no"
        fi
        ;;

      "CX Chromosome Plot")
        if whiptail --yesno "Generate chromosome-wide CX plots?" 10 60; then
          SCRIPT1_CX_ChrPlot="yes"
        else
          SCRIPT1_CX_ChrPlot="no"
        fi
        ;;

      "TEs distance and size")
        if whiptail --yesno "Analyze TEs total methylation by size and distance from centromere?" 10 60; then
          SCRIPT1_TEs_distance_n_size="yes"
        else
          SCRIPT1_TEs_distance_n_size="no"
        fi
        ;;

      "Total methylation annotations")
        if whiptail --yesno "Analyze total methylation per genic annotations?" 10 60; then
          SCRIPT1_total_meth_ann="yes"
        else
          SCRIPT1_total_meth_ann="no"
        fi
        ;;

      "TF motifs")
        if whiptail --yesno "Analyze transcription factor motifs?" 10 60; then
          SCRIPT1_TF_motifs="yes"
        else
          SCRIPT1_TF_motifs="no"
        fi
        ;;

      "Functional groups")
        if whiptail --yesno "Analyze functional groups genes overlap DMRs?" 10 60; then
          SCRIPT1_func_groups="yes"
        else
          SCRIPT1_func_groups="no"
        fi
        ;;

      "Disable DMRs analysis")
        if whiptail --yesno "Disable DMRs analysis?" 10 60; then
          SCRIPT1_disable_DMRs="yes"
        else
          SCRIPT1_disable_DMRs="no"
        fi
        ;;

      "Strand-specific DMRs")
        if whiptail --yesno "Analyze strand-specific DMRs?" 10 60; then
          SCRIPT1_strand_DMRs="yes"
        else
          SCRIPT1_strand_DMRs="no"
        fi
        ;;

      "DMVs analysis")
        if whiptail --yesno "Analyze differentially methylated vallies (1kbp)?" 10 60; then
          SCRIPT1_DMVs="yes"
        else
          SCRIPT1_DMVs="no"
        fi
        ;;

      "dH analysis")
        if whiptail --yesno "Analyze dH (delta H)?" 10 60; then
          SCRIPT1_delta_H="yes"
        else
          SCRIPT1_delta_H="no"
        fi
        ;;

      "Annotation file (gtf/gff/csv)")
        SCRIPT1_annotation_file=$(whiptail --inputbox "Path to gene annotation file (leave bundled value for the selected reference bundle)" 10 78 "$SCRIPT1_annotation_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_annotation_file")
        ;;

      "Gene descriptions (txt/csv)")
        SCRIPT1_description_file=$(whiptail --inputbox "Path to gene-description table" 10 78 "$SCRIPT1_description_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_description_file")
        ;;

      "TE annotation file (txt)")
        SCRIPT1_TEs_file=$(whiptail --inputbox "Path to TE annotation file" 10 78 "$SCRIPT1_TEs_file" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_TEs_file")
        ;;

      "Reference bundle (yaml)")
        SCRIPT1_reference_bundle=$(whiptail --inputbox "Path to species/assembly reference bundle YAML" 10 78 "$SCRIPT1_reference_bundle" 3>&1 1>&2 2>&3 || echo "$SCRIPT1_reference_bundle")
        ;;

      *)
        # Unknown / spacer items (if you add headers later)
        ;;
    esac
  done
}


###################
# Gather run_bismark.sh
###################
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Bismark" ]]; then
    # Initialize parameters with defaults
    SCRIPT_BIS_genome="$SCRIPT_BIS_DEFAULT_genome"
    SCRIPT_BIS_ncores="$SCRIPT_BIS_DEFAULT_ncores"

    # Directly go to the parameters selection menu
    edit_script_bis_parameters || exit 1
fi

###################
# Gather Methylome.Plants.sh
###################
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.Plants" ]]; then
    # Initialize parameters with defaults
    SCRIPT1_minProportionDiff_CG="$SCRIPT1_DEFAULT_minProportionDiff_CG"
    SCRIPT1_minProportionDiff_CHG="$SCRIPT1_DEFAULT_minProportionDiff_CHG"
    SCRIPT1_minProportionDiff_CHH="$SCRIPT1_DEFAULT_minProportionDiff_CHH"
    SCRIPT1_binSize="$SCRIPT1_DEFAULT_binSize"
    SCRIPT1_minCytosinesCount="$SCRIPT1_DEFAULT_minCytosinesCount"
    SCRIPT1_minReadsPerCytosine="$SCRIPT1_DEFAULT_minReadsPerCytosine"
    SCRIPT1_pValueThreshold="$SCRIPT1_DEFAULT_pValueThreshold"
    SCRIPT1_n_cores="$SCRIPT1_DEFAULT_n_cores"
    SCRIPT1_pca="$SCRIPT1_DEFAULT_pca"
    SCRIPT1_QC_plots="$SCRIPT1_DEFAULT_QC_plots"
    SCRIPT1_total_methylation="$SCRIPT1_DEFAULT_total_methylation"
    SCRIPT1_CX_ChrPlot="$SCRIPT1_DEFAULT_CX_ChrPlot"
    SCRIPT1_TEs_distance_n_size="$SCRIPT1_DEFAULT_TEs_distance_n_size"
    SCRIPT1_total_meth_ann="$SCRIPT1_DEFAULT_total_meth_ann"
    SCRIPT1_TF_motifs="$SCRIPT1_DEFAULT_TF_motifs"
    SCRIPT1_func_groups="$SCRIPT1_DEFAULT_func_groups"
    SCRIPT1_GO_analysis="$SCRIPT1_DEFAULT_GO_analysis"
    SCRIPT1_KEGG_pathways="$SCRIPT1_DEFAULT_KEGG_pathways"
    SCRIPT1_file_type="$SCRIPT1_DEFAULT_file_type"
    SCRIPT1_img_type="$SCRIPT1_DEFAULT_img_type"
    SCRIPT1_annotation_file="$SCRIPT1_DEFAULT_annotation_file"
    SCRIPT1_description_file="$SCRIPT1_DEFAULT_description_file"
    SCRIPT1_TEs_file="$SCRIPT1_DEFAULT_TEs_file"
    SCRIPT1_reference_bundle="$SCRIPT1_DEFAULT_reference_bundle"
    SCRIPT1_disable_DMRs="$SCRIPT1_DEFAULT_disable_DMRs"
    SCRIPT1_strand_DMRs="$SCRIPT1_DEFAULT_strand_DMRs"
    SCRIPT1_DMVs="$SCRIPT1_DEFAULT_DMVs"
    SCRIPT1_delta_H="$SCRIPT1_DEFAULT_delta_H"
    SCRIPT1_TEs_metaplots="$SCRIPT1_DEFAULT_TEs_metaplots"
    SCRIPT1_Genes_metaplots="$SCRIPT1_DEFAULT_Genes_metaplots"
    SCRIPT1_Gene_features_metaplots="$SCRIPT1_DEFAULT_Gene_features_metaplots"
    SCRIPT1_bin_size_features="$SCRIPT1_DEFAULT_bin_size_features"
    SCRIPT1_metaPlot_random_genes="$SCRIPT1_DEFAULT_metaPlot_random_genes"

    # Directly go to the parameters selection menu
    edit_script1_parameters || exit 1
fi

###########
# RUN SCRIPTS
###########

# Construct a message listing the chosen scripts
chosen_message=""
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Bismark" ]]; then
    chosen_message+="'Bismark' "
fi
if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.Plants" ]]; then
    chosen_message+="'Methylome.Plants' "
fi

# Trim trailing space
chosen_message=$(echo "$chosen_message" | sed 's/[[:space:]]*$//')

# If somehow no scripts are chosen (shouldn't happen due to earlier checks), handle gracefully
if [ -z "$chosen_message" ]; then
    chosen_message="No scripts selected"
fi

# Display yes/no dialog
if (whiptail --title "All done!" --yesno "You have chosen to run: $chosen_message.\n\nWould you like to proceed?" 12 70); then

  cd "$Methylome_At_path"

  # Bismark pipeline for 'cx_report' files
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Bismark" ]]; then
    echo "Running run_bismark.sh..."
    SAMPLES_FILE_CX=$(bash "$SCRIPT_BIS_PATH" -s "$SAMPLES_FILE" -g "$SCRIPT_BIS_genome" -n "$SCRIPT_BIS_ncores" -o "${Methylome_At_path}/bismark_CX_reports" --cx --mat | tail -n1)
  else
    SAMPLES_FILE_CX="$SAMPLES_FILE"
  fi

  cd "$Methylome_At_path"

  # Methylome.Plants pipeline invocation
  if [[ " ${SELECTED_SCRIPTS[*]} " =~ "Methylome.Plants" ]]; then
    echo "Running Methylome.Plants.sh..."
    bash "$SCRIPT1_PATH" \
      --samples_file "$SAMPLES_FILE_CX" \
      --minProportionDiff_CG "$SCRIPT1_minProportionDiff_CG" \
      --minProportionDiff_CHG "$SCRIPT1_minProportionDiff_CHG" \
      --minProportionDiff_CHH "$SCRIPT1_minProportionDiff_CHH" \
      --binSize "$SCRIPT1_binSize" \
      --minCytosinesCount "$SCRIPT1_minCytosinesCount" \
      --minReadsPerCytosine "$SCRIPT1_minReadsPerCytosine" \
      --pValueThreshold "$SCRIPT1_pValueThreshold" \
      --n_cores "$SCRIPT1_n_cores" \
      $( [ "$SCRIPT1_pca" = "yes" ] && echo "--pca" ) \
      $( [ "$SCRIPT1_QC_plots" = "no" ] && echo "--QC_off" ) \
      $( [ "$SCRIPT1_total_methylation" = "yes" ] && echo "--total_methylation" ) \
      $( [ "$SCRIPT1_CX_ChrPlot" = "yes" ] && echo "--CX_ChrPlot" ) \
      $( [ "$SCRIPT1_TEs_distance_n_size" = "yes" ] && echo "--TEs_distance_n_size" ) \
      $( [ "$SCRIPT1_total_meth_ann" = "yes" ] && echo "--total_meth_ann" ) \
      $( [ "$SCRIPT1_TF_motifs" = "yes" ] && echo "--TF_motifs" ) \
      $( [ "$SCRIPT1_func_groups" = "yes" ] && echo "--func_groups" ) \
      $( [ "$SCRIPT1_GO_analysis" = "yes" ] && echo "--GO_analysis" ) \
      $( [ "$SCRIPT1_KEGG_pathways" = "yes" ] && echo "--KEGG_pathways" ) \
      --file_type "$SCRIPT1_file_type" \
      --image_type "$SCRIPT1_img_type" \
      --annotation_file "$SCRIPT1_annotation_file" \
      --description_file "$SCRIPT1_description_file" \
      --TEs_file "$SCRIPT1_TEs_file" \
      --reference_bundle "$SCRIPT1_reference_bundle" \
      $( [ "$SCRIPT1_TEs_metaplots" = "yes" ] && echo "--MP_TEs" ) \
      $( [ "$SCRIPT1_Genes_metaplots" = "yes" ] && echo "--MP_Genes" ) \
      $( [ "$SCRIPT1_Gene_features_metaplots" = "yes" ] && echo "--MP_Gene_features" ) \
      --MP_features_bin_size "$SCRIPT1_bin_size_features" \
      --metaPlot_random "$SCRIPT1_metaPlot_random_genes" \
      $( [ "$SCRIPT1_disable_DMRs" = "yes" ] && echo "--DMRs_off" ) \
      $( [ "$SCRIPT1_strand_DMRs" = "yes" ] && echo "--strand_DMRs" ) \
      $( [ "$SCRIPT1_DMVs" = "yes" ] && echo "--DMVs" ) \
      $( [ "$SCRIPT1_delta_H" = "yes" ] && echo "--dH" )
  fi
else

  echo "You chose not to run the scripts. Exiting."
  exit 0
fi
