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
REFERENCE_WIZARD_PATH="./scripts/reference_wizard/reference_wizard.R"
REFERENCE_CACHE_DIR="./reference_bundles/cache"
REFERENCE_GENERATED_DIR="./reference_bundles/generated"
R_SCRIPT_BIN="${R_SCRIPT_BIN:-Rscript}"

normalize_ui_path() {
  local value="$1"
  local first last

  # Trim pasted whitespace, then remove one matching pair of shell quotes.
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  if [ "${#value}" -ge 2 ]; then
    first="${value:0:1}"
    last="${value: -1}"
    if { [ "$first" = '"' ] && [ "$last" = '"' ]; } ||
       { [ "$first" = "'" ] && [ "$last" = "'" ]; }; then
      value="${value:1:${#value}-2}"
    fi
  fi

  printf '%s' "$value"
}

normalize_ui_path_var() {
  local variable_name="$1"
  printf -v "$variable_name" '%s' "$(normalize_ui_path "${!variable_name}")"
}

################################
# REFERENCE SETUP / WIZARD
################################

reference_optional_path() {
  local title="$1"
  local prompt="$2"
  local current="$3"
  whiptail --title "$title" --inputbox "$prompt\n\nLeave empty to skip this resource." \
    13 86 "$current" 3>&1 1>&2 2>&3
}

infer_primary_seqlevels() {
  local source_path="$1"
  local source_kind="$2"
  local values=""

  if [ -n "$source_path" ] && [ -f "$source_path" ]; then
    case "$source_kind" in
      fasta)
        if [[ "$source_path" =~ \.fai$ ]]; then
          values=$(awk 'BEGIN {ORS=","} {print $1}' "$source_path" | sed 's/,$//')
        else
          values=$(grep '^>' "$source_path" | sed 's/^>//; s/[[:space:]].*$//' | paste -sd, -)
        fi
        ;;
      sizes)
        values=$(awk -F'[\t,]' 'NR == 1 && tolower($1) == "seqname" {next} {print $1}' "$source_path" | paste -sd, -)
        ;;
      annotation)
        values=$(awk -F'[\t,]' '$1 !~ /^#/ && NF > 1 {print $1}' "$source_path" | awk '!seen[$0]++' | paste -sd, -)
        ;;
    esac
  fi

  printf '%s' "$values"
}

apply_reference_capabilities() {
  local bundle_path="$1"
  local capabilities

  capabilities=$("$R_SCRIPT_BIN" -e '
    args <- commandArgs(TRUE)
    source(file.path(args[2], "scripts/reference_bundle.R"))
    bundle <- read_reference_bundle(args[1])
    caps <- bundle_capabilities(bundle)
    cat(paste(names(caps)[caps], collapse = ","))
  ' "$bundle_path" "$Methylome_At_path" 2>/dev/null) || return 1

  reference_has_capability() {
    [[ ",$capabilities," == *",$1,"* ]]
  }

  if ! reference_has_capability transposable_elements; then
    SCRIPT1_TEs_distance_n_size="no"
    SCRIPT1_TEs_metaplots="no"
  fi
  reference_has_capability tfbs || SCRIPT1_TF_motifs="no"
  reference_has_capability functional_groups || SCRIPT1_func_groups="no"
  reference_has_capability go || SCRIPT1_GO_analysis="no"
  reference_has_capability kegg || SCRIPT1_KEGG_pathways="no"
  if ! reference_has_capability gene_annotation; then
    SCRIPT1_Genes_metaplots="no"
    SCRIPT1_Gene_features_metaplots="no"
  fi
  if ! reference_has_capability gene_annotation &&
     ! reference_has_capability transposable_elements; then
    SCRIPT1_total_meth_ann="no"
  fi
}

install_orgdb_package() {
  local package="$1"
  local log_file
  log_file=$(mktemp)
  whiptail --title "Installing OrgDb" --infobox \
    "Installing $package with BiocManager. This can take several minutes..." 8 76
  if "$R_SCRIPT_BIN" -e \
    'pkg <- commandArgs(TRUE)[1]; if (!requireNamespace("BiocManager", quietly=TRUE)) stop("BiocManager is not installed"); BiocManager::install(pkg, ask=FALSE, update=FALSE)' \
    "$package" >"$log_file" 2>&1; then
    rm -f "$log_file"
    return 0
  fi
  whiptail --title "OrgDb installation failed" --textbox "$log_file" 24 100
  rm -f "$log_file"
  return 1
}

select_orgdb_source() {
  WIZARD_ORGDB=""
  if ! whiptail --title "GO annotation" --yesno \
    "Configure GO enrichment with a Bioconductor OrgDb package?\n\nYou can instead add a custom gene-to-GO table in the optional-files step." \
    13 82; then
    return 0
  fi

  local search_term list_file error_file
  search_term=$(whiptail --title "Find an OrgDb" --inputbox \
    "Search by organism name or package name:" 10 78 "$WIZARD_ORGANISM" \
    3>&1 1>&2 2>&3) || return 0
  list_file=$(mktemp)
  error_file=$(mktemp)
  mkdir -p "$REFERENCE_CACHE_DIR"
  whiptail --title "OrgDb discovery" --infobox \
    "Searching installed and available Bioconductor OrgDb packages..." 8 76
  "$R_SCRIPT_BIN" "$REFERENCE_WIZARD_PATH" list-orgdb \
    --query "$search_term" --include-available \
    --cache "$REFERENCE_CACHE_DIR/orgdb_packages.tsv" \
    >"$list_file" 2>"$error_file" || true

  local -a choices
  choices=("none" "Do not use an OrgDb" ON "manual" "Enter an OrgDb package name manually" OFF)
  declare -A orgdb_status
  local package description status
  while IFS=$'\t' read -r package description status; do
    [ -z "$package" ] && continue
    choices+=("$package" "$description" OFF)
    orgdb_status["$package"]="$status"
  done <"$list_file"

  local selected
  selected=$(whiptail --title "GO / OrgDb organism" --radiolist \
    "Select the GO organism database. The wizard will test its gene identifiers against the GTF." \
    26 104 16 "${choices[@]}" 3>&1 1>&2 2>&3) || selected="none"

  if [ "$selected" = "manual" ]; then
    selected=$(whiptail --title "OrgDb package" --inputbox \
      "Enter the package name, for example org.At.tair.db:" 10 78 "" \
      3>&1 1>&2 2>&3) || selected=""
  fi

  if [ -n "$selected" ] && [ "$selected" != "none" ]; then
    if [ "${orgdb_status[$selected]}" = "available" ] || \
       ! "$R_SCRIPT_BIN" -e 'pkg <- commandArgs(TRUE)[1]; quit(status=if (requireNamespace(pkg, quietly=TRUE)) 0 else 1)' \
         "$selected" >/dev/null 2>&1; then
      if whiptail --title "Install OrgDb" --yesno \
        "$selected is not installed. Install it now with BiocManager?" 10 76; then
        install_orgdb_package "$selected" || selected=""
      else
        selected=""
      fi
    fi
  else
    selected=""
  fi

  WIZARD_ORGDB="$selected"
  rm -f "$list_file" "$error_file"
}

select_kegg_source() {
  WIZARD_KEGG=""
  if ! whiptail --title "KEGG annotation" --yesno \
    "Configure KEGG pathway analysis for this organism?" 10 72; then
    return 0
  fi

  local search_term list_file error_file
  search_term=$(whiptail --title "Find a KEGG organism" --inputbox \
    "Search the KEGG organism list by scientific name or code:" \
    10 82 "$WIZARD_ORGANISM" 3>&1 1>&2 2>&3) || return 0
  list_file=$(mktemp)
  error_file=$(mktemp)
  mkdir -p "$REFERENCE_CACHE_DIR"
  whiptail --title "KEGG discovery" --infobox \
    "Retrieving and filtering the KEGG organism list..." 8 72
  "$R_SCRIPT_BIN" "$REFERENCE_WIZARD_PATH" list-kegg \
    --query "$search_term" --cache "$REFERENCE_CACHE_DIR/kegg_organisms.tsv" \
    --limit 100 >"$list_file" 2>"$error_file" || true

  local -a choices
  choices=("none" "Do not configure KEGG" ON "manual" "Enter a KEGG organism code manually" OFF)
  local code description
  while IFS=$'\t' read -r code description; do
    [ -z "$code" ] && continue
    choices+=("$code" "$description" OFF)
  done <"$list_file"

  local selected
  selected=$(whiptail --title "KEGG organism" --radiolist \
    "Select the KEGG organism independently from the GO database." \
    26 100 16 "${choices[@]}" 3>&1 1>&2 2>&3) || selected="none"
  if [ "$selected" = "manual" ]; then
    selected=$(whiptail --title "KEGG organism code" --inputbox \
      "Enter the KEGG code, for example ath, osa, sly, or zma:" \
      10 76 "" 3>&1 1>&2 2>&3) || selected=""
  fi
  [ "$selected" = "none" ] && selected=""
  WIZARD_KEGG="$selected"
  rm -f "$list_file" "$error_file"
}

collect_optional_reference_files() {
  WIZARD_DESCRIPTIONS=""
  WIZARD_TE=""
  WIZARD_STABLE_GBM=""
  WIZARD_DYNAMIC_GBM=""
  WIZARD_CENTROMERES=""
  WIZARD_HETEROCHROMATIN=""
  WIZARD_TFBS=""
  WIZARD_GENE_TO_GO=""
  WIZARD_GENE_SETS=""
  WIZARD_KEGG_ID_MAP=""
  WIZARD_ALIASES=""

  local selected
  selected=$(whiptail --title "Optional reference resources" --checklist \
    "Select every resource you have. Missing resources simply disable their dependent analyses." \
    30 106 16 \
    "descriptions" "Gene-description table" OFF \
    "te" "TE annotation (BED/GFF3/GTF/CSV/TSV)" OFF \
    "stable_gbm" "Stable gbM gene-ID list" OFF \
    "dynamic_gbm" "Dynamic gbM gene-ID list" OFF \
    "centromeres" "Centromere BED" OFF \
    "heterochromatin" "Heterochromatin BED" OFF \
    "tfbs" "Transcription-factor binding sites BED" OFF \
    "gene_to_go" "Custom gene_id / go_id table" OFF \
    "gene_sets" "Functional gene sets (TSV or GMT)" OFF \
    "kegg_id_map" "Custom gene_id / kegg_id table" OFF \
    "aliases" "Sequence alias / canonical table" OFF \
    3>&1 1>&2 2>&3) || selected=""

  if [[ "$selected" == *'"descriptions"'* ]]; then
    WIZARD_DESCRIPTIONS=$(reference_optional_path "Gene descriptions" "Path to the description table" "") || WIZARD_DESCRIPTIONS=""
  fi
  if [[ "$selected" == *'"te"'* ]]; then
    WIZARD_TE=$(reference_optional_path "TE annotation" "Path to the TE annotation" "") || WIZARD_TE=""
  fi
  if [[ "$selected" == *'"stable_gbm"'* ]]; then
    WIZARD_STABLE_GBM=$(reference_optional_path "Stable gbM" "Path to the stable-gbM gene-ID list" "") || WIZARD_STABLE_GBM=""
  fi
  if [[ "$selected" == *'"dynamic_gbm"'* ]]; then
    WIZARD_DYNAMIC_GBM=$(reference_optional_path "Dynamic gbM" "Path to the dynamic-gbM gene-ID list" "") || WIZARD_DYNAMIC_GBM=""
  fi
  if [[ "$selected" == *'"centromeres"'* ]]; then
    WIZARD_CENTROMERES=$(reference_optional_path "Centromeres" "Path to the centromere BED file" "") || WIZARD_CENTROMERES=""
  fi
  if [[ "$selected" == *'"heterochromatin"'* ]]; then
    WIZARD_HETEROCHROMATIN=$(reference_optional_path "Heterochromatin" "Path to the heterochromatin BED file" "") || WIZARD_HETEROCHROMATIN=""
  fi
  if [[ "$selected" == *'"tfbs"'* ]]; then
    WIZARD_TFBS=$(reference_optional_path "TFBS" "Path to the TFBS BED file" "") || WIZARD_TFBS=""
  fi
  if [[ "$selected" == *'"gene_to_go"'* ]]; then
    WIZARD_GENE_TO_GO=$(reference_optional_path "Gene to GO" "Path to the table containing gene_id and go_id" "") || WIZARD_GENE_TO_GO=""
  fi
  if [[ "$selected" == *'"gene_sets"'* ]]; then
    WIZARD_GENE_SETS=$(reference_optional_path "Functional gene sets" "Path to the two-column TSV or GMT file" "") || WIZARD_GENE_SETS=""
  fi
  if [[ "$selected" == *'"kegg_id_map"'* ]]; then
    WIZARD_KEGG_ID_MAP=$(reference_optional_path "KEGG ID mapping" "Path to the table containing gene_id and kegg_id" "") || WIZARD_KEGG_ID_MAP=""
  fi
  if [[ "$selected" == *'"aliases"'* ]]; then
    WIZARD_ALIASES=$(reference_optional_path "Sequence aliases" "Path to the table containing alias and canonical" "") || WIZARD_ALIASES=""
  fi

  local path_variable
  for path_variable in \
    WIZARD_DESCRIPTIONS WIZARD_TE WIZARD_STABLE_GBM WIZARD_DYNAMIC_GBM \
    WIZARD_CENTROMERES WIZARD_HETEROCHROMATIN WIZARD_TFBS WIZARD_GENE_TO_GO \
    WIZARD_GENE_SETS WIZARD_KEGG_ID_MAP WIZARD_ALIASES; do
    normalize_ui_path_var "$path_variable"
  done
}

run_reference_wizard() {
  WIZARD_ORGANISM=$(whiptail --title "Plant organism" --inputbox \
    "Scientific/display name of the plant:" 10 80 "" 3>&1 1>&2 2>&3) || return 1
  [ -z "$WIZARD_ORGANISM" ] && return 1
  WIZARD_ASSEMBLY=$(whiptail --title "Genome assembly" --inputbox \
    "Assembly name or accession. This must match the annotation and methylation alignment:" \
    11 86 "" 3>&1 1>&2 2>&3) || return 1
  [ -z "$WIZARD_ASSEMBLY" ] && return 1
  WIZARD_ANNOTATION=$(whiptail --title "Gene annotation" --inputbox \
    "Path to the GTF, GFF3, GFF, or normalized CSV annotation:" \
    10 86 "" 3>&1 1>&2 2>&3) || return 1
  normalize_ui_path_var WIZARD_ANNOTATION
  [ -z "$WIZARD_ANNOTATION" ] && return 1

  select_orgdb_source
  select_kegg_source

  local default_length_mode="auto"
  local default_fasta=""
  if [ -n "${SCRIPT_BIS_genome:-}" ] && [ "$SCRIPT_BIS_genome" != "TAIR10" ]; then
    default_length_mode="fasta"
    default_fasta="$SCRIPT_BIS_genome"
  fi
  WIZARD_LENGTH_MODE=$(whiptail --title "Chromosome lengths" --radiolist \
    "Choose the source for exact chromosome lengths. Automatic uses CHRLENGTHS from the selected OrgDb when available." \
    18 100 5 \
    "auto" "Automatic: selected OrgDb CHRLENGTHS" $([ "$default_length_mode" = "auto" ] && echo ON || echo OFF) \
    "fasta" "Reference FASTA or FASTA index (.fai)" $([ "$default_length_mode" = "fasta" ] && echo ON || echo OFF) \
    "sizes" "Two-column chromosome-size table" OFF \
    "none" "No exact lengths; use observed analysis ranges" OFF \
    3>&1 1>&2 2>&3) || WIZARD_LENGTH_MODE="auto"
  WIZARD_FASTA=""
  WIZARD_CHROMOSOME_SIZES=""
  if [ "$WIZARD_LENGTH_MODE" = "fasta" ]; then
    WIZARD_FASTA=$(whiptail --title "Reference FASTA" --inputbox \
      "Path to the FASTA or .fai file:" 10 86 "$default_fasta" \
      3>&1 1>&2 2>&3) || WIZARD_FASTA=""
    normalize_ui_path_var WIZARD_FASTA
  elif [ "$WIZARD_LENGTH_MODE" = "sizes" ]; then
    WIZARD_CHROMOSOME_SIZES=$(whiptail --title "Chromosome sizes" --inputbox \
      "Path to a TSV containing seqname and length columns:" 10 86 "" \
      3>&1 1>&2 2>&3) || WIZARD_CHROMOSOME_SIZES=""
    normalize_ui_path_var WIZARD_CHROMOSOME_SIZES
  fi

  local primary_default=""
  if [ -n "$WIZARD_FASTA" ]; then
    primary_default=$(infer_primary_seqlevels "$WIZARD_FASTA" "fasta")
  elif [ -n "$WIZARD_CHROMOSOME_SIZES" ]; then
    primary_default=$(infer_primary_seqlevels "$WIZARD_CHROMOSOME_SIZES" "sizes")
  else
    primary_default=$(infer_primary_seqlevels "$WIZARD_ANNOTATION" "annotation")
  fi
  WIZARD_PRIMARY_SEQLEVELS=$(whiptail --title "Primary chromosomes" --inputbox \
    "Comma-separated primary sequence names to analyze.\n\nThe default is inferred from the length source when possible. Remove unplaced/scaffold sequences such as ChrUn if they are not in the FASTA." \
    14 100 "$primary_default" 3>&1 1>&2 2>&3) || WIZARD_PRIMARY_SEQLEVELS="$primary_default"

  WIZARD_CHLOROPLAST=""
  WIZARD_MITOCHONDRIAL=""
  if whiptail --title "Organelle sequences" --yesno \
    "Configure chloroplast or mitochondrial sequence names for QC and ordering?" 10 78; then
    WIZARD_CHLOROPLAST=$(reference_optional_path "Chloroplast" "Canonical chloroplast sequence name(s), comma-separated" "") || WIZARD_CHLOROPLAST=""
    WIZARD_MITOCHONDRIAL=$(reference_optional_path "Mitochondria" "Canonical mitochondrial sequence name(s), comma-separated" "") || WIZARD_MITOCHONDRIAL=""
  fi

  collect_optional_reference_files

  local organism_slug assembly_slug default_output output_path
  organism_slug=$(printf '%s' "$WIZARD_ORGANISM" | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z0-9]/_/g; s/__*/_/g; s/^_//; s/_$//')
  assembly_slug=$(printf '%s' "$WIZARD_ASSEMBLY" | sed 's/[^A-Za-z0-9._-]/_/g')
  mkdir -p "$REFERENCE_GENERATED_DIR"
  default_output="$REFERENCE_GENERATED_DIR/${organism_slug}_${assembly_slug}.yaml"
  output_path=$(whiptail --title "Save generated reference" --inputbox \
    "The wizard configuration will be saved and can be reused later:" \
    11 90 "$default_output" 3>&1 1>&2 2>&3) || return 1
  output_path=$(normalize_ui_path "$output_path")
  if [ -f "$output_path" ]; then
    whiptail --title "Replace generated reference?" --yesno \
      "$output_path already exists. Replace this generated configuration?" 10 92 || return 1
  fi

  local -a args
  args=(generate --organism "$WIZARD_ORGANISM" --assembly "$WIZARD_ASSEMBLY" \
    --annotation "$WIZARD_ANNOTATION" --output "$output_path")
  [ -n "$WIZARD_ORGDB" ] && args+=(--orgdb-package "$WIZARD_ORGDB" --orgdb-keytype auto)
  [ -n "$WIZARD_KEGG" ] && args+=(--kegg-organism "$WIZARD_KEGG")
  [ -n "$WIZARD_FASTA" ] && args+=(--fasta "$WIZARD_FASTA")
  [ -n "$WIZARD_CHROMOSOME_SIZES" ] && args+=(--chromosome-sizes "$WIZARD_CHROMOSOME_SIZES")
  [ -n "$WIZARD_PRIMARY_SEQLEVELS" ] && args+=(--primary-seqlevels "$WIZARD_PRIMARY_SEQLEVELS")
  [ "$WIZARD_LENGTH_MODE" = "none" ] && args+=(--disable-orgdb-lengths)
  [ -n "$WIZARD_DESCRIPTIONS" ] && args+=(--descriptions "$WIZARD_DESCRIPTIONS")
  [ -n "$WIZARD_TE" ] && args+=(--te "$WIZARD_TE")
  [ -n "$WIZARD_STABLE_GBM" ] && args+=(--stable-gbm "$WIZARD_STABLE_GBM")
  [ -n "$WIZARD_DYNAMIC_GBM" ] && args+=(--dynamic-gbm "$WIZARD_DYNAMIC_GBM")
  [ -n "$WIZARD_CENTROMERES" ] && args+=(--centromeres "$WIZARD_CENTROMERES")
  [ -n "$WIZARD_HETEROCHROMATIN" ] && args+=(--heterochromatin "$WIZARD_HETEROCHROMATIN")
  [ -n "$WIZARD_TFBS" ] && args+=(--tfbs "$WIZARD_TFBS")
  [ -n "$WIZARD_GENE_TO_GO" ] && args+=(--gene-to-go "$WIZARD_GENE_TO_GO")
  [ -n "$WIZARD_GENE_SETS" ] && args+=(--gene-sets "$WIZARD_GENE_SETS")
  [ -n "$WIZARD_KEGG_ID_MAP" ] && args+=(--kegg-id-map "$WIZARD_KEGG_ID_MAP")
  [ -n "$WIZARD_ALIASES" ] && args+=(--seqname-aliases "$WIZARD_ALIASES")
  [ -n "$WIZARD_CHLOROPLAST" ] && args+=(--chloroplast-seqlevels "$WIZARD_CHLOROPLAST")
  [ -n "$WIZARD_MITOCHONDRIAL" ] && args+=(--mitochondrial-seqlevels "$WIZARD_MITOCHONDRIAL")

  local summary_file error_file
  summary_file=$(mktemp)
  error_file=$(mktemp)
  whiptail --title "Reference validation" --infobox \
    "Inspecting annotation fields, identifiers, chromosome names, and optional resources..." 8 86
  if ! "$R_SCRIPT_BIN" "$REFERENCE_WIZARD_PATH" "${args[@]}" \
    >"$summary_file" 2>"$error_file"; then
    whiptail --title "Reference setup failed" --textbox "$error_file" 28 108
    rm -f "$summary_file" "$error_file"
    return 1
  fi

  SCRIPT1_reference_bundle=$(readlink -f "$output_path")
  SCRIPT1_reference_mode="wizard"
  SCRIPT1_annotation_file=""
  SCRIPT1_description_file=""
  SCRIPT1_TEs_file=""
  [ -n "$WIZARD_ORGDB$WIZARD_GENE_TO_GO" ] && SCRIPT1_GO_analysis="yes" || SCRIPT1_GO_analysis="no"
  [ -n "$WIZARD_KEGG" ] && SCRIPT1_KEGG_pathways="yes" || SCRIPT1_KEGG_pathways="no"
  [ -n "$WIZARD_TE" ] || { SCRIPT1_TEs_distance_n_size="no"; SCRIPT1_TEs_metaplots="no"; }
  [ -n "$WIZARD_TFBS" ] || SCRIPT1_TF_motifs="no"
  [ -n "$WIZARD_GENE_SETS" ] || SCRIPT1_func_groups="no"
  apply_reference_capabilities "$SCRIPT1_reference_bundle" || true

  local summary_text
  summary_text=$(sed $'s/\t/: /' "$summary_file")
  whiptail --title "Reference ready" --msgbox \
    "$summary_text\n\nThe generated bundle will be saved with the analysis and can be reused." \
    26 104
  rm -f "$summary_file" "$error_file"
}

configure_reference() {
  local current_mode="${SCRIPT1_reference_mode:-wizard}"
  local selected
  selected=$(whiptail --title "Plant reference setup" --radiolist \
    "The wizard is recommended. YAML bundles remain available for reproducible or advanced setups." \
    17 96 4 \
    "wizard" "Create a reference from GTF/GFF and optional files" $([ "$current_mode" = "wizard" ] && echo ON || echo OFF) \
    "bundle" "Load an existing reference-bundle YAML" $([ "$current_mode" = "bundle" ] && echo ON || echo OFF) \
    "tair10" "Use the bundled Arabidopsis TAIR10 reference" $([ "$current_mode" = "tair10" ] && echo ON || echo OFF) \
    3>&1 1>&2 2>&3) || return 1

  case "$selected" in
    wizard)
      run_reference_wizard
      ;;
    bundle)
      local bundle_path
      bundle_path=$(whiptail --title "Existing reference bundle" --inputbox \
        "Path to the species/assembly YAML bundle:" 10 88 "$SCRIPT1_reference_bundle" \
        3>&1 1>&2 2>&3) || return 1
      bundle_path=$(normalize_ui_path "$bundle_path")
      [ ! -f "$bundle_path" ] && {
        whiptail --title "Reference bundle" --msgbox "Bundle not found: $bundle_path" 10 88
        return 1
      }
      SCRIPT1_reference_bundle=$(readlink -f "$bundle_path")
      SCRIPT1_reference_mode="bundle"
      SCRIPT1_annotation_file=""
      SCRIPT1_description_file=""
      SCRIPT1_TEs_file=""
      if ! apply_reference_capabilities "$SCRIPT1_reference_bundle"; then
        whiptail --title "Reference bundle" --msgbox \
          "Could not inspect capabilities in: $SCRIPT1_reference_bundle" 10 92
        return 1
      fi
      ;;
    tair10)
      SCRIPT1_reference_bundle=$(readlink -f "$SCRIPT1_DEFAULT_reference_bundle")
      SCRIPT1_reference_mode="tair10"
      SCRIPT1_annotation_file=""
      SCRIPT1_description_file=""
      SCRIPT1_TEs_file=""
      apply_reference_capabilities "$SCRIPT1_reference_bundle" || true
      ;;
  esac
}

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
if ! SAMPLES_FILE=$(whiptail --title "samples_file" --inputbox \
  "Enter the path to the samples table file:" \
  10 70 \
  "" \
  3>&1 1>&2 2>&3); then
  echo "samples_file is required. Exiting."
  exit 1
fi
SAMPLES_FILE=$(normalize_ui_path "$SAMPLES_FILE")
if [ -z "$SAMPLES_FILE" ]; then
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
      SCRIPT_BIS_DEFAULT_genome=$(normalize_ui_path "$SCRIPT_BIS_DEFAULT_genome")
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
            --menu "Select a parameter to change or proceed with current settings." 40 108 33 \
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
            "Reference setup"          "$(fmt "$SCRIPT1_reference_mode" "$SCRIPT1_reference_bundle")" \
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

      "Reference setup")
        configure_reference
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
    SCRIPT1_reference_mode="wizard"
    SCRIPT1_disable_DMRs="$SCRIPT1_DEFAULT_disable_DMRs"
    SCRIPT1_strand_DMRs="$SCRIPT1_DEFAULT_strand_DMRs"
    SCRIPT1_DMVs="$SCRIPT1_DEFAULT_DMVs"
    SCRIPT1_delta_H="$SCRIPT1_DEFAULT_delta_H"
    SCRIPT1_TEs_metaplots="$SCRIPT1_DEFAULT_TEs_metaplots"
    SCRIPT1_Genes_metaplots="$SCRIPT1_DEFAULT_Genes_metaplots"
    SCRIPT1_Gene_features_metaplots="$SCRIPT1_DEFAULT_Gene_features_metaplots"
    SCRIPT1_bin_size_features="$SCRIPT1_DEFAULT_bin_size_features"
    SCRIPT1_metaPlot_random_genes="$SCRIPT1_DEFAULT_metaPlot_random_genes"

    # Configure the plant reference before selecting dependent analyses.
    configure_reference || exit 1

    # Continue to the analysis parameter menu.
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
