#!/bin/bash

# Anaconda packages to install
packages=("r-curl" "r-rcurl" "zlib" "r-textshaping" "harfbuzz" "fribidi" "freetype" "libpng" "pkg-config" "libxml2" "r-xml" "bioconductor-rsamtools" "r-svglite") # "r-devtools"

#---------------------------------------------------------------------------#

# Default values for optional arguments
CHECK=false
PERMISSION=false

# Function to display help text
usage() {
    echo ""
    echo "Usage: $0 [options]"
    echo ""
    echo "Optional arguments:"
    echo "  --check         Check which R and Anaconda packages are installed"
    echo "  --permission    Ensure file permissions"
    echo ""
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --check) CHECK="true"; shift ;;
        --permission) PERMISSION="true"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
    shift
done

#---------------------------------------------------------------------------#

# Initialize conda
eval "$(conda shell.bash hook)"

#---------------------------------------------------------------------------#

if [ "$CHECK" == "true" ]; then
    # Activate Methylome.At environment
    if [ "$CONDA_DEFAULT_ENV" != "Methylome.At_env" ]; then
        eval "$(conda shell.bash hook)"
        conda activate Methylome.At_env
        if [ "$CONDA_DEFAULT_ENV" != "Methylome.At_env" ]; then
            echo "Error: Failed to activate the 'Methylome.At_env' Conda environment."
            echo "Please activate it manually using 'conda activate Methylome.At_env' and rerun the script."
            exit 1
        fi
    fi
    
    # Check installed R packages
    Rscript scripts/install_R_packages.R true
    
    # Check installed conda-forge packages
    echo -e "\n\n \tChecking installed Conda packages from conda-forge/bioconda:"
    for pkg in "r-base" "${packages[@]}"; do
        pkg_info=$(conda list "$pkg")
        channel=$(echo "$pkg_info" | awk '/^'"$pkg"'[[:space:]]/ {print $NF}')
        if [ "$channel" = "conda-forge" ] || [ "$channel" = "bioconda" ]; then
            echo -e "*\tinstalled $pkg: yes"
        else
            echo -e "*\tinstalled $pkg: no"
        fi
    done
    exit 0
fi

#---------------------------------------------------------------------------#

# Ensure unix line endings (can also try: sed -i 's/\r$//' "$sample_table")
dos2unix ./Methylome.A*.sh
dos2unix ./scripts/*.sh

# Permission
chmod +x ./Methylome.A*.sh
chmod +x ./scripts/*.sh

if [ "$PERMISSION" == "true" ]; then
    exit 0
fi

#---------------------------------------------------------------------------#

# Determine the directory of 'Methylome.At'
script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
cd "$script_dir"

# Generate log file with a timestamp
log_file="setup_env.log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "" >> "$log_file"

# Create and activate the conda environment
echo "Creating Conda environment..." >> "$log_file"
if conda create --name Methylome.At_env -c conda-forge -c bioconda r-base=4.4.3 "${packages[@]}" -y && conda activate Methylome.At_env; then
    echo "Conda environment created and activated successfully." >> "$log_file"
    echo "Conda environment name: $CONDA_DEFAULT_ENV" >> "$log_file"
else
    echo "Failed to create and activate Conda environment." >> "$log_file"
    exit 1
fi

echo "" >> "$log_file"

#---------------------------------------------------------------------------#

# Check if packages are installed from conda-forge
echo "Verifying that packages are installed from conda-forge..." >> "$log_file"
error=0
for pkg in "r-base" "${packages[@]}"; do
    # Get the package info
    pkg_info=$(conda list "$pkg")
    # Extract the channel from the package info
    channel=$(echo "$pkg_info" | awk '/^'"$pkg"'[[:space:]]/ {print $NF}')
    if [ "$channel" = "conda-forge" ] || [ "$channel" = "bioconda" ]; then
        echo "* installed $pkg: yes"
    else
        echo "* installed $pkg: no"
    fi
done

echo "" >> "$log_file"

#---------------------------------------------------------------------------#

# Run the R script to install additional R packages
echo ""
echo "Install R packages..." >> "$log_file"
Rscript scripts/install_R_packages.R false 2>> "$log_file"