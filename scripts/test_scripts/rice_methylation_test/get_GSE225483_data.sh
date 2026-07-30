#!/usr/bin/env bash

set -euo pipefail

cd /home/yoyerush/yo/methylome_plants_test_300726/rice_GSE225483
mkdir cov_files
mkdir genome_full
mkdir cx_report

##########################################################################

cd cov_files

# samples=(GSM8189821 GSM8189822 GSM8189823 GSM8189824 GSM8189825 GSM8189826 GSM8189827 GSM8189828)
# names=(osino80_rep1 osino80_rep2 osino80_rep3 osino80_rep4 WT_rep1 WT_rep2 WT_rep3 WT_rep4)
# 
# for i in "${!samples[@]}"; do
#   gsm="${samples[$i]}"
#   name="${names[$i]}"
#   wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm::-3}nnn/${gsm}/suppl/" -A "*.bismark%5Fbt2%5Fpe%2Ededuplicated%2Ebismark%2Ecov%2Egz" -O "${name}.bismark.cov.gz"
# done

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE225483&format=file" -O GSE225483_RAW.tar
# tar -xvf GSE225483_RAW.tar

cd ../

##########################################################################



##########################################################################

# coverage2cytosine \
#   --CX \
#   --gzip \
#   --genome_folder genome_full \
#   --output_dir cx_report \
#   -o WT_rep1_full_genome \
#   "cov_files/${COV_FILE}"
