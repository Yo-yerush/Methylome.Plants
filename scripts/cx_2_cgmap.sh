#!/bin/bash

# Check arguments
if [ $# -ne 2 ]; then
  echo "Usage: $0 <input.CX_report> <output.CGmap>"
  echo "
Column Description of 'CX' report:
1) Chromosome
2) Position
3) Strand
4) Count methylated
5) Count unmethylated
6) C-context (CG, CHG, CHH)
7) Trinucleotide context (CCG, CTG, CAG, etc.)

Column Description of 'CGmap' format:
1) Chromosome
2) Strand
3) Position
4) C-context (CG, CHG, CHH)
5) Dinucleotide context (CC, CT, CA, etc.)
6) Methylation level (3 decimal places)
7) Count methylated
8) Count total
  "
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Check if input file is gzipped and unzip if necessary
if [[ "${INPUT_FILE}" == *.gz ]]; then
    gunzip -c "${INPUT_FILE}" > tmp_cx2cgmap_gunzipped_file
    INPUT_FILE="tmp_cx2cgmap_gunzipped_file"
fi

awk 'BEGIN {
  FS="\t"; OFS="\t"
}
{
  chr = $1
  pos = $2
  strand = $3
  methyl = $4 + 0
  unmethyl = $5 + 0
  context = $6
  trinuc = $7
  dinuc = substr(trinuc, 1, 2)
  total = methyl + unmethyl
  if (total == 0) { level = "0.000" }
  else { level = sprintf("%.3f", methyl / total) }
  print chr, strand, pos, context, dinuc, level, methyl, total
}' "${INPUT_FILE}" > "${OUTPUT_FILE}"

# Remove temporary file if created
if [[ "${INPUT_FILE}" == "tmp_cx2cgmap_gunzipped_file" ]]; then
    rm tmp_cx2cgmap_gunzipped_file
fi

echo "Conversion complete."
echo "Output written to: ${OUTPUT_FILE}"
