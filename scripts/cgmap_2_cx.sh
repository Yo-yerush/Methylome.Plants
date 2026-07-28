#!/bin/bash

# Check arguments
if [ $# -ne 2 ]; then
  echo "Usage: $0 <input.cgmap> <output.CX_report>"
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
6) Methylation level
7) Count methylated
8) Count total
  "
  exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Check if input file is gzipped and unzip if necessary
if [[ "${INPUT_FILE}" == *.gz ]]; then
    gunzip -c "${INPUT_FILE}" > tmp_cgmap2cx_gunzipeed_file
    INPUT_FILE="tmp_cgmap2cx_gunzipeed_file"
fi

awk 'BEGIN {
  FS="\t"; OFS="\t"
}
{
  if ($2 == "C") { strand = "+" } else { strand = "-" }
  unmethyl = $8 - $7
  print $1, $3, strand, $7, unmethyl, $4, $5
}' "${INPUT_FILE}" > "${OUTPUT_FILE}"

# Check if input file is gzipped remove 'tmp' file
if [[ "${INPUT_FILE}" == tmp_cgmap2cx_gunzipeed_file ]]; then
    rm tmp_cgmap2cx_gunzipeed_file
fi

echo "Conversion complete."
echo "Output written to: ${OUTPUT_FILE}"
