set -euo pipefail

usage() {
cat <<'EOF'
-----------------------------------------------------------
./bedmethyl_2_cx.sh  - convert DeepSignal-plant bedMethyl to CX_report

-----------------------------------------------------------
USAGE EXAMPLES
  Without trinucleotide column:
  ./bedmethyl_2_cx.sh -i deepsig_output.bed -o test

  With trinucleotide column:
  ./bedmethyl_2_cx.sh -i deepsig_output.bed -t /path/to/genome_dir/ -o test

-----------------------------------------------------------
OPTIONS
  -i  <bed>         input bedMethyl file (required)
  -o  <prefix>      output prefix (default: input filename 'name.bed')
  -t  <genome_dir>  add trinucleotide column (if needed). insert genome directory path.
                    note: installed Bismark spftware and a genome file are required.

  -h | --help

-----------------------------------------------------------
EOF
}

# parse arguments
bed=""; prefix=""; genome_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i)  bed="$2";          shift 2 ;;
    -o)  prefix="$2";       shift 2 ;;
    -t)          
        # make sure the next token exists
        if [[ $# -lt 2 || "$2" == -* ]]; then
            echo "ERROR: -t option needs a genome directory path" >&2
            exit 1
        fi
        genome_dir="$2";    shift 2 ;;
    -h|--help)  usage; exit 0 ;;
    *)   echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -z $bed ]] && { echo "ERROR: -i <bed> is required"; exit 1; }
[[ -z $prefix ]] && prefix="${bed%.bed}"


CX="$prefix.CX_report.txt"

#############################

if [[ -z $genome_dir ]]; then
    awk 'BEGIN{OFS="\t"}{
            cov  = $10;                             # column 10 = total coverage
            meth = int(($11/100)*cov + 0.5);        # column 11 = % methylated
            unm  = cov - meth;                      # unmethylated reads
            pos  = $2 + 1;                          # 1-based coordinate
            print $1, pos, $6, meth, unm, $4, "."   # chr pos strand meth unm context.
         }'  "$bed"  >  "$CX"
    echo "* Done (without trinucleotide contexts):  $CX"

else
    COV="tmp_coverage_$prefix.txt"

    echo "* Converting $bed to Bismark coverage file..."
    awk '
    BEGIN { OFS = "\t" }
    {
        cov  = $10;                     # total reads
        pct  = $11;                     # % methylated (float)
        meth = int(pct/100*cov + 0.5);  # round to nearest
        unm  = cov - meth;
        pos  = $2 + 1;                  # 1-based for coverage format
        print $1, pos, pos, int(pct + 0.5), meth, unm
    }
    ' "$bed" > "$COV"

    trap 'echo "
    Remove coverage temporary file..."; rm -f "$COV"' EXIT

    echo "* Converting coverage to CX_report file..."
    echo "* Running 'coverage2cytosine' to produce CX_report..."
    coverage2cytosine --CX --gzip --genome_folder "$genome_dir" -o "$CX" "$COV"
    echo "* Remove Bismark coverage file..."
    rm "$COV"
    echo "* Done (with trinucleotide contexts):  $CX"
fi