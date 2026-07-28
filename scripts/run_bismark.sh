#!/bin/bash

usage_yo="
###############################################################################
YO - 260525
Bismark WGBS pipeline

------------------------------------------------------------

Usage:
------
run_bismark.sh [-s <required>] [-g TAIR10] [options]

Options:
--------
-s, --samples   Tab-delimited two-column file: sample-name <TAB> fastq-path
-g, --genome    FASTA of the reference genome [default: TAIR10]
-o, --outdir    Output directory [default: ./bismark_results]
-n, --ncores    Number of cores (max). multiples of 4 recommended [default: 8]
-m, --mem       Buffer size for 'bismark_methylation_extractor' [default: 8G]
--cx            Produce and keep only '_CX_report.txt.gz' file
--mat           Produce samples table (.txt) for 'Methylome.At' pipeline
--indx          Keep the genome index directory (applies only if --cx is on)
--sort          Sort & index BAM files (applies only if --cx is off)
--strand        Keep top/bottom strand (OT/OB) files [remove in default]
--um            Produce and keep only unmapped files (as FASTQ)
--help

------------------------------------------------------------

Example:
--------

Create a sample table file (example):
-------------------------------------
mt_1    PATH/TO/FILE/mt1_R1.fastq
mt_1    PATH/TO/FILE/mt1_R2.fastq
mt_2    PATH/TO/FILE/mt2_R1.fastq
mt_2    PATH/TO/FILE/mt2_R2.fastq
wt_1    PATH/TO/FILE/wt1_R1.fastq
wt_1    PATH/TO/FILE/wt1_R2.fastq
wt_2    PATH/TO/FILE/wt2_R1.fastq
wt_2    PATH/TO/FILE/wt2_R2.fastq

Run:
----
$ ./run_bismark.sh -s samples_table.txt -g TAIR10 -n 30 --cx --mat
###############################################################################
"

####################
### default values
sample_table=
genome_file_name="TAIR10"
output_path="./bismark_results"
output_suffix="wgbs_bismark_$(date +%d%m%y)"
n_cores=8
buffer_size=8G
keep_cx=false
methAt_samples=false
keep_indx=false
sort_bam=false
keep_strand=false
keep_unmapped=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -s | --samples)
            sample_table=$2
            shift 2
        ;;
        -g | --genome)
            genome_file_name=$2
            shift 2
        ;;
        -o | --outdir)
            output_path=$2
            shift 2
        ;;
        -n | --ncores)
            n_cores=$2
            shift 2
        ;;
        -m | --mem)
            buffer_size=$2
            shift 2
        ;;
        --cx)
            keep_cx=true
            shift
        ;;
        --mat)
            methAt_samples=true
            shift
        ;;
        --indx)
            keep_indx=true
            shift
        ;;
        --sort)
            sort_bam=true
            shift
        ;;
        --strand)
            keep_strand=true
            shift
        ;;
        --um)
            keep_unmapped=true
            shift
        ;;
        -h | --help)
            echo "$usage_yo"
            exit 0
        ;;
        *)
            echo "Unknown option: $1"
            exit 1
        ;;
    esac
done

####################

# check if 'sample_table' file exists
if [[ ! -f "$sample_table" ]]; then
    echo "Error: Sample table file '$sample_table' does not exist."
    exit 1
fi

# ensure sample table has unix line endings (can also try: sed -i 's/\r$//' "$sample_table")
dos2unix "$sample_table" 2>/dev/null

# check if 'genome_file_name' file exists
if [[ ! -f "$genome_file_name" && "$genome_file_name" != "TAIR10" ]]; then
    echo "Error: Genome file '$genome_file_name' does not exist."
    exit 1
fi

# check for duplicate arguments
if [[ "$sort_bam" == "true" && "$keep_cx" == "true" ]]; then
    echo "Error: Can't use both '--cx' and '--sort' arguments"
    exit 1
fi

# n-cores for bismark
if [ $n_cores -gt 2 ]; then
    n_cores_2=$((n_cores / 4))
else
    n_cores_2=1
fi

####################
### read sample names and fastq file paths as an array
mapfile -t sample_name < <(awk '!seen[$1]++ {print $1}' "$sample_table")
mapfile -t R1_fastq_path < <(awk '$2 ~ /(_R1_|_R1\.fq|_R1\.fastq|_1\.fq|_1\.fastq)/ {print $2}' "$sample_table")
mapfile -t R2_fastq_path < <(awk '$2 ~ /(_R2_|_R2\.fq|_R2\.fastq|_2\.fq|_2\.fastq)/ {print $2}' "$sample_table")

paired_end_sequence=true
((${#R2_fastq_path[@]} == 0)) && paired_end_sequence=false
if [[ "$paired_end_sequence" == "false" ]]; then
    # For single-end, use all rows in the samples table as read file paths
    mapfile -t R1_fastq_path < <(awk '{print $2}' "$sample_table")
fi

####################

ori_path=$(pwd)
mkdir -p $output_path
cd $output_path
output_path=$(pwd)

### tmp file for analysis
mkdir -p $output_path/tmp
cd $output_path/tmp

### Generate log file with a timestamp
log_file="../${output_suffix}.log"
echo "**  $(date +"%d-%m-%y %H:%M")" > "$log_file"
echo "**  samples: ${sample_name[@]}" >> "$log_file"
#echo "paired-end sequence: $paired_end_sequence" >> "$log_file"
if [[ "$keep_cx" == "true" ]]; then
    echo "**  keep just 'CX_report' file" >> "$log_file"
fi
if [[ "$sort_bam" == "true" ]]; then
    echo "**  sort bam file" >> "$log_file"
fi
if [[ "$keep_strand" == "true" ]]; then
    echo "**  keep also top/bottom strand (OT/OB) files" >> "$log_file"
fi
if [[ "$keep_unmapped" == "true" ]]; then
    echo "**  keep just the unmapped files ('.fastq')" >> "$log_file"
fi

echo "" >> "$log_file"

####################
### index the genom
mkdir -p $output_path/genome_indx

# Get the genome file from TAIR10 or from file path
if [[ "$genome_file_name" == "TAIR10" ]]; then
    echo "download TAIR10 FASTA" >> "$log_file"
    echo -e "\n*************\n\ndownload the default TAIR10 reference FASTA file\n\n*************\n\n"
    wget -O "${output_path}/genome_indx/TAIR10_chr_all.fa.gz" "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz"
    genome_b_name="TAIR10_chr_all.fa.gz"
else
    cp $genome_file_name $output_path/genome_indx
    # Rename genome file if needed
    genome_b_name=$(basename "$genome_file_name")
    genome_new_path=$output_path/genome_indx/$genome_b_name
    if [[ "$genome_b_name" == *.fas ]]; then
        mv "$genome_new_path" "${genome_new_path%.fas}.fa"
        echo "* rename genome file: '$(basename ${genome_new_path%.fas}.fa)'" >> "$log_file"
    elif [[ "$genome_new_path" == *.fas.gz ]]; then
        mv "$genome_new_path" "${genome_new_path%.fas.gz}.fa.gz"
        echo "* rename genome file: '$(basename ${genome_new_path%.fas.gz}.fa.gz)'" >> "$log_file"
    fi
fi

echo "indexing genome file: '$genome_b_name'" >> "$log_file"
bismark_genome_preparation $output_path/genome_indx


####################
### create samples table file for Methylome.At
if [[ "$methAt_samples" == "true" ]]; then
    samples_table_tmp="${output_path}/tmp/S_T_$(date +"%y%m%d%H%M%S").tmp"
    > "$samples_table_tmp"
fi

echo "" >> "$log_file"
echo "-----------------------------------" >> "$log_file"

####################
# main loop
for ((u = 0; u < ${#sample_name[@]}; u++)); do
    echo "**  $(date +"%d-%m-%y %H:%M")" >> "$log_file"
    echo "" >> "$log_file"

    i="${sample_name[$u]}"
    R1_i="${R1_fastq_path[$u]}"
    R2_i="${R2_fastq_path[$u]}"

    echo "Processing sample: $i" >> "$log_file"
    mkdir -p "$output_path/$i"

    # # # # # # # # # # # #
    ### mapping to genome
    if [[ "$paired_end_sequence" == "false" ]]; then
        Rs_type="se"
        echo "mapping to genome for single-end sequence:" >> "$log_file"
        echo "read1 file: $R1_i" >> "$log_file"
        if [[ "$keep_unmapped" == "true" ]]; then
            bismark --bowtie2 --parallel "$n_cores_2" $output_path/genome_indx "$R1_i" -o $output_path/"$i" --prefix "$i" --unmapped
        else
            bismark --bowtie2 --parallel "$n_cores_2" $output_path/genome_indx "$R1_i" -o $output_path/"$i" --prefix "$i"
        fi
    else
        Rs_type="pe"
        echo "mapping to genome for paired-end sequence:" >> "$log_file"
        echo "* read1 file: '$(basename "$R1_i")'" >> "$log_file"
        echo "* read2 file: '$(basename "$R2_i")'" >> "$log_file"
        if [[ "$keep_unmapped" == "true" ]]; then
            bismark --bowtie2 --parallel "$n_cores_2" $output_path/genome_indx -1 "$R1_i" -2 "$R2_i" -o $output_path/"$i" --prefix "$i" --unmapped
        else
            bismark --bowtie2 --parallel "$n_cores_2" $output_path/genome_indx -1 "$R1_i" -2 "$R2_i" -o $output_path/"$i" --prefix "$i"
        fi
    fi
    if [[ "$keep_unmapped" == "true" ]]; then
        rm $output_path/"$i"/"$i"*.bam
        rm $output_path/"$i"/"$i"*_report.txt
        if [[ "$paired_end_sequence" == "true" ]]; then
            mv $output_path/"$i"/*unmapped_reads_1.fq.gz $output_path/"$i"/"$i"_unmapped_reads_1.fq.gz # rename
            mv $output_path/"$i"/*unmapped_reads_2.fq.gz $output_path/"$i"/"$i"_unmapped_reads_2.fq.gz # rename
            zcat $output_path/"$i"/"$i"_unmapped_reads_1.fq.gz $output_path/"$i"/"$i"_unmapped_reads_2.fq.gz > $output_path/"$i"/"$i"_unmapped_paired_reads.fq # concatenate
            gzip $output_path/"$i"/"$i"_unmapped_paired_reads.fq
        else
            mv $output_path/"$i"/*unmapped_reads.fq.gz $output_path/"$i"/"$i"_unmapped_reads.fq.gz # rename
        fi
    else
        mv $output_path/"$i"/"$i"*.bam $output_path/"$i"/"$i"_bismark_"$Rs_type".bam # rename
        mv $output_path/"$i"/"$i"*_report.txt $output_path/"$i"/"$i"_bismark_"$Rs_type"_report.txt # rename
    fi

    if [[ "$keep_unmapped" == "false" ]]; then
        # # # # # # # # # # # #
        ### methylation calling
        echo "" >> "$log_file"
        echo "methylation calling..." >> "$log_file"
        mkdir -p $output_path/"$i"/methylation_extractor

        if [[ "$keep_cx" == "true" ]]; then
            # run 'methylation_extractor' and keep 'CX_report' file only
            bismark_methylation_extractor --CX --cytosine_report --parallel "$n_cores_2" --buffer_size "$buffer_size" --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_bismark_"$Rs_type".bam

            gzip $output_path/"$i"/methylation_extractor/*.CX_report.txt
            mv  $output_path/"$i"/methylation_extractor/*.CX_report.txt.gz $output_path
            rm -r -- "$output_path/$i"

        else
            # run 'methylation_extractor' and keep files *without 'CX_report' file*
            bismark_methylation_extractor --CX --bedGraph --cytosine_report --parallel "$n_cores_2" --buffer_size "$buffer_size" --genome_folder $output_path/genome_indx -o $output_path/"$i"/methylation_extractor $output_path/"$i"/"$i"_bismark_"$Rs_type".bam

            gzip $output_path/"$i"/methylation_extractor/*.CX_report.txt

            if [[ "$keep_strand" == "false" ]]; then
                rm -v -- "$output_path/$i"/methylation_extractor/{CHG,CHH,CpG}_{CTOB,CTOT,OB,OT}_*
            fi
        fi

        # # # # # # # # # # # #
        # sort bam files (can use in IGV software to watch the reads over the genome)
        if [[ "$sort_bam" == "true" && "$keep_cx" == "false" ]]; then
            samtools sort $output_path/"$i"/"$i"_bismark_"$Rs_type".bam -o $output_path/"$i"/"$i"_sorted.bam
            samtools index $output_path/"$i"/"$i"_sorted.bam
        fi

        # # # # # # # # # # # #
        # samples table for Methylome.At
        if [[ "$methAt_samples" == "true" ]]; then
            i_unique=$(printf '%s\n' "$i" | sed 's/[._][0-9]*$//')
            if [[ "$keep_cx" == "true" ]]; then
                echo -e "$i_unique"\\t"$output_path"/"$i"_bismark_"$Rs_type".CX_report.txt.gz >> "$samples_table_tmp"
            else
                echo -e "$i_unique"\\t"$output_path"/"$i"/methylation_extractor/"$i"_bismark_"$Rs_type".CX_report.txt.gz >> $samples_table_tmp
            fi
        fi
    fi
    echo "" >> "$log_file"
    echo "-----------------------------------" >> "$log_file"
done

if [[ "$methAt_samples" == "true" ]]; then
    samples_var=$(cut -f1 "$samples_table_tmp" | uniq)
    control_s=$(echo "$samples_var" | head -n1)
    treatment_s=$(echo "$samples_var" | head -n2 | tail -n1)
    st_out_path=$(readlink -f "${output_path}/../samples_CX_table")
    st_out_name="${st_out_path}/samples_file_${treatment_s}_vs_${control_s}.txt"
    mkdir -p "$st_out_path"
    cat "$samples_table_tmp" > "$st_out_name"
    rm "$samples_table_tmp"
    echo "$st_out_name" # print the table path to the next script
fi

if [[ "$keep_cx" == "true" && "$keep_indx" == "false" ]]; then
    rm -r -- $output_path/genome_indx
fi

echo "**  $(date +"%d-%m-%y %H:%M")" >> "$log_file"
cd $ori_path
rm -r $output_path/tmp