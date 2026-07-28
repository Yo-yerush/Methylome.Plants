# Methylome.At

Methylome.At is a comprehensive, R-based pipeline for *Arabidopsis thaliana* that processes post-alignment **WGBS** or **Nanopore** sequencing data for CG, CHG and CHH DNA methylation contexts, identifies differentially methylated regions (DMRs, using [DMRcaller](https://github.com/nrzabet/DMRcaller) package) to replicates/single samples data, integrates multiple genomic resources for functional interpretation, and generates extensive visualizations and annotations to advance understanding of plant epigenetic regulation.

---

```mermaid
%%{init: {'theme':'redux-dark', 'themeVariables': { 'fontFamily':'Georgia, Times New Roman, serif', 'fontSize':'60px'}, 'flowchart': {'nodeSpacing':15, 'rankSpacing': 15}}}%%


flowchart LR
    %% Inputs
    subgraph INPUT[" "]
        direction TB
        SF@{ shape: tag-rect, label: "Samples File"}
        FASTQ@{ shape: processes, label: "Raw FASTQ"}
        TAIR@{ shape: div-rect, label: "Genome annotation"}
        TEs@{ shape: div-rect, label: "TEs annotation"}
        GeneDesc@{ shape: div-rect, label: "Gene descriptions"}
    end
    
    %% Pre-processing
    subgraph PREPROC[" "]
        direction LR
        Trim@{ shape: bow-rect, label: "Trim reads"}
        Bismark@{ shape: bow-rect, label: "Bismark alignment"}
        MethEx@{ shape: bow-rect, label: "Methylation extractor"}
    end
    
    %% CX Reports
    CX@{ shape: processes, label: "*CX_report* files"}
    
    %% Loading Data
    LOAD@{ shape: win-pane, label: "Load methylation data"}
    
    %% QC
    QC@{ shape: tag-rect, label: "Chloroplast C→T conversion"}
        
    %% Visualizations
    subgraph VIZ[" "]
        direction TB
        TotalMeth@{ shape: rounded, label: "Total methylation"}
        PCA@{ shape: rounded, label: "PCA of samples"}
        ChrPlots@{ shape: rounded, label: "ChrPlots genome-wide"}
        SIZE@{ shape: rounded, label: "Δ <i>vs.</i> TE size"}
        DISTANCE@{ shape: rounded, label: "Δ <i>vs.</i> TE distance"}
        %%%% TE delta analyses
        %%subgraph TEdelta[" "]
        %%    direction TB
        %%    SIZE[Δ-methylation <i>vs.</i> TE size]
        %%    DISTANCE[Δ-methylation by TE distance-from-centromere]
        %%end
    end

    %% DMR Calling
    %%DMRCALL[Call DMRs<br/><i>Single/replicate samples</i><br/>]
    DMRCALL@{ shape: win-pane, label: "DMR calling"}
    
    %% DMR Annotation
    subgraph ANNOT[" "]

        %% DMR vissualizations
        subgraph DIST[" "]
            direction TB    
            DMRdist@{ shape: rounded, label: "DMR distribution"}
            DMRpie@{ shape: rounded, label: "Gain/Loss"}
            DMRchrplot@{ shape: rounded, label: "ChrPlots DMRs"}
        end

        direction TB
        DMRannot@{ shape: win-pane, label: "DMR annotation"}
        GbDMR@{ shape: win-pane, label: "<i>Promoter/Gene</i> overlapping DMRs"}
        GbFeature@{ shape: tag-rect, label: "<i>Promoter/CDS/Intron/ 5'UTR/3'UTR/TEG/pseudogene</i>"}
        teDMR@{ shape: win-pane, label: "Transposable element DMRs"}
        GROUP@{ shape: tag-rect, label: "Functional groups"}        
        subgraph TEdown[" "]
            direction TB
            SFf@{ shape: bow-rect, label: "Super-family frequency"}
            SFd@{ shape: bow-rect, label: "Super-family distribution"}
        end

    end
    
    %% Meta Plots
    subgraph META[" "]
        direction TB
        Meta@{ shape: tag-rect, label: "MetaPlots"}    
        %% dH Analysis
        dH@{ shape: tag-rect, label: "ΔH"}
    end
    
    %% Functional Analysis
    subgraph func[ ]
        direction LR
        GO@{ shape: tag-rect, label: "GO enrichment"}
        KEGG@{ shape: tag-rect, label: "KEGG pathways"}
    end

    %% Flow
    FASTQ ------------> |"<i>optional</i>"| PREPROC
    Trim ======> Bismark ======> MethEx ==========> CX
    SF e1@==> CX
    %%LOAD ==> TEdelta
    CX ========> LOAD
    LOAD e2@==> DMRCALL
    LOAD --> |"<i>Visualizations</i>"| VIZ
    LOAD --> |"<i>optional</i>"| META
    LOAD ==> |"<i>QC</i>"| QC
    TEs -.-> teDMR
    TAIR e6@-.-> LOAD
    TEs -.-> VIZ

    %%TAIR ==> GbFeature
    DMRCALL --> |"<i>Visualizations</i>"| DIST
    DMRCALL e3@==> DMRannot
    DMRannot e4@==> GbDMR ===> GbFeature --> |"<i>optional</i>"| func
    %% DMRannot ==> wig@{ shape: tag-rect, label: "*.wig* files output"}
    GbDMR ==> GROUP
    TAIR e7@-.-> DMRannot
    GeneDesc e8@-.-> DMRannot
    DMRannot e5@==> teDMR ==> TEdown

    e1@{ animate: true }
    e2@{ animate: true }
    e3@{ animate: true }
    e4@{ animate: true }
    e5@{ animate: true }
    e6@{ animate: true }
    e7@{ animate: true }
    e8@{ animate: true }
````


---

## What the pipeline produces

For each contrast (treatment vs control), the main workflow can generate:

- **Chloroplast conversion rate**
- **PCA** (replicates only; CG/CHG/CHH + all contexts)
- **Genome-wide methylation levels** (+ eu/heterochromatin partitions)
- **Chromosome methylation profiles** (genome-wide tracks; including per-replicate “subCX” plots when available)
- **DMR calling** (CG/CHG/CHH; beta regression for replicates, Fisher’s exact test for single samples)
- **DMR chromosome maps** + **circular density (circos-like) plot**
- **DMR annotation** to:
  - Genes
  - Promoters
  - Gene features (CDS / introns / 5'UTR / 3'UTR)
  - Transposable elements (TEs)
- **TE additional summaries**:
  - TE super-family frequency
  - Coding genes vs TE overlap summaries
  - **TE Δmethylation vs TE size** (scatter; per context)
  - **TE Δmethylation vs distance from centromere** (scatter; windowed, e.g. 1 Mbp bins)
- **Transcription factors motif analysis**
- **GO enrichment**
- **KEGG enrichment**
- **Meta-plots** (Genes / TEs / Gene-features; plus delta meta-plots)
- **bigWig export** for genome browsers (DMR/Δ tracks)
- **Analyze differentially methylated vallies (1kbp)**
- **dH / surprisal module** (deltaH folder with summary plots + annotations)

---

## System requirements

### Conda (recommended)

- Linux environment (WSL works)
- [Conda / Miniconda](https://docs.conda.io/en/latest/miniconda.html) ([download](https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh)) 
- [`whiptail`](https://linux.die.net/man/1/whiptail) (for UI tutorial) (UI mode)

### Local R environment

- R ≥ 4.4.0
- System dependencies for some R packages (fonts/harfbuzz/freetype/xml, etc.)
- R packages:

```text
dplyr
tidyr
ggplot2
data.table
lattice
PeakSegDisk
geomtextpath
parallel
BiocManager
RColorBrewer
circlize
cowplot
knitr
kableExtra
DMRcaller
rtracklayer
topGO
KEGGREST
Rgraphviz
org.At.tair.db
GenomicFeatures
plyranges
```

---

## Installation

### 1) Download the source code

```bash
git clone https://github.com/Yo-yerush/Methylome.At.git
cd ./Methylome.At
```

### 2) Setup the conda environment

**Using the built-in setup script**:

```bash
chmod +x ./setup_env.sh
./setup_env.sh
```

- Use `--check` flag to check if **R** and **conda** pckages are installed

- Use `--permission` flag to ensure permission of scripts (`dos2unix`, `chmod`) **without installation**


**Or manually** (example):

```bash
packages=("r-curl" "r-rcurl" "zlib" "r-textshaping" "harfbuzz" "fribidi" "freetype" "libpng" "pkg-config" "libxml2" "r-xml" "bioconductor-rsamtools" "r-svglite") 

conda create --name Methylome.At_env
conda activate Methylome.At_env
conda install -c conda-forge -c bioconda r-base=4.4.2 ${packages[@]}

Rscript scripts/install_R_packages.R

chmod +x ./Methylome.At_UI.sh
chmod +x ./scripts/*.At.sh
```

---

## Input files

### 1) Samples table file
- `tab` or `,` delimited
- no header

Each line is:

```text
<group_name>    /path/to/methylation_calls
```

Example:

```text
wt      /data/wt_rep1.CX_report.txt
wt      /data/wt_rep2.CX_report.txt
mto1    /data/mto1_rep1.CX_report.txt
mto1    /data/mto1_rep2.CX_report.txt
mto1    /data/mto1_rep3.CX_report.txt
```

- Repeats of the same `group_name` are treated as **replicates**.
- If each group has only one file, the pipeline runs in **single-sample** mode (DMR test changes accordingly).

### 2) Supported methylation call formats

Methylome.At supports:

- **Bismark `CX_report`** (WGBS)
- **Nanopore `bedMethyl`** (recommended to generate using a plant-aware caller such as deepsignal-plant; trinucleotide column is optional)
- **CGmap (`.CGmap`)** (CGmapTools output)

All will convert to `CX_report` file which contain `tab`-delimited, no header. See [columns definition](https://support.illumina.com/help/BaseSpace_App_MethylSeq_help/Content/Vault/Informatics/Sequencing_Analysis/Apps/swSEQ_mAPP_MethylSeq_CytosineReport.htm)

```text
Chr1     3563    +       0       6       CHG     CCG
Chr1     3564    +       5       2       CG      CGA
Chr1     3565    -       2       3       CG      CGG
Chr1     3571    -       0       5       CHH     CAA
Chr1     3577    +       1       5       CHH     CTA
```

You can either:
- let the pipeline auto-detect the format, or
- set `--file_type` explicitly (`CX_report`, `bedMethyl`, `CGmap`)

#### Convert CGmap → CX_report (optional)

```bash
./scripts/cgmap_2_cx.sh /path/to/input.CGmap /path/to/output_CX_report.txt
```

#### Convert bedMethyl → CX_report

```bash
# Without trinucleotide column:
./scripts/bedmethyl_2_cx.sh -i /path/to/input.bed -o output_prefix

# With trinucleotide column (genome dir needed):
./scripts/bedmethyl_2_cx.sh -i /path/to/input.bed -t /path/to/genome_dir/ -o output_prefix
```
* *genome file as `.fasta` or `.fa`*
* *trinucleotide context are **not required** for `Methylome.At` pipeline*
  
### 3) Annotation and description files

By default, Methylome.At expects:
- a genome annotation file or table (`gtf`/`gff`/`gff3`/`csv`)
- a gene description table (adds functional descriptions to outputs)
- a TE annotation file (TAIR10 “Transposable Elements” style)

If you provide custom files, ensure they contain the columns required by the annotation scripts.
(If you are unsure, start with the default files shipped in `annotation_files/` and only then replace them.)

---

## Running Methylome.At

### UI mode

```bash
./Methylome.At_UI.sh
```

### Manual mode

#### Main pipeline

```bash
./scripts/Methylome.At.sh /path/to/samples_table.txt
```

#### Usage:

```text
$ ./scripts/Methylome.At.sh --help

Usage: ./scripts/Methylome.At.sh [samples_file] [options]

Required argument:
  --samples_file                Path to samples file [required]

Optional arguments:
  --minReadsPerCytosine         Minimum reads per cytosine [default: 6]
  --n_cores                     Number of cores [default: 8]
  --image_type                  Output images format [default: 'pdf']
  --file_type                   Post-alignment file type - 'CX_report', 'bedMethyl' and 'CGmap' [default: 'CX_report' OR determine automatically]
  --annotation_file             Genome Annotation file [default: Methylome.At annotations file (TAIR10 based)]
  --description_file            Description file [default: Methylome.At description file]
  --TEs_file                    Transposable Elements file [default: TAIR10 'Transposable Elements' annotations]
  --Methylome_At_path           Path to Methylome.At [default: /home/yoyerush]

DMRs analysis arguments:
  --minProportionDiff_CG        Minimum proportion difference for CG [default: 0.4]
  --minProportionDiff_CHG       Minimum proportion difference for CHG [default: 0.2]
  --minProportionDiff_CHH       Minimum proportion difference for CHH [default: 0.1]
  --binSize                     DMRs bin size [default: 100]
  --minCytosinesCount           Minimum cytosines count [default: 4]
  --pValueThreshold             P-value (padj) threshold [default: 0.05]

  --pca                         Perform PCA for total methylation levels
  --total_methylation           Total methylation bar-plot
  --CX_ChrPlot                  total methylation chromosome plot
  --TEs_distance_n_size         Analyze TEs total methylation by size and distance from centromere
  --total_meth_ann              Total methylation per genic annotations
  --TF_motifs                   Transcription factors motif analysis
  --func_groups                 Functional groups genes overlap DMRs
  --GO_analysis                 Perform GO analysis over DMRs
  --KEGG_pathways               Perform KEGG pathways analysis over DMRs
  --all_analyses                Enable all analyses (sets all analysis flags above [pca --> KEGG_pathways])

  --DMRs_off                    Disable the main DMRs analysis workflow
  --strand_DMRs                 Analyze strand-specific (+/-) DMRs
  --DMVs                        Analyze differentially methylated vallies (1kbp)
  --dH                          instead of DMRs worflow (calculated by ratios [p]), analyze SurpDMRs:
                                delta-H = -(p * log2(p) + (1 - p) * log2(1 - p))

MetaPlots analysis arguments:
  --MP_TEs                      Analyze of TEs metaPlots
  --MP_Genes                    Analyze of Genes-body metaPlots
  --MP_Gene_features            Analyze Gene Features metaPlots
  --all_metaplots               Enable all metaPlots analyses
  --MP_features_bin_size        Bin-size (set only for 'Gene_features' analysis!) [default: 10]
  --metaPlot_random             Number of random genes/TEs for metaPlots. 'all' for all the coding-genes and TEs [default: 10000]
```

---


## Output folders (per contrast)

A typical output tree under `results/<treatment>_vs_<control>/`:

```text
<contrast>/
  conversion_rate.csv
  <contrast>_report.html               (auto-generated HTML report)
  *.log                                 (pipeline log)
  total_methylation_analysis/
    PCA_plots/                          (replicates only)
    methylation_levels/
    ChrPlot_CX/
      subCX/                            (replicate-level tracks)
    TE_size_n_distance/
    total_methylation_annotations/      (per-annotation methylation tables)
    TF_motifs/                          (transcription factor motif plots)
  DMR_analysis/
    DMRs_CG_<contrast>.csv
    DMRs_CHG_<contrast>.csv
    DMRs_CHH_<contrast>.csv
    DMRs_Density_<contrast>.*           (circular density plot)
    gain_OR_loss/
    ChrPlot_DMRs/
    genome_annotation/
      CG/ CHG/ CHH/                    (annotation tables + summary figures)
      TEs_addiotionnal_results/
        super_family_frequency/
      functional_groups/                (functional group overlap plots)
    DMRs_bigWig/
  Strand_Asymmetry_DMRs/                (optional)
    CG/ CHG/ CHH/                       (symmetric, hemi, conflicting tables)
  DMV_analysis/                         (optional, 1 kbp valleys)
  MetaPlots/                            (optional)
    Genes/ TEs/ Gene_features/
  deltaH/                               (optional)
    genome_annotation/
  GO_analysis/                          (optional)
  KEGG_pathway/                         (optional)
```

---

## Reports

- **log**: The pipeline writes log output during the run. See [example](https://raw.githubusercontent.com/Yo-yerush/Methylome.At/refs/heads/main/output_example/Methylome.At_log_file.log) `.log` file.
  
- **HTML report**: *When the pipeline finishes, it automatically produces a `<contrast>_report.html` file with an analyses checklist, configurations, results summary (including plots and tables), colorized log, and session info.*
Render it manually:
```r
rmarkdown::render(
  "scripts/Methylome.At_report.Rmd",
  params = list(var1 = "wt", var2 = "mt1",
                Methylome.At_path = "."),
  output_file = "mt1_vs_wt_report.html"
)
```

---

## (Optional) From FASTQ to CX_report: run_bismark.sh (WGBS)

If you start from raw WGBS FASTQ files, you can generate Bismark methylation calls and CX_report outputs using:
```text
$ ./scripts/run_bismark.sh --help

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
```

#### Requirements
This script assumes you have these available in your environment:
- bismark (and its aligner dependency, typically bowtie2)
- samtools

#### Input
- Reference genome file (as `.fasta` or `.fa`)
- Samples table for the `.fastq` files (tab-delimited, 2 columns, no header)

Format:
> sample_name   /path/to/sample_R.fastq

Paired-end example (sample name appears twice: R1 + R2):
```text
mt_1    PATH/TO/FILE/mt1_R1.fastq
mt_1    PATH/TO/FILE/mt1_R2.fastq
mt_2    PATH/TO/FILE/mt2_R1.fastq
mt_2    PATH/TO/FILE/mt2_R2.fastq
wt_1    PATH/TO/FILE/wt1_R1.fastq
wt_1    PATH/TO/FILE/wt1_R2.fastq
wt_2    PATH/TO/FILE/wt2_R1.fastq
wt_2    PATH/TO/FILE/wt2_R2.fastq
```

#### Run

- *Use `-g TAIR10` for the standart reference genome of Arabidopsis (auto-download FASTA)*
- *Add `--mat` to create a ready-to-use samples table for Methylome.At*
- *Add `--cx` to produce and Keep only `*_CX_report.txt.gz` files*

```bash
./scripts/run_bismark.sh -s samples_table.txt -g TAIR10 -n 8 --cx --mat
```

---

## Troubleshooting / common issues

- **UI issues (whiptail booleans toggling incorrectly):** if multiple menu items flip together, it usually means the *tag/value* fields were reused; ensure each menu item has a unique tag (the left column), and only display the current value in the right column.
- **Single-sample runs:** PCA is skipped; DMR testing uses Fisher’s exact test instead of beta regression.
- **Nanopore bedMethyl:** if you request trinucleotide context, you must provide the genome directory (`-t`).

---

## License

This project is licensed under the [MIT License](https://github.com/Yo-yerush/Methylome.At/blob/main/LICENSE).