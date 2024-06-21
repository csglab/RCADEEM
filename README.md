# RCADEEM: Recognition Code-Assisted Discovery of regulatory Elements by Expectation-Maximization

RCADEEM uses a hidden Markov model (HMM) to represent multiple, alternative DNA motifs, each corresponding to the binding preference of a zinc finger array.

## Requirements
- Unix-compatible OS.
- [Python 3.6 or later](https://www.python.org/).
- [Bedtools](https://bedtools.readthedocs.io/en/latest/index.html).
- [R version 3.0.1 or later](http://www.r-project.org/).
 
- R libraries:
[assertr](https://cran.r-project.org/web/packages/assertr/vignettes/assertr.html), [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [ggh4x](https://github.com/teunbrand/ggh4x), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggseqlogo](https://github.com/omarwagih/ggseqlogo), [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html), [gplots](https://cran.r-project.org/web/packages/gplots/index.html), [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html), [matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html), [pROC](https://cran.r-project.org/web/packages/pROC/index.html), [randomForest](http://cran.r-project.org/web/packages/randomForest/index.html), [readr](https://cran.r-project.org/web/packages/readr/index.html), [reshape](https://cran.r-project.org/web/packages/reshape/index.html), [stringdist](https://cran.r-project.org/web/packages/stringdist/index.html), [stringi](https://cran.r-project.org/web/packages/stringi/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), [circlize](https://cran.r-project.org/web/packages/circlize/index.html), and [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html).

- [GNU-compatible MAKE software](https://gcc.gnu.org/).
- [MEME Suite](http://meme.nbcr.net/meme/downloads.html).
 
## Installation
- Step 1. Run the `make` command.
- Step 2. Change the value of line # of the “RCADEEM” script to the path to the executable MEME files (`fasta-center` and `fasta-dinucleotide-shuffle`). Alternatively, you can provide the path via the `--meme_lib_exec_meme_dir` argument.

## Usage
- `-h, --help` Show this help message and exit.

- `--meme_lib_exec_meme_dir MEME_DIR` 

- `--task [TASK ...]` Task: RCADEEM and/or HEATMAP (separated by a comma) the former runs the main algoritm and the latter visualizes the binding modes as a heatmap of the original sequences.

- `--job_ID JOB_ID` Experiment ID. ex: CTCF_rep2_HEK293

- `--out_dir OUT_DIR` Output directory.

- `--C2H2_ZFP_fasta ZFP_FA` Fasta file of the C2H2-ZF protein sequence.

- `--target_bed TARGET_BED` BED file specifying the sequences of interest, requires --genome_fasta argument.

- `--target_fasta TARGET_FASTA` FASTA file with the sequences of interest, mutually exclusive with --target_bed and --genome_fasta.

- `--genome_fasta GENOME_FA` See --target_bed.

- `--chr_sizes CHR_SIZES` Chromosome sizes for --genome_fasta.

- `--cutoff CUTOFF` In the heatmap, exclude binding modes with HMM scores with fewer than --minzise entries that pass this cutoff.

- `--minsize MINSIZE` See cuttoff.

- `--bw BW_FILE` BIGWIG file, adds a panel to the main heatmap (e.g. ChIP-seq counts). For multiple files, separate them by a comma (also applies to --bw_labels and --bw_units).

- `--bw_labels BW_LABELS` Labels for the --bw heatmap panel.

- `--bw_units BW_UNITS` Units for the --bw heatmap panel.

- `--ref_repeats REF_REPEATS` File with repeat element coordinates, adds a panel to the main heatmap.

- `--rcadeem_range RC_RANGE` Range for the sequences to be included for RCADEEM, measure from the center of the fasta or bed files (--target_fasta or --target_bed).

- `--bed_has_score` The fourth column in the input bed files is a numeric score that will be shown in the aligned heatmap.

- `--footprint_file FOOTPRINT` BIGWIG with a Footprint signal. Adds a panel to the main heatmap.














