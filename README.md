# RCADEEM: Recognition Code-Assisted Discovery of regulatory Elements by Expectation-Maximization

RCADEEM uses a hidden Markov model (HMM) to represent multiple, alternative DNA motifs, each corresponding to the binding preference of a zinc finger array.

## Requirements
- Unix-compatible OS.
- [Python 3.6 or later.](https://www.python.org/).
- [Bedtools 2.30 or later.](https://bedtools.readthedocs.io/en/latest/index.html).
- [R version 3.0.1 or later.](http://www.r-project.org/).
 
- R libraries:
[assertr](https://cran.r-project.org/web/packages/assertr/vignettes/assertr.html), [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html), [data.table](https://cran.r-project.org/web/packages/data.table/index.html), [ggh4x](https://github.com/teunbrand/ggh4x), [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [ggseqlogo](https://github.com/omarwagih/ggseqlogo), [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html), [gplots](https://cran.r-project.org/web/packages/gplots/index.html), [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html), [matrixStats](https://cran.r-project.org/web/packages/matrixStats/index.html), [optparse](https://cran.r-project.org/web/packages/optparse/index.html), [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html), [pROC](https://cran.r-project.org/web/packages/pROC/index.html), [randomForest](http://cran.r-project.org/web/packages/randomForest/index.html), [readr](https://cran.r-project.org/web/packages/readr/index.html), [reshape](https://cran.r-project.org/web/packages/reshape/index.html), [stringdist](https://cran.r-project.org/web/packages/stringdist/index.html), [stringi](https://cran.r-project.org/web/packages/stringi/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), [circlize](https://cran.r-project.org/web/packages/circlize/index.html), and [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html).

- [GNU-compatible MAKE software](https://gcc.gnu.org/).
- [MEME Suite](http://meme.nbcr.net/meme/downloads.html).
 
## Installation

```bash
git clone https://github.com/csglab/RCADEEM.git
```
After cloning:
1. You can add the line export PATH=${cloning_directory}/RCADEEM:$PATH to your .bashrc file.
2. `cd` to the repository directory and run the `make` command.
3. Change the value of line `91` of the `RCADEEM` script to the path to the executable MEME files (`fasta-center` and `fasta-dinucleotide-shuffle`). Alternatively, you can provide the path via the `--meme_lib_exec_meme_dir` argument.

## Demo
The demo script is `RCADEEM_demo.sh` and the input files are in `./data/demo_CTCF/IN/`. The demo run time for task `RCADEEM` is ~2 hours and ~10 minutes for task `HEATMAP`. 

```bash
## Using fasta sequences directly.
RCADEEM --task RCADEEM,HEATMAP \
          --job_ID CTCF_demo_from_fasta \
          --out_dir ./data/demo_CTCF/OUT \
          --C2H2_ZFP_fasta ./data/demo_CTCF/IN/CTCF_protein_sequence.fa \
          --target_fasta ./data/demo_CTCF/IN/GSM1407629.top500summits.500bp.fasta
```

```bash
## Using the bed file and the genome fasta (requires --chr_sizes). 
RCADEEM --task RCADEEM,HEATMAP \
          --job_ID CTCF_demo_from_bed \
          --out_dir ./data/demo_CTCF/OUT \
          --C2H2_ZFP_fasta ./data/demo_CTCF/IN/CTCF_protein_sequence.fa \
          --target_bed ./data/demo_CTCF/IN/target_CTCF_coef_br_top_100.bed \
          --genome_fasta ${GENOME} \
          --chr_sizes ${CHR_SIZES}
```


## Input files

- `--C2H2_ZFP_fasta` Fasta file of the C2H2-ZF protein sequence.
- `--target_fasta` FASTA file with the sequences of interest, mutually exclusive with `--target_bed` and `--genome_fasta`.
- `--target_bed` BED file specifying the sequences of interest, requires `--genome_fasta` argument.
- `--genome_fasta`  Fasta file of the reference genome, requires `--chr_sizes`.


## **Citation**

Jolma, A., Hernandez-Corchado, A., Yang, A. W. H., Fathi, A., Laverty, K. U., Brechalov, A., Razavi, R., Albu, M., Zheng, H., Kulakovskiy, I. V., Najafabadi, H. S., & Hughes, T. R. (2024). GHT-SELEX demonstrates unexpectedly high intrinsic sequence specificity and complex DNA binding of many human transcription factors. BioRxiv, 2024.11.11.618478. https://doi.org/10.1101/2024.11.11.618478





Heatmaps are generated with the ComplexHeatmap and circlize packages. If you use them in published research, please cite:

- Gu, Z. Complex Heatmap Visualization. iMeta 2022.
or
- Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
    genomic data. Bioinformatics 2016.
and
- Gu, Z. circlize implements and enhances circular visualization
  in R. Bioinformatics 2014.
