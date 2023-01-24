
<!-- badges: start -->

[![R](https://github.com/Sudaraka88/BacGWES/workflows/r.yml/badge.svg)](https://github.com/Sudaraka88/BacGWES/actions)
<!-- badges: end -->

# BacGWES

## Genomewide Epistasis Analysis on Bacteria

## Installation

`BacGWES` is currently available on github. It can be installed with
`devtools`

``` r
install.packages("devtools")
devtools::install_github("Sudaraka88/BacGWES")
```

## Quick Start

To run BacGWES, all you need is a fasta alignment and the corresponding
genbank annotation file. A toy dataset is provided in the package itself
to confirm that everything is setup correctly.

> **Note** For GWES analysis, it is recommended to use an alignment with
> a large number of sequences (at least \>500). Such alignments can be
> several gigabytes in size and not suitable to bundle into an R
> package.

The toy alignment comprises the first 50kb of 400 *S. pnuemoniae*
isolates (randomy selected from the 616 reported
<a href="https://www.nature.com/articles/ng.2625" target="_blank">here</a>).
The reads were aligned against the
<a href="https://www.ncbi.nlm.nih.gov/nuccore/NC_011900.1" target="_blank">ATCC
700669 reference genome</a> using
<a href="https://github.com/tseemann/snippy" target="_blank">snippy</a>.
Since this is only a part of an alignment, set \<check_gbk_fasta_lengths
= F\>, forcing BacGWES to ignore the mismatch in sequence lengths
between the genbank reference and the fasta alignment.

``` r
# devtools::install_github("Sudaraka88/BacGWES")
library(BacGWES)
dset <- "sample"
aln_path <- system.file("extdata", "sample.aln.gz", package = "BacGWES")
gbk_path <- system.file("extdata", "sample.gbk", package = "BacGWES")
snp_filt_method = "relaxed"
BacGWES(dset = dset, aln_path = aln_path, gbk_path = gbk_path, check_gbk_fasta_lengths = F,
        num_clusts_cds = 2, SnpEff_Annotate = F, snp_filt_method = snp_filt_method, 
        tanglegram_break_segments = 1)
```
