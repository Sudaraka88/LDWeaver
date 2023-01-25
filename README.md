
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

To run BacGWES, all you need is a fasta alignment (can be gz compressed)
and the corresponding genbank annotation file with the reference
sequence. A toy dataset is provided in the package itself to confirm
that everything is setup correctly.

> **Note** For GWES analysis, it is recommended to use an alignment with
> a large number of sequences (\>500). Such alignments can be several
> gigabytes in size and not suitable to bundle into an R package.

The toy alignment comprises the first 50kb of 400 *S. pnuemoniae*
isolates (randomy selected from the 616 reported
<a href="https://www.nature.com/articles/ng.2625" target="_blank">here</a>).
The reads were aligned against the
<a href="https://www.ncbi.nlm.nih.gov/nuccore/NC_011900.1" target="_blank">ATCC
700669 reference genome</a> using
<a href="https://github.com/tseemann/snippy" target="_blank">snippy</a>.
Since this is only a part of an alignment, set \<check_gbk_fasta_lengths
= F\>, forcing BacGWES to ignore the mismatch in sequence lengths
between the genbank reference and the fasta alignment. Several
additional options are also used, for information on all available
options run: `help(package="BacGWES")`.

``` r
# devtools::install_github("Sudaraka88/BacGWES")
library(BacGWES)
dset <- "sample"
aln_path <- system.file("extdata", "sample.aln.gz", package = "BacGWES")
gbk_path <- system.file("extdata", "sample.gbk", package = "BacGWES")
snp_filt_method = "relaxed"
BacGWES(dset = dset, aln_path = aln_path, gbk_path = gbk_path, check_gbk_fasta_lengths = F,
        num_clusts_CDS = 2, SnpEff_Annotate = F, snp_filt_method = snp_filt_method, discard_MI_threshold_lr = 0.11)
```

## Basic Outputs

If the function performs as expceted, all generated output will be saved
to a folder called “sample” in the current working directory, which can
be queried using: `getwd()`. In the example, following files will be
generated in “sample”:

- Figures

  1.  cX_fit.png - shows the distribution and modelling of the
      background linkage disequilibrium (estimated using weighted Mutual
      Information) vs. bp-separation within each cluster (two clusters
      in this example because \<num_clust_CDS=2\>)
  2.  CDS_clustering.png - shows the genome partitioning (2 clusters in
      this case), based on the CDS diversity (compared to the reference
      sequence)
  3.  sr_gwes_clust.png - short-range GWES plot for each cluster (2 in
      this case)
  4.  sr_gwes_combi.png - combined short-range GWES plot (for links with
      bp positions spanning two clusters, the max srp_value is used)
  5.  lr_gwes.png - Long range GWES plot (similar to the output from
      <a href="https://github.com/santeripuranen/SpydrPick" target="_blank">SpydrPick</a>)

- Outputs

  6.  sr_links.tsv - tab separated file containing details on
      short-range links (i.e. links \<= sr_dist bp apart)
  7.  lr_links.tsv - tab separated file containing details on long-range
      links (i.e. links \> sr_dist bp apart)

> **Note** The default sr_dist value in BacGWES is 20000bp, this can be
> modified using the \<sr_dist\> option.

- Temporary files - used to avoid costly recomputations. These files
  should be deleted before repeating an analysis using the same \<dset\>
  name.

  5.  cds_var.rds - list comprising alignment diversity information
  6.  hdw.rds - named vector comprising Hamming distance weights for
      each sequence
  7.  parsed_gbk.rds - GenBankRecord of the genbank annotation data
  8.  snp_ACGTN.rds - list comprising sparse SNP data from the alignment

## Performing Annotations

Additionally, BacGWES has an interface to perform detailed annotations
using
<a href="https://pcingola.github.io/SnpEff/" target="_blank">SnpEff</a>.
Once downloaded, set the two options: \<SnpEff_Annotate = T\> and
\<snpeff_jar_path = path_to\_\[snpEff.jar\]\_file\>.

> **Note** Since the genbank annotation file is provided, it is
> unnecessary to download and setup any snpEff databases, only the
> \<snpEff.jar\> file is required to perform annotations.

This will create additional outputs in the \<dset\> folder.

- Outputs

  9.  annotations.tsv - tab separated file containing full snpEff
      annotations on each site associated with a short-range GWES link
      with srp_max \> srp_cutoff
  10. sr_links_annotated.tsv - tab separated file similar to
      sr_links.tsv with additional snpEff annotation and allele
      distribution information
  11. tophits.tsv - tab separated file containing the top <max_tophits>
      links. Several filters are applied to extract the top links from
      sr_links_annotated.tsv
  12. Tanglegram - folder compirising html tanglegrams to easily
      visualise links and the corresponding genomic regions
  13. GWESExplorer - folder containing the outputs necessary to explore
      short-range GWES links using
      <a href="https://github.com/jurikuronen/GWES-Explorer" target="_blank">GWESExplorer</a>)

> **Note** The default srp_cutoff is 3 (i.e. p=0.001). Short-range links
> with p\>0.001 are automatically discarded, this can be modified using
> the \<srp_cutoff\> option. The default max_tophits value is 250, this
> can be modified using the \<max_tophits\> option.

- Temporary files created for snpEff annotations

  14. snpEff_data - data folder for snpEff
  15. snpEff.config - configuration file for snpEff
  16. sr_annotataed_stats.genes.txt - annotations and statistics in tab
      separated format
  17. sr_annotated_stats.html - annotations and statistics in html
      format
  18. sr_snps.vcf, sr_snps_ann.vcf - input and output from the snpEff
      annotation pipeline

## Advanced Options

While `BacGWES::BacGWES()` one-liner should work for most analyses, it
is possible to write a customised pipeline using the functions available
in the package. For a full list of available functions, run:
`help(package="BacGWES")`.

## Detailed Workthrough using Real Data
