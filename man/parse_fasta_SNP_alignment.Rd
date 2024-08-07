% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractSNPs.R
\name{parse_fasta_SNP_alignment}
\alias{parse_fasta_SNP_alignment}
\title{extractSNPs_from_SNP_MSA}
\usage{
parse_fasta_SNP_alignment(
  aln_path,
  pos,
  gap_freq = 0.15,
  maf_freq = 0.01,
  method = "default",
  mega_dset = F
)
}
\arguments{
\item{aln_path}{path to multi fasta SNP alignment}

\item{pos}{numeric vector containing positions in the original alignment}

\item{gap_freq}{sites with a gap frequency >gap_greq will be dropped (default = 0.15)}

\item{maf_freq}{sites with a minor allele frequency <maf_freq will be dropped (default = 0.01)}

\item{method}{specify the filtering method 'relaxed' or 'default' (default = 'default')}

\item{mega_dset}{set to TRUE for mega scale datasets (default = F)}
}
\value{
R list with SNPs in sparse format and additional parameters
}
\description{
Function to extract SNPs from a SNP only multi-fasta alignment. This will be the output from running snp-sites.
}
\examples{
\dontrun{
aln_path <- system.file("extdata", "sample.aln", package = "LDWeaver")
snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
}
}
