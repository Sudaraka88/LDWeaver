#' extractSNPs
#'
#' Function to extract SNPs from a multi-fasta alignment.
#'
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#'
#' @param aln_path path to multi fasta alignment
#' @param gap_freq sites with a gap frequency >gap_greq will be dropped (default = 0.15)
#' @param maf_freq sites with a minor allele frequency <maf_freq will be dropped (default = 0.01)
#' @param method specify the filtering method 'relaxed' or 'default' (default = 'default')
#'
#' @return R list with SNPs in sparse format and additional parameters
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "sample.aln", package = "BacGWES")
#' snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
#' }
#' @export
parseFastaAlignment <- function(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default"){
  # Check inputs
  if(!file.exists(aln_path)) stop(paste("Can't locate file", aln_path))

  if(method == "default") {
    filter = 0
  } else if(method == "relaxed") {
    filter = 1
  } else {
    warning("Unkown filtering method, using default...")
    filter = 0
  }

  t0 = Sys.time()
  # snp.data <- getSNP_Sites(aln_path, filter, gap_freq, maf_freq)
  snp.data <- .getACGTN_Sites(aln_path, filter, gap_freq, maf_freq)

  g = snp.data$params$seq.length
  nseq = snp.data$params$num.seqs
  # print(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))

  if(g==-1) stop("Error! sequences are of different lengths!")
  if(nseq==0) stop("File does not contain any sequences!")


  seq.names <-  gsub("^>","",snp.data$params$seq.names)
  nsnp = snp.data$params$num.snps

  snp.matrix_A <- Matrix::sparseMatrix(i=snp.data$i_A,
                                       j=snp.data$j_A,
                                       x=as.logical(snp.data$x_A),
                                       dims = c(nseq, nsnp),
                                       dimnames = list(seq.names, snp.data$pos))

  snp.matrix_C <- Matrix::sparseMatrix(i=snp.data$i_C,
                                       j=snp.data$j_C,
                                       x=as.logical(snp.data$x_C),
                                       dims = c(nseq, nsnp),
                                       dimnames = list(seq.names, snp.data$pos))

  snp.matrix_G <- Matrix::sparseMatrix(i=snp.data$i_G,
                                       j=snp.data$j_G,
                                       x=as.logical(snp.data$x_G),
                                       dims = c(nseq, nsnp),
                                       dimnames = list(seq.names, snp.data$pos))

  snp.matrix_T <- Matrix::sparseMatrix(i=snp.data$i_T,
                                       j=snp.data$j_T,
                                       x=as.logical(snp.data$x_T),
                                       dims = c(nseq, nsnp),
                                       dimnames = list(seq.names, snp.data$pos))

  snp.matrix_N <- Matrix::sparseMatrix(i=snp.data$i_N,
                                       j=snp.data$j_N,
                                       x=as.logical(snp.data$x_N),
                                       dims = c(nseq, nsnp),
                                       dimnames = list(seq.names, snp.data$pos))

  uqe = apply(snp.data$ACGTN_table>0, 1, function(x) as.numeric(x>0))
  r = rowSums(uqe)

  print(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A), snp.matrix_C=Matrix::t(snp.matrix_C),
              snp.matrix_G=Matrix::t(snp.matrix_G), snp.matrix_T=Matrix::t(snp.matrix_T),
              snp.matrix_N=Matrix::t(snp.matrix_N), g = g, nsnp = nsnp, nseq = nseq,
              seq.names = seq.names, r = r, uqe = uqe, POS = as.numeric(snp.data$pos)))
}
