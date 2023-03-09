#' extractSNPs
#'
#' Function to extract SNPs from a multi-fasta alignment.
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix sparseMatrix t
#'
#' @param aln_path path to multi fasta alignment
#' @param gap_freq sites with a gap frequency >gap_greq will be dropped (default = 0.15)
#' @param maf_freq sites with a minor allele frequency <maf_freq will be dropped (default = 0.01)
#' @param method specify the filtering method 'relaxed' or 'default' (default = 'default')
#'
#' @return R list with SNPs in sparse format and additional parameters
#' @useDynLib LDWeaver
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "sample.aln", package = "LDWeaver")
#' snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
#' }
#' @export
parse_fasta_alignment <- function(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default"){
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

  # t0 = Sys.time() # This timer is not accurate, return function has computations to save memory, time externally instead

  # snp.data <- getSNP_Sites(aln_path, filter, gap_freq, maf_freq)
  # snp.data <- .getACGTN_Sites(aln_path, filter, gap_freq, maf_freq)
  snp.param <- .extractAlnParam(aln_path, filter, gap_freq, maf_freq)

  # g = snp.data$params$seq.length
  # nseq = snp.data$params$num.seqs
  # g = snp.param$seq.length; snp.param$seq.length = NULL
  # nseq = snp.param$num.seqs; snp.param$num.seqs = NULL
  # nsnp = snp.param$num.snps; snp.param$num.snps = NULL
  # print(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))

  if(snp.param$seq.length==-1) stop("Error! sequences are of different lengths!")
  if(snp.param$num.seqs==0) stop("File does not contain any sequences!")
  if(snp.param$num.snps==0) stop("File does not contain any SNPs")

  snp.data <- .extractSNPs(aln_path, snp.param$num.seqs, snp.param$num.snps, snp.param$pos)
  # seq.names <-  gsub("^>","",snp.data$params$seq.names)
  seq.names <- gsub("^>","",snp.data$seq.names); snp.data$seq.names = NULL
  uqe = apply(snp.data$ACGTN_table>0, 1, function(x) as.numeric(x>0)); snp.data$ACGTN_table = NULL

  # nsnp = snp.data$params$num.snps


  snp.matrix_A <- Matrix::sparseMatrix(i=snp.data$i_A,
                                       j=snp.data$j_A,
                                       x=as.logical(snp.data$x_A),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL

  snp.matrix_C <- Matrix::sparseMatrix(i=snp.data$i_C,
                                       j=snp.data$j_C,
                                       x=as.logical(snp.data$x_C),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL

  snp.matrix_G <- Matrix::sparseMatrix(i=snp.data$i_G,
                                       j=snp.data$j_G,
                                       x=as.logical(snp.data$x_G),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL

  snp.matrix_T <- Matrix::sparseMatrix(i=snp.data$i_T,
                                       j=snp.data$j_T,
                                       x=as.logical(snp.data$x_T),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL

  snp.matrix_N <- Matrix::sparseMatrix(i=snp.data$i_N,
                                       j=snp.data$j_N,
                                       x=as.logical(snp.data$x_N),
                                       dims = c(snp.param$num.seqs, snp.param$num.snps),
                                       dimnames = list(seq.names, snp.param$pos))
  snp.data = NULL


  # cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A), snp.matrix_C=Matrix::t(snp.matrix_C),
              snp.matrix_G=Matrix::t(snp.matrix_G), snp.matrix_T=Matrix::t(snp.matrix_T),
              snp.matrix_N=Matrix::t(snp.matrix_N), g = snp.param$seq.length, nsnp = snp.param$num.snps,
              nseq = snp.param$num.seqs, seq.names = seq.names, r = rowSums(uqe), uqe = uqe, POS = snp.param$pos))
}
