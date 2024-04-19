#' extractSNPs_from_megaMSA
#'
#' Function to extract SNPs from mega scale multi fasta-alignments.
#' Depending on alignment size, this function can take a very long time. Best to run snp-sites externally and use the snp-only alignment with parse_fasta_SNP_alignment_mega)
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
parse_fasta_alignment_mega <- function(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default"){

  if(!requireNamespace("spam") & !requireNamespace("spam64")){
    message("This feature requires spam and spam64 packages.")
    return(invisible())
  } else {
    message("Using spam64 with spam.force64=TRUE")
    options(spam.force64 = TRUE)
  }

  # Check inputs
  aln_path = normalizePath(aln_path) # C safety (~ character causes crash)
  if(!file.exists(aln_path)) stop(paste("Can't locate file", aln_path))

  if(method == "default") {
    filter = 0
  } else if(method == "relaxed") {
    filter = 1
  } else {
    warning("Unkown filtering method, using default...")
    filter = 0
  }

  # extract align parameters and filtering...
  snp.param <- .extractAlnParam(aln_path, filter, gap_freq, maf_freq)

  if(snp.param$seq.length==-1) stop("Error! sequences are of different lengths!")
  if(snp.param$num.seqs==0) stop("File does not contain any sequences!")
  if(snp.param$num.snps==0) stop("File does not contain any SNPs")

  snp.data <- .extractSNPs(aln_path, snp.param$num.seqs, snp.param$num.snps, snp.param$pos)
  seq.names <- gsub("^>","",snp.data$seq.names); snp.data$seq.names = NULL
  uqe = apply(snp.data$ACGTN_table>0, 1, function(x) as.numeric(x>0)); snp.data$ACGTN_table = NULL

  snp.matrix_A_spam <- spam::spam(list(i=snp.data$i_A, j=snp.data$j_A, values=as.logical(snp.data$x_A)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL

  snp.matrix_C_spam <- spam::spam(list(i=snp.data$i_C, j=snp.data$j_C, values=as.logical(snp.data$x_C)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL

  snp.matrix_G_spam <- spam::spam(list(i=snp.data$i_G, j=snp.data$j_G, values=as.logical(snp.data$x_G)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL

  snp.matrix_T_spam <- spam::spam(list(i=snp.data$i_T, j=snp.data$j_T, values=as.logical(snp.data$x_T)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL

  snp.matrix_N_spam <- spam::spam(list(i=snp.data$i_N, j=snp.data$j_N, values=as.logical(snp.data$x_N)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data = NULL

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A_spam), snp.matrix_C=Matrix::t(snp.matrix_C_spam),
              snp.matrix_G=Matrix::t(snp.matrix_G_spam), snp.matrix_T=Matrix::t(snp.matrix_T_spam),
              snp.matrix_N=Matrix::t(snp.matrix_N_spam), g = snp.param$seq.length, nsnp = snp.param$num.snps,
              nseq = snp.param$num.seqs, seq.names = seq.names, r = rowSums(uqe), uqe = uqe, POS = snp.param$pos))

}

#' extractSNPs_from_SNP_megaMSA
#'
#' Function to extract SNPs from a mega scale SNP only multi-fasta alignment. This will be the output from running snp-sites.
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix sparseMatrix t
#'
#' @param aln_path path to multi fasta SNP alignment
#' @param pos numeric vector containing positions in the original alignment
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
parse_fasta_SNP_alignment_mega <- function(aln_path, pos, gap_freq = 0.15, maf_freq = 0.01, method = "default"){

  if(!requireNamespace("spam") & !requireNamespace("spam64")){
    message("This feature requires spam and spam64 packages.")
    return(invisible())
  } else {
    message("Using spam64 with spam.force64=TRUE")
    options(spam.force64 = TRUE)
  }

  # Check inputs
  aln_path = normalizePath(aln_path) # C safety (~ character causes crash)
  if(!file.exists(aln_path)) stop(paste("Can't locate file", aln_path))

  if(method == "default") {
    filter = 0
  } else if(method == "relaxed") {
    filter = 1
  } else {
    warning("Unkown filtering method, using default...")
    filter = 0
  }

  # update .extractAlnParam to say we are working with a SNP alignment
  snp.param <- .extractAlnParam(aln_path, filter, gap_freq, maf_freq)

  if(snp.param$seq.length==-1) stop("Error! sequences are of different lengths!")
  if(snp.param$num.seqs==0) stop("File does not contain any sequences!")
  if(snp.param$num.snps==0) stop("File does not contain any SNPs")

  # we should do a sanity check here to ensure POS can be extracted
  if(length(pos) != snp.param$seq.length) stop("Error! Number of positions do not match the fasta sequence length")

  # need to extract SNPs using the snp.param positions
  snp.data <- .extractSNPs(aln_path, snp.param$num.seqs, snp.param$num.snps, snp.param$pos)

  # now we can update snp.param$pos to real ones
  snp.param$pos = as.integer(pos[snp.param$pos]) # to account for possible dropped snps

  # seq.names <-  gsub("^>","",snp.data$params$seq.names)
  seq.names <- gsub("^>","",snp.data$seq.names); snp.data$seq.names = NULL
  uqe = apply(snp.data$ACGTN_table>0, 1, function(x) as.numeric(x>0)); snp.data$ACGTN_table = NULL # uqe is not sensitive to SNP positions


  snp.matrix_A_spam <- spam::spam(list(i=snp.data$i_A, j=snp.data$j_A, values=as.logical(snp.data$x_A)),
                                  nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL

  snp.matrix_C_spam <- spam::spam(list(i=snp.data$i_C, j=snp.data$j_C, values=as.logical(snp.data$x_C)),
                                  nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL

  snp.matrix_G_spam <- spam::spam(list(i=snp.data$i_G, j=snp.data$j_G, values=as.logical(snp.data$x_G)),
                                  nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL

  snp.matrix_T_spam <- spam::spam(list(i=snp.data$i_T, j=snp.data$j_T, values=as.logical(snp.data$x_T)),
                                  nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL

  snp.matrix_N_spam <- spam::spam(list(i=snp.data$i_N, j=snp.data$j_N, values=as.logical(snp.data$x_N)),
                                  nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
  snp.data = NULL

  # g is wrong, change to NULL

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A_spam), snp.matrix_C=Matrix::t(snp.matrix_C_spam),
              snp.matrix_G=Matrix::t(snp.matrix_G_spam), snp.matrix_T=Matrix::t(snp.matrix_T_spam),
              snp.matrix_N=Matrix::t(snp.matrix_N_spam), g = NULL, nsnp = snp.param$num.snps,
              nseq = snp.param$num.seqs, seq.names = seq.names, r = rowSums(uqe), uqe = uqe, POS = snp.param$pos))
}
