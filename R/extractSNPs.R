#' extractSNP_from_MSA
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
#' @param mega_dset set to TRUE for mega scale datasets (default = F)
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
parse_fasta_alignment <- function(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default", mega_dset = F){

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



  if(mega_dset){ # Using SPAM
    if(!requireNamespace("spam") & !requireNamespace("spam64")){
      message("This feature requires spam and spam64 packages.")
      return(invisible())
    } else {
      # We need to make sure we are using spam64, set it quietly
      if(!getOption("spam.force64")) options(spam.force64 = T)
      snp.matrix_A <- spam::spam(list(i=snp.data$i_A, j=snp.data$j_A, values=as.logical(snp.data$x_A)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL

      snp.matrix_C <- spam::spam(list(i=snp.data$i_C, j=snp.data$j_C, values=as.logical(snp.data$x_C)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL

      snp.matrix_G <- spam::spam(list(i=snp.data$i_G, j=snp.data$j_G, values=as.logical(snp.data$x_G)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL

      snp.matrix_T <- spam::spam(list(i=snp.data$i_T, j=snp.data$j_T, values=as.logical(snp.data$x_T)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL

      snp.matrix_N <- spam::spam(list(i=snp.data$i_N, j=snp.data$j_N, values=as.logical(snp.data$x_N)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)

      # snp.matrix_A <- spam::spam(list(i=snp.data$i_A, j=snp.data$j_A, values= snp.data$x_A),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL
      #
      # snp.matrix_C <- spam::spam(list(i=snp.data$i_C, j=snp.data$j_C, values=snp.data$x_C),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL
      #
      # snp.matrix_G <- spam::spam(list(i=snp.data$i_G, j=snp.data$j_G, values=snp.data$x_G),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL
      #
      # snp.matrix_T <- spam::spam(list(i=snp.data$i_T, j=snp.data$j_T, values=snp.data$x_T),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL
      #
      # snp.matrix_N <- spam::spam(list(i=snp.data$i_N, j=snp.data$j_N, values=snp.data$x_N),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)

      snp.data = NULL
    }
  }
  else {
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

    # snp.matrix_A <- Matrix::sparseMatrix(i=snp.data$i_A,
    #                                      j=snp.data$j_A,
    #                                      x=snp.data$x_A,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL
    #
    # snp.matrix_C <- Matrix::sparseMatrix(i=snp.data$i_C,
    #                                      j=snp.data$j_C,
    #                                      x=snp.data$x_C,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL
    #
    # snp.matrix_G <- Matrix::sparseMatrix(i=snp.data$i_G,
    #                                      j=snp.data$j_G,
    #                                      x=snp.data$x_G,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL
    #
    # snp.matrix_T <- Matrix::sparseMatrix(i=snp.data$i_T,
    #                                      j=snp.data$j_T,
    #                                      x=snp.data$x_T,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL
    #
    # snp.matrix_N <- Matrix::sparseMatrix(i=snp.data$i_N,
    #                                      j=snp.data$j_N,
    #                                      x=snp.data$x_N,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    snp.data = NULL

  }

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A), snp.matrix_C=Matrix::t(snp.matrix_C),
              snp.matrix_G=Matrix::t(snp.matrix_G), snp.matrix_T=Matrix::t(snp.matrix_T),
              snp.matrix_N=Matrix::t(snp.matrix_N), g = snp.param$seq.length, nsnp = snp.param$num.snps,
              nseq = snp.param$num.seqs, seq.names = seq.names, r = rowSums(uqe), uqe = uqe, POS = snp.param$pos))
}


#' extractSNPs_from_SNP_MSA
#'
#' Function to extract SNPs from a SNP only multi-fasta alignment. This will be the output from running snp-sites.
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix sparseMatrix t
#'
#' @param aln_path path to multi fasta SNP alignment
#' @param pos numeric vector containing positions in the original alignment
#' @param gap_freq sites with a gap frequency >gap_greq will be dropped (default = 0.15)
#' @param maf_freq sites with a minor allele frequency <maf_freq will be dropped (default = 0.01)
#' @param method specify the filtering method 'relaxed' or 'default' (default = 'default')
#' @param mega_dset set to TRUE for mega scale datasets (default = F)
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
parse_fasta_SNP_alignment <- function(aln_path, pos, gap_freq = 0.15, maf_freq = 0.01, method = "default", mega_dset = F){

  ## Update to accept a SNP alignment with a positions files (output from snp-sites or https://github.com/sudaraka88/FastaR)
  # We also input pos <numeric vector> here 2023/10/12

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

  if(mega_dset){ # Using SPAM
    if(!requireNamespace("spam") & !requireNamespace("spam64")){
      message("This feature requires spam and spam64 packages.")
      return(invisible())
    } else {

      # We need to make sure we are using spam64, set it quietly
      if(!getOption("spam.force64")) options(spam.force64 = T)

      snp.matrix_A <- spam::spam(list(i=snp.data$i_A, j=snp.data$j_A, values=as.logical(snp.data$x_A)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL

      snp.matrix_C <- spam::spam(list(i=snp.data$i_C, j=snp.data$j_C, values=as.logical(snp.data$x_C)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL

      snp.matrix_G <- spam::spam(list(i=snp.data$i_G, j=snp.data$j_G, values=as.logical(snp.data$x_G)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL

      snp.matrix_T <- spam::spam(list(i=snp.data$i_T, j=snp.data$j_T, values=as.logical(snp.data$x_T)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL

      snp.matrix_N <- spam::spam(list(i=snp.data$i_N, j=snp.data$j_N, values=as.logical(snp.data$x_N)),
                                 nrow = snp.param$num.seqs, ncol = snp.param$num.snps)

      # snp.matrix_A <- spam::spam(list(i=snp.data$i_A, j=snp.data$j_A, values=snp.data$x_A),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL
      #
      # snp.matrix_C <- spam::spam(list(i=snp.data$i_C, j=snp.data$j_C, values=snp.data$x_C),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL
      #
      # snp.matrix_G <- spam::spam(list(i=snp.data$i_G, j=snp.data$j_G, values=snp.data$x_G),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL
      #
      # snp.matrix_T <- spam::spam(list(i=snp.data$i_T, j=snp.data$j_T, values=snp.data$x_T),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      # snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL
      #
      # snp.matrix_N <- spam::spam(list(i=snp.data$i_N, j=snp.data$j_N, values=snp.data$x_N),
      #                            nrow = snp.param$num.seqs, ncol = snp.param$num.snps)
      snp.data = NULL
    }
  }
  else {
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

    # snp.matrix_A <- Matrix::sparseMatrix(i=snp.data$i_A,
    #                                      j=snp.data$j_A,
    #                                      x=snp.data$x_A,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_A = snp.data$j_A = snp.data$x_A = NULL
    #
    # snp.matrix_C <- Matrix::sparseMatrix(i=snp.data$i_C,
    #                                      j=snp.data$j_C,
    #                                      x=snp.data$x_C,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_C = snp.data$j_C = snp.data$x_C = NULL
    #
    # snp.matrix_G <- Matrix::sparseMatrix(i=snp.data$i_G,
    #                                      j=snp.data$j_G,
    #                                      x=snp.data$x_G,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_G = snp.data$j_G = snp.data$x_G = NULL
    #
    # snp.matrix_T <- Matrix::sparseMatrix(i=snp.data$i_T,
    #                                      j=snp.data$j_T,
    #                                      x=snp.data$x_T,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    # snp.data$i_T = snp.data$j_T = snp.data$x_T = NULL
    #
    # snp.matrix_N <- Matrix::sparseMatrix(i=snp.data$i_N,
    #                                      j=snp.data$j_N,
    #                                      x=snp.data$x_N,
    #                                      dims = c(snp.param$num.seqs, snp.param$num.snps),
    #                                      dimnames = list(seq.names, snp.param$pos))
    snp.data = NULL
  }
  # g is wrong, change to NULL

  return(list(snp.matrix_A=Matrix::t(snp.matrix_A), snp.matrix_C=Matrix::t(snp.matrix_C),
              snp.matrix_G=Matrix::t(snp.matrix_G), snp.matrix_T=Matrix::t(snp.matrix_T),
              snp.matrix_N=Matrix::t(snp.matrix_N), g = NULL, nsnp = snp.param$num.snps,
              nseq = snp.param$num.seqs, seq.names = seq.names, r = rowSums(uqe), uqe = uqe, POS = snp.param$pos))
}
