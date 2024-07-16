#' estimate_Hamming_distance_weights
#'
#' Function to estimate the similarity between sequences based on the Hamming distance
#'
#' @importMethodsFrom MatrixExtra crossprod
#' @importFrom Matrix colSums
#' @importFrom methods as
#'
#' @param snp.dat output from parsing the multi fasta alignment using LDWeaver::parse_fasta_alignment()
#' @param threshold Hamming distance similarity threshold (default = 0.1, i.e. 10\%)
#' @param mega_dset set to TRUE for mega scale datasets (default = F)
#'
#' @return Vector with Hamming distance weights for each sequence that can be used to reweight Mutual Information estimates
#'
#' @examples
#' \dontrun{
#' hdw <- estimate_Hamming_distance_weights(snp.dat, threshold=0.1)
#' }
#' @export
estimate_Hamming_distance_weights = function(snp.dat, threshold = 0.1, mega_dset = F){

  t0 = Sys.time()
  thresh =  as.integer(snp.dat$nsnp*threshold)

  if(mega_dset){ # Using SPAM
    if(!requireNamespace("spam") & !requireNamespace("spam64")){
      message("This feature requires spam and spam64 packages.")
      return(invisible())
    } else {

      # shared.snps.spam = spam::crossprod(snp.dat$snp.matrix_A)
      # shared.snps.spam = shared.snps.spam + spam::crossprod(snp.dat$snp.matrix_C)
      # shared.snps.spam = shared.snps.spam + spam::crossprod(snp.dat$snp.matrix_G)
      # shared.snps.spam = shared.snps.spam + spam::crossprod(snp.dat$snp.matrix_T)
      # shared.snps.spam = shared.snps.spam + spam::crossprod(snp.dat$snp.matrix_N)

      shared.snps.spam = crossprod(snp.dat$snp.matrix_A)
      shared.snps.spam = shared.snps.spam + crossprod(snp.dat$snp.matrix_C)
      shared.snps.spam = shared.snps.spam + crossprod(snp.dat$snp.matrix_G)
      shared.snps.spam = shared.snps.spam + crossprod(snp.dat$snp.matrix_T)
      shared.snps.spam = shared.snps.spam + crossprod(snp.dat$snp.matrix_N)

      ## Is this a bug in the spam package? This should work:
      # hdw.spam = 1/( spam::colSums( (snp.dat$nsnp - shared.snps.spam) < thresh) + 1)
      # After last run line above, both "shared.snps.spam" and "shared.snps" have the same values

      # Have to manually adjust the values in shared.snps.spam@entries as below to get spam::colSums to get the same output
      shared.snps.spam@entries = snp.dat$nsnp - shared.snps.spam@entries # Get the residual
      shared.snps.spam@entries = sapply(shared.snps.spam@entries, function(x) ifelse( (x < thresh), 1, 0)) # Check if it's < thresh, allocate a numeric value
      # It is possible spam requires the values to be numeric and not logical for spam::colSums to work - check later (2024/04/22)
      hdw = 1/( spam::colSums(shared.snps.spam) + 1)

    }
  } else {
    snpmat_t = as(snp.dat$snp.matrix_A, 'dgCMatrix')
    # shared.snps = MatrixExtra::crossprod((snpmat_t))
    shared.snps = MatrixExtra::crossprod((snpmat_t))

    snpmat_t = as(snp.dat$snp.matrix_C, 'dgCMatrix')
    shared.snps = shared.snps + crossprod((snpmat_t))

    snpmat_t = as(snp.dat$snp.matrix_G, 'dgCMatrix')
    shared.snps = shared.snps + crossprod((snpmat_t))

    snpmat_t = as(snp.dat$snp.matrix_T, 'dgCMatrix')
    shared.snps = shared.snps + crossprod((snpmat_t))

    snpmat_t = as(snp.dat$snp.matrix_N, 'dgCMatrix')
    shared.snps = shared.snps + crossprod((snpmat_t))

    hdw = 1/( Matrix::colSums( (snp.dat$nsnp - shared.snps) < thresh) + 1)
  }

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  return(hdw)
}
