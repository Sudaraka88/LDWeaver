#' estimate_Hamming_distance_weights
#'
#' Function to estimate the similarity between sequences based on the Hamming distance
#'
#' @importFrom MatrixExtra crossprod
#' @importFrom Matrix colSums
#' @importFrom methods as
#'
#' @param snp.dat output from parsing the multi fasta alignment using LDWeaver::parse_fasta_alignment()
#' @param threshold Hamming distance similarity threshold (default = 0.1, i.e. 10\%)
#'
#' @return Vector with Hamming distance weights for each sequence that can be used to reweight Mutual Information estimates
#'
#' @examples
#' \dontrun{
#' hdw <- estimate_Hamming_distance_weights(snp.dat, threshold=0.1)
#' }
#' @export
estimate_Hamming_distance_weights = function(snp.dat, threshold = 0.1){

  t0 = Sys.time()
  thresh =  as.integer(snp.dat$nsnp*threshold)

  snpmat_t = as(snp.dat$snp.matrix_A, 'lgeMatrix') # This crashes in linuxMint for some reason!
  shared.snps = Matrix::crossprod(snpmat_t)

  snpmat_t = as(snp.dat$snp.matrix_A, 'lgeMatrix')
  shared.snps = MatrixExtra::crossprod((snpmat_t))

  snpmat_t = as(snp.dat$snp.matrix_C, 'lgeMatrix')
  shared.snps = shared.snps + MatrixExtra::crossprod((snpmat_t))

  snpmat_t = as(snp.dat$snp.matrix_G, 'lgeMatrix')
  shared.snps = shared.snps + MatrixExtra::crossprod((snpmat_t))

  snpmat_t = as(snp.dat$snp.matrix_T, 'lgeMatrix')
  shared.snps = shared.snps + MatrixExtra::crossprod((snpmat_t))

  snpmat_t = as(snp.dat$snp.matrix_N, 'lgeMatrix')
  shared.snps = shared.snps + MatrixExtra::crossprod((snpmat_t))

  hdw = 1/Matrix::colSums((snp.dat$nsnp - shared.snps) < thresh)
  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  return(hdw)
}

# # Write a slower method to handle larger dsets - avoid for now
# perform_pop_struct_correction = function(snp.dat, threshold = 0.1){
#   thresh =  as.integer(snp.dat$nsnp*threshold)
#
#   t0 =Sys.time()  # will take a minute+ (depending on the dataset)
#   hdw = foreach(i = 1:nseq, .combine = "c") %dopar% hdwC_M(seqs,i-1, thresh, snp.dat$nseq)
#   hdw = 1/hdw
#   print(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
#   return(hdw)
# }
