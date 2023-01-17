#' estimate_variation_in_CDS
#'
#' Function to estimate the variation within each coding region, the output from this function
#' can be used to segment the genome into diversity-based clusters.
#'
#' @import genbankr
#' @importFrom GenomicRanges start width end
#' @import Matrix
#' @importFrom data.table data.table setattr %between% .I
#'
#' @param gbk output from parsing the genbank file using BacGWES::parseGenBankFile()
#' @param snp.dat output from parsing the multi fasta alignment using BacGWES::parseFastaAlignment()
#' @param cores specify the number of cores to use
#'
#' @return R list with CDS variation and allele distribution details
#'
#' @examples
#' \dontrun{
#' snp.dat <- estimate_variation_in_CDS(gbk, snp.dat, ncores=10)
#' }
#' @export
estimate_variation_in_CDS = function(gbk, snp.dat, ncores){
  # This method is only approximate, but much MUCH faster and easier on resources
  # TODO: Include the higher accuracy function
  t0 = Sys.time()
  cds_reg = genbankr::cds(gbk)
  starts = GenomicRanges::start(cds_reg)
  widths = GenomicRanges::width(cds_reg)
  ends = GenomicRanges::end(cds_reg)
  ncds = length(starts)
  var_estimate = rep(NA, ncds)
  # convert ref to a CharacterVector
  ref = unlist(unname(strsplit(as.character(genbankr::getSeq(gbk)), '')))[snp.dat$POS]

  variation = matrix(c(Matrix::rowSums(snp.dat$snp.matrix_A),
                       Matrix::rowSums(snp.dat$snp.matrix_C),
                       Matrix::rowSums(snp.dat$snp.matrix_G),
                       Matrix::rowSums(snp.dat$snp.matrix_T),
                       Matrix::rowSums(snp.dat$snp.matrix_N)),
                     ncol = snp.dat$nsnp, byrow = T)

  # Generate a reference masking mx with 0 at reference allele
  reference = matrix(rep(1, 5*snp.dat$nsnp), nrow = 5); ACGTN2num(reference, ref, ncores)

  variation_wo_ref = variation*reference # Remove ref. allele sums

  # side-job - prepare ALT allele for VCF format (required for snpEff)
  alpha = c("A", "C", "G", "T", "*")
  alt = apply(variation_wo_ref > 0, 2, function(x) paste(alpha[which(x)], sep = ',', collapse = ','))

  snp_var = Matrix::colSums(variation_wo_ref) # Variation at each SNP

  POS_dt = data.table::data.table(w = snp.dat$POS)
  data.table::setattr(POS_dt, "sorted", "w")
  for(cds in 1:ncds){
    # pos_idx = POS_dt[, .I[POS_dt$w %between% c(starts[cds], ends[cds])]]
    pos_idx = POS_dt[, .I[POS_dt$w %between% c(starts[cds], ends[cds])]]
    if(length(pos_idx) > 0) var_estimate[cds] = sum(snp_var[pos_idx])/widths[cds]
    # Divide by width to normalise between CDS (i.e. long CDS with fewer SNPs have smaller variability)
  }
  cds_idx = !is.na(var_estimate)

  rownames(variation) = c("A", "C", "G", "T", "N")

  print(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
  return(list(var_estimate = var_estimate[cds_idx], cds_start = starts[cds_idx], cds_end = ends[cds_idx],
              ref = ref, alt = alt, allele_table = variation))
}
