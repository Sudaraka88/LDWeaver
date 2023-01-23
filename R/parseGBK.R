#' parseGenBankFile
#'
#' Function to parse the genbank file for the fasta alignment.
#'
#' @importFrom genbankr import getSeq
#'
#' @param gbk_path path to genbank file
#' @param g sequence length, available from the BacGWES::parse_fasta_alignment() output
#'
#' @return GenBankRecord object
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "sample.aln", package = "BacGWES")
#' snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
#' gbk_path <- system.file("extdata", "sample.gbk", package = "BacGWES")
#' gbk <- parseGenBankFile(gbk_path, snp.dat$g)
#' }
#' @export
parse_genbank_file = function(gbk_path, g){
  t0 = Sys.time()
  # Check inputs
  if(!file.exists(gbk_path)) stop(paste("Can't locate file", gbk_path))

  gbk = suppressWarnings(genbankr::import(gbk_path))
  cat(paste("Successfully read gbk file:", gbk_path, "in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
  refseq = genbankr::getSeq(gbk)
  if(length(refseq) != 1){
    cat("The GBK file should contain the reference sequence!\n")
    return(-1)
  }
  if(length(refseq[[1]]) != g){
    cat("Reference sequence length does not match with fasta file length\n")
    return(-1)
  }
  # genbankr::seqinfo(gbk)
  return(gbk)
}
