#' parse_genbank_file
#'
#' Function to parse the genbank file for the fasta alignment.
#'
#' @importFrom genbankr readGenBank getSeq
#'
#' @param gbk_path path to genbank file
#' @param g sequence length, available from the LDWeaver::parse_fasta_alignment() output (default = NULL),
#' required if <length_check = T>
#' @param length_check specify whether to check if fasta and gbk sequence lengths are equal (default = T)
#'
#' @return GenBankRecord object
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "sample.aln", package = "LDWeaver")
#' snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
#' gbk_path <- system.file("extdata", "sample.gbk", package = "LDWeaver")
#' gbk <- parse_genbank_file(gbk_path, snp.dat$g)
#' }
#' @export
parse_genbank_file = function(gbk_path, g = NULL, length_check = T){
  t0 = Sys.time()
  # Check inputs
  if(length_check){ # perform the length check
    if(is.null(g)){
      # then g cannot be null
      stop("g must be provided to perform length check!\n")
      # return(-1)
    }
  }

  if(!file.exists(gbk_path)) stop(paste("Can't locate file", gbk_path))

  # gbk = suppressWarnings(genbankr::import(gbk_path))
  gbk = suppressWarnings(genbankr::readGenBank(gbk_path))
  refseq = genbankr::getSeq(gbk)

  if(length(refseq) != 1){
    stop("The GBK file should contain the reference sequence!\n")
    # return(-1)
  }
  # the length check is good to check if the alignment matches with the gbk, setting it to F will stop it
  ref_g = length(refseq[[1]]) # length of the reference sequence

  if(length_check){ # perform the length check
    if(ref_g != g){ # perform check
      stop("Genbank reference sequence length mismatches with the fasta alignment!\n")
      # return(-1)
    }

  } else {
    if(!is.null(g)){
      if(ref_g != g){
        warning("Fasta length does not match the genbank reference sequence length!\n")
      }
    } else {
      warning("Similarity between the genbank reference and fasta sequences NOT checked, ignore if <pos> was provided...\n")
    }
  }
  cat(paste("Successfully read gbk file:", gbk_path, "in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))

  # genbankr::seqinfo(gbk)
  return(list(gbk = gbk,
              ref_g = ref_g))
}
