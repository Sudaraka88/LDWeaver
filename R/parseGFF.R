#' parse_gff_file
#'
#' Function to parse the gff3 file for the fasta alignment.
#'
#'
#' @param gff3_path path to gff3_annotation file. This function relies on ape::read.gff() to read in the file
#' @param ref_fasta_path path to Reference fasta file. The file MUST be in fasta format and contain exactly one sequence!
#' @param perform_length_check specify whether to check if the gff3 annotations agree with the reference length (default = T)
#'
#' @return gff annotation and reference sequence as a list object
#'
#' @examples
#' \dontrun{
#' gff3_path <- "<path_to_gff3_file>"
#' ref_fasta_path <- "<path_to_reference_fasta>"
#' ref_gff <- parse_gff_file(gff3_path, ref_fasta_path)
#' }
#' @export
parse_gff_file = function(gff3_path, ref_fasta_path, perform_length_check = T){
  ref = read_ReferenceFasta(ref_fasta_path)
  gff = read_GFF3_Annotation(gff3_path)

  if(perform_length_check){
    # TODO: we can probably be more helpful than this!
    if(min(c(gff$start, gff$end)) < 0) stop("Invalid start position found!")
    if(max(c(gff$start, gff$end)) > ref$g) stop("Invalid stop position found!")
    if(any(gff$end < gff$start)) stop("Invalid start-stop pair found!")
  }

  return(list(gff = gff, ref = ref$ref, ref_name = ref$ref_name, g = ref$g, gff_path = gff3_path, ref_path = ref_fasta_path))

}
