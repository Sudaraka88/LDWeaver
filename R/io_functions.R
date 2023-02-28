#' readTopHits
#'
#' @param top_hits_path path to saved top_hits file
#'
#' @return Top hits as a data.frame
#'
#' @examples
#' \dontrun{
#' top_hits_path = "<path_to_tophits_file>"
#' tophits = readTopHits(top_hits_path)
#' }
#' @export
readTopHits = function(top_hits_path){
  tophits = read.table(top_hits_path, sep = "\t", header = T, quote = "", comment.char = "")
  return(tophits)
}
