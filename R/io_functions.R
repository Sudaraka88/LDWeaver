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

#' readLongRangeLinks
#'
#' @param lr_links_path path to saved lr_links.tsv or spydrpick file
#' @param links_from_spydrpick are the links computed using spydrpick (default = F)
#' @param sr_dist links less than <sr_dist> apart will be dropped from the data.frame (default = 20000)
#'
#' @return Long range links as a data.frame
#'
#' @examples
#' \dontrun{
#' lr_links_path = "<path_to_lrlinks.tsv_file>"
#' lr_links = readLongRangeLinks(lr_links_path)
#' }
#' @export
readLongRangeLinks = function(lr_links_path, links_from_spydrpick = F, sr_dist = 20000){
  if(!links_from_spydrpick){
    lr_links = read.table(lr_links_path, sep = "\t", header = F, quote = "", comment.char = "")
    colnames(lr_links) = c("pos1", "pos2", "c1", "c2", "len", "MI")
  } else {
    stop("Spydrpick links not yet supported")
  }

  drops = lr_links$len < sr_dist
  if(any(drops)) lr_links = lr_links[!drops, ]

  return(lr_links)
}

#' readShortRangeLinks
#'
#' @param sr_links_path path to saved sr_links.tsv file
#'
#' @return Long range links as a data.frame
#'
#' @examples
#' \dontrun{
#' sr_links_path = "<path_to_sr_links.tsv_file>"
#' sr_links = readShortRangeLinks(sr_links_path)
#' }
#' @export
readShortRangeLinks = function(sr_links_path){
  sr_links = read.table(sr_links_path, sep = "\t", header = F, quote = "", comment.char = "")
  colnames(sr_links) = c("clust_c", "pos1", "pos2", "clust1", "clust2", "len", "MI", "srp_max", "ARACNE")

  return(sr_links)
}

#' runARACNE
#'
#' @param links_to_check data.frame comprising subset of links to run ARACNE on
#' @param links_full data.frame with all links
#' @return Long range links as a data.frame
#'
#' @examples
#' \dontrun{
#' sr_links_path = "<path_to_sr_links.tsv_file>"
#' sr_links = readShortRangeLinks(sr_links_path)
#' sr_links_red = sr_links[which(sr_links$srp > 3), ]
#' sr_links_red$ARACNE = runARACNE(sr_links_red, sr_links)
#' }
runARACNE = function(links_to_check, links_full){
  t0 = Sys.time()
  cat(paste("Running ARACNE\n"))
  pos_mat = matrix(c(links_full$pos1, links_full$pos2), nrow = nrow(links_full)) # for the reduced link set
  MIs = matrix(links_full$MI)

  ## mx of values to check
  nlinks = nrow(links_to_check)
  pos_mat_chk = matrix(c(links_to_check$pos1, links_to_check$pos2), nrow = nlinks)
  MIs_chk = matrix(links_to_check$MI)

  ARACNE = rep(T, nlinks) # unchekable links must be marked T
  t0 = Sys.time()
  pb = utils::txtProgressBar(min = 1, max = nlinks, initial = 1)

  pX_ = 0

  for(i in 1:nlinks){
    utils::setTxtProgressBar(pb,i)
    pX = pos_mat_chk[i,1]; # X = which(.compareToRow(POS, pX)) #which(POS %in% pX)
    pZ = pos_mat_chk[i,2]; #Z = which(.compareToRow(POS, pZ)) #which(POS %in% pZ)

    if(pX != pX_){ # skip this check if pX is unchanged
      idX = which(.compareToRow(pos_mat, pX)) #decodeIndex(index[[X]])
      matX = c(rbind(pos_mat[idX,1], pos_mat[idX,2])); matX = matX[matX != pX]
      pX_ = pX
    }
    idZ = which(.compareToRow(pos_mat, pZ)) #decodeIndex(index[[Z]])
    matZ = c(rbind(pos_mat[idZ,1], pos_mat[idZ,2])); matZ = matZ[matZ != pZ]

    comXZ = Rfast2::Intersect(matX, matZ) # either of matX or matZ should change for a new link, cannot skip!

    if(length(comXZ) > 0){ # These are the only links
      MI0X = MIs[idX[.vecPosMatch(comXZ, matX)], 1]
      MI0Z = MIs[idZ[.vecPosMatch(comXZ, matZ)], 1]
      ARACNE[i] = .compareTriplet(MI0X, MI0Z, MIs_chk[i,1])
    }
  }
  close(pb)
  cat(paste("\nDone in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  table(ARACNE)
  return(ARACNE)
}


