#' read_TopHits
#'
#' @param top_hits_path path to saved top_hits file
#'
#' @return Top hits as a data.frame
#'
#' @examples
#' \dontrun{
#' top_hits_path = "<path_to_tophits_file>"
#' tophits = read_TopHits(top_hits_path)
#' }
#' @export
read_TopHits = function(top_hits_path){
  tophits = read.table(top_hits_path, sep = "\t", header = T, quote = "", comment.char = "")
  return(tophits)
}

#' read_LongRangeLinks
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
#' lr_links = read_LongRangeLinks(lr_links_path)
#' }
#' @export
read_LongRangeLinks = function(lr_links_path, links_from_spydrpick = F, sr_dist = 20000){
  if(!links_from_spydrpick){
    lr_links = read.table(lr_links_path, sep = "\t", header = F, quote = "", comment.char = "")
    colnames(lr_links) = c("pos1", "pos2", "c1", "c2", "len", "MI")
  } else {
    # stop("Spydrpick links not yet supported")
    lr_links = read.table(lr_links_path, sep = " ", header = F, quote = "", comment.char = "")
    if(ncol(lr_links) == 5) colnames(lr_links) = c("pos1", "pos2", "len", "ARACNE", "MI")
    if(ncol(lr_links) == 4) colnames(lr_links) = c("pos1", "pos2", "len", "MI")
  }

  drops = lr_links$len < sr_dist
  if(any(drops)) lr_links = lr_links[!drops, ]

  return(lr_links)
}

#' read_ShortRangeLinks
#'
#' @param sr_links_path path to saved sr_links.tsv file
#'
#' @return Long range links as a data.frame
#'
#' @examples
#' \dontrun{
#' sr_links_path = "<path_to_sr_links.tsv_file>"
#' sr_links = read_ShortRangeLinks(sr_links_path)
#' }
#' @export
read_ShortRangeLinks = function(sr_links_path){
  sr_links = read.table(sr_links_path, sep = "\t", header = F, quote = "", comment.char = "")
  colnames(sr_links) = c("clust_c", "pos1", "pos2", "clust1", "clust2", "len", "MI", "srp_max", "ARACNE")

  return(sr_links)
}

#' read_AnnotatedLinks
#'
#' @param annotated_links_path path to saved <sr> or <lr>_annotated_links.tsv file
#'
#' @return Annotated links as a data.frame
#'
#' @examples
#' \dontrun{
#' sr_annotated_links_path = "<path_to_annotated_sr_links.tsv_file>"
#' sr_links_ann = read_AnnotatedLinks(sr_annotated_links_path)
#' }
#' @export
read_AnnotatedLinks = function(annotated_links_path){
  annotated_links = read.table(annotated_links_path, sep = "\t", quote = "", comment.char = "", header = T)
  return(annotated_links)
}


#' runARACNE
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
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
  cat(paste("Running ARACNE on", nrow(links_to_check), "links... \n"))
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

    # comXZ = Rfast2::Intersect(matX, matZ) # either of matX or matZ should change for a new link, cannot skip!
    # comXZ = matZ[LDWeaver:::.fast_intersect(matX, matZ)+1] # replacement function
    comXZ = .fast_intersect(matX, matZ) # replacement function

    ##############
    # temporary check block to evaluate .fast_intersect accuracy
    # comXZ2 = intersect(matX, matZ)
    # if(length(comXZ) == 0 & length(comXZ) == 0){ # This is fine
    #
    # } else {
    #   if(length(comXZ) != length(comXZ2)){ # values chosen wrongly
    #     warning(paste("length mismatch: .fast_intersect", comXZ, "intersect", comXZ2))
    #   } else { # this is fine
    #     if(all(comXZ != sort(comXZ2))){
    #       warning(paste("val mismatch: .fast_intersect", comXZ, "intersect", comXZ2))
    #       cat("matX", matX)
    #       cat("matZ", matZ)
    #     }
    #   }
    # }
    ##############

    if(length(comXZ) > 0){ # These are the only links

      MI0X = MIs[idX[.vecPosMatch(comXZ, matX)], 1]
      MI0Z = MIs[idZ[.vecPosMatch(comXZ, matZ)], 1]
      ARACNE[i] = .compareTriplet(MI0X, MI0Z, MIs_chk[i,1])
    }
  }
  close(pb)
  cat(paste("\nDone in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  # table(ARACNE)
  return(ARACNE)
}

#' read_ReferenceFasta
#'
#' @param ref_fasta_path path to Reference fasta file. The file MUST be in fasta format and contain exactly one sequence!
#'
#' @return parsed reference sequence with metadata in list format
#'
#' @examples
#' \dontrun{
#' ref_fa = "<path_to_ref_fasta"
#' reference = read_ReferenceFasta(ref_fa)
#' }
read_ReferenceFasta = function(ref_fasta_path){
  # This function is usually super quick, a timer is unnecessary
  ref_fasta_path = normalizePath(ref_fasta_path)
  if(!file.exists(ref_fasta_path)) stop(paste(ref_fasta_path, "not found!"))
  parsedRef = .extractRef(ref_fasta_path)

  g = parsedRef$seq.length
  if(g <= 0) stop("empty sequence!")

  ref = unname(unlist(strsplit(parsedRef$seq, split = "", fixed = T)))
  if(length(ref) != g) stop("sequence length mismatch!")

  ref_name = parsedRef$seq.name
  if(is.null(ref_name)) stop("empty sequence!")

  return(list(ref = ref,
              ref_name = ref_name,
              g = g))
}

#' read_GFF3_Annotation
#'
#' @importFrom ape read.gff
#'
#' @param gff3_path path to gff3_annotation file. This function relies on ape::read.gff() to read in the file,
#' refer to ?ape::read.gff for additional help
#'
#' @return parsed annotation file in list format
#'
#' @examples
#' \dontrun{
#' gff3_path = "<path_to_gff3_annotation"
#' gff = read_GFF3_Annotation(gff3_path)
#' }
read_GFF3_Annotation = function(gff3_path){
  gff3_path = normalizePath(gff3_path)
  if(!file.exists(gff3_path)) stop(paste(gff3_path, "not found!"))

  gff = ape::read.gff(file = gff3_path)

  return(gff)
}

#' cleanup
#'
#' Function to clean up the files in <dset> folder, this will organise the saved files in to a more sensible folder structure.
#' If clean up was performed already, nothing will be changed.
#'
#' @param dset <dset> folder that requires cleanup.
#' @param delete_after_moving specify whether the organised files should be deleted from the original location (default = F).
#' If False, all organised files will be moved to a directory called OLD
#'
#'
#' @examples
#' \dontrun{
#' gff3_path = "<path_to_gff3_annotation"
#' gff = read_GFF3_Annotation(gff3_path)
#' }
#' @export
cleanup = function(dset, delete_after_moving = F){
  cat("Cleaning up...\n")
  dset = normalizePath(dset)
  if(!file.exists(dset)) stop(paste("Dataset:", dset , "not found!"))

  mv_success = c()
  files = dir(dset)

  ##### Additional Outputs #####
  # New fit data 20231102
  idx = grep("^c*[:0-999:]_fit_data.rds$", files)
  if(length(idx) > 0){
    fldr = file.path(dset, "Fit")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  idx = c(grep("cds_var.rds", files), grep("hdw.rds", files), grep("parsed_gbk.rds", files), grep("parsed_gff3.rds", files), grep("snp_ACGTN.rds", files))
  if(length(idx) > 0){
    fldr = file.path(dset, "Additional_Outputs")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  #### Fit folder ####
  idx = c(grep("^c*[:0-999:]_fit.png$", files),  grep("CDS_clustering.png", files))
  if(length(idx) > 0){
    fldr = file.path(dset, "Fit")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  #### Annotated_links folder ####
  idx = grep("*_links_annotated.tsv", files)
  if(length(idx) > 0){
    fldr = file.path(dset, "Annotated_links")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  #### GWESPlots ####
  idx = grep("*_gwes(.+)png", files)
  if(length(idx) > 0){
    fldr = file.path(dset, "GWESPlots")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  #### Tophits ####
  idx = c(grep("*_tophits.tsv", files), grep("*_network_plot.png", files))
  if(length(idx) > 0){
    fldr = file.path(dset, "Tophits")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  #### GWESExplorer folder ####
  idx = c(grep("*_GWESExplorer", files))
  if(length(idx) > 0){
    fldr = file.path(dset, "GWESExplorer")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }

  #### Temp folder ####
  idx = c(grep("snpEff", files), grep("*.vcf", files), grep("*annotations.tsv", files), grep("*_links.tsv", files), grep("LDW_run_*", files))
  if(length(idx) > 0){
    fldr = file.path(dset, "Temp")
    cleanup_support(files = file.path(dset, files[idx]), fldr)
    mv_success = c(mv_success, idx)
  }


  #### Relocate or delete moved files ####
  moved_idx = sort(unique(mv_success))
  if(length(moved_idx) > 0){
    if(!delete_after_moving){
      fldr = file.path(dset, "OLD")
      if(!file.exists(fldr)) dir.create(fldr)
      chk = file.copy(file.path(dset, files[moved_idx]), fldr, overwrite = T, recursive = T)
    }
    unlink(file.path(dset, files[moved_idx]), recursive = T)
  }
}

#' cleanup_support
#'
#' File operations for the cleanup function
#'
#' @param files files to copy <from>
#' @param fldr folder to copy <to>
#'
#' @return parsed annotation file in list format
#'
#' @examples
#' \dontrun{
#' gff3_path = "<path_to_gff3_annotation"
#' gff = read_GFF3_Annotation(gff3_path)
#' }
cleanup_support = function(files, fldr){
  if(!file.exists(fldr)) dir.create(fldr)
  chk = file.copy(files, fldr, overwrite = F, recursive = T)
  if(any(chk==FALSE)){ # failure to copy
    wrn_chk = files[!chk]
    for(x in wrn_chk) cat("Not overwriting:", x, "\n")
  }
}

#' snpdat_to_fa
#'
#' Create a SNP only fasta or tsv file from snp.dat. The user has the option to include only chosen positions.
#'
#' @param snp.dat output from parsing the multi fasta alignment using LDWeaver::parse_fasta_alignment(), generally saved to Additional_Outputs/snp_ACGTN.rds
#' @param aln_path path to save output fasta file or tsv file
#' @param pos_path path to write the SNP positions, only required for format = "fasta" (default = NULL)
#' @param pos set of loci to be included in the saved alignment (default = NULL - save everything). <pos> should be a numeric vector comprising a subset of sites in snp.dat$POS and cannot contain duplicates
#' @param format Output can be presented as a standard fasta file, or an optional tsv (default = "fasta")
#'
#' @export
snpdat_to_fa = function(snp.dat, aln_path, pos_path = NULL, pos = NULL, format = "fasta"){
  t0 = Sys.time()
  # sanity check
  if(format != "fasta" & format != "tsv") {
    warning(paste("Format", format, "unsupported, has to be: <fasta> or <tsv>. Changed to default <fasta>"))
    format = "fasta"
  }

  if(format == "fasta" & is.null(pos_path)) stop("Saving in fasta format requires a path for the pos file <pos_path>")

  snps_idx = c()
  if(is.null(pos)){ # use all positions
    snps_idx = 1:length(snp.dat$POS)
    pos = snp.dat$POS
  } else { # search through snp.dat$POS for pos value and extract index
    pos = sort(pos)
    if(any(duplicated(sort(pos)))) stop("Duplicated entries found in pos")
    duplicated(pos)
    for(i in 1:length(pos)){
      idx = which(snp.dat$POS %in% pos[i])
      if(length(idx) != 1) stop(paste("pos=", pos[i], "cannot be extracted from snp.dat"))
      snps_idx = c(snps_idx, idx)
    }
    # snps_idx = sapply(pos, function(x) which(snp.dat$POS %in% x))
  }

  cat("Converting... ")
  fasta = matrix(rep(NA, snp.dat$nseq*length(snps_idx)), nrow = length(snps_idx))
  tidx = Matrix::which(snp.dat$snp.matrix_A[snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "A"
  tidx = Matrix::which(snp.dat$snp.matrix_C[snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "C"
  tidx = Matrix::which(snp.dat$snp.matrix_G[snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "G"
  tidx = Matrix::which(snp.dat$snp.matrix_T[snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "T"
  tidx = Matrix::which(snp.dat$snp.matrix_N[snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "N"
  fasta = t(fasta)

  cat("Writing... ")
  # write fasta+pos file
  if(format == "fasta"){
    for(i in 1:snp.dat$nseq){
      write.table(paste(">", snp.dat$seq.names[i], sep = ""), aln_path, quote = F, col.names = F, row.names = F, append = T)
      write.table(paste(fasta[i,], collapse = ""), aln_path, quote = F, col.names = F, row.names = F, append = T)
    }
    write.table(x = pos, file = pos_path, quote = F, col.names = F, row.names = F)
  } else if(format == "tsv"){
    # write tsv file
    rownames(fasta) = snp.dat$seq.names
    colnames(fasta) = pos
    write.table(fasta, aln_path, sep = "\t", quote = F)

  } else {
    stop(paste("unknown save format", format)) # can't come here
  }
  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))

}

#' generate_Links_SNPS_fasta
#'
#' Create a link SNP only fasta file. This file is required to generate (detailed) tree plots
#'
#' @param snp.dat output from parsing the multi fasta alignment using LDWeaver::parse_fasta_alignment(), generally saved to Additional_Outputs/snp_ACGTN.rds
#' @param aln_path path to save output fasta file
#' @param pos_path path to write the SNP positions
#' @param lr_tophits_path path to the lr_tophits file (typically Tophits/lr_tophits.tsv).
#' @param lr_annotated_links_path path to the lr_annotated_links file (typically Annotated_links/lr_links_annotated.tsv).
#' @param sr_tophits_path path to the lr_tophits file (typically Tophits/sr_tophits.tsv).
#' @param sr_annotated_links_path path to the sr_annotated_links file (typically Annotated_links/sr_links_annotated.tsv).

#' @export
generate_Links_SNPS_fasta = function(snp.dat, aln_path, pos_path, lr_tophits_path = NULL, lr_annotated_links_path = NULL,
                                     sr_tophits_path = NULL, sr_annotated_links_path = NULL){

  if(is.null(lr_tophits_path) & is.null(lr_annotated_links_path) & is.null(sr_tophits_path) & is.null(sr_annotated_links_path)) stop("At least one links file must be provided")

  pos = c()
  if(!is.null(lr_tophits_path)) {
    temp = LDWeaver::read_TopHits(normalizePath(lr_tophits_path))
    pos = c(pos, temp$pos1, temp$pos2)
  }
  if(!is.null(sr_tophits_path)) {
    temp = LDWeaver::read_TopHits(normalizePath(sr_tophits_path))
    pos = c(pos, temp$pos1, temp$pos2)
  }
  if(!is.null(lr_annotated_links_path)) {
    temp = LDWeaver::read_AnnotatedLinks(normalizePath(lr_annotated_links_path))
    pos = c(pos, temp$pos1, temp$pos2)
  }
  if(!is.null(sr_annotated_links_path)) {
    temp = LDWeaver::read_AnnotatedLinks(normalizePath(sr_annotated_links_path))
    pos = c(pos, temp$pos1, temp$pos2)
  }

  pos = sort(pos); pos = pos[!duplicated(pos)]

  LDWeaver::snpdat_to_fa(snp.dat = snp.dat, aln_path = aln_path, pos_path = pos_path, pos = pos, format = "fasta")


}
