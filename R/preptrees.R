#' view_tree
#'
#' Function to visualise allele distributions and metadata aligned below the tree. The tree file must be readable with ape::read.tree() and <tip.labels> must match the sequence names in the
#' provided fasta file(s). Besides the tree file, this function requires minimum three inputs: (1) a links file or data.frame, (2) a SNP fasta alignment and (3) a position <text> file with each row containing
#' the genomic positions of the fasta alignment. The links file is straightforward: it can be either one of <Tophits> or the <Annotated_links> files, both short and long range links
#' can be provided. A subset of manually chosen links can also be provided as a data.frame. The fasta/position files can be generated using LDWeaver::generate_Links_SNPS_fasta().
#' If metadata is available, they can be visualised along with the allele data in a separate panel.
#'
#' @importFrom ape read.tree
#' @importFrom phytools midpoint.root
#' @importFrom ggtree ggtree gheatmap
#' @importFrom ggplot2 ggsave scale_fill_viridis_d
#' @importFrom ggnewscale new_scale_fill
#'
#' @param tree_path Path to tree file (must be readable using ape::read.tree). tip.labels must exactly match the sequence names in the fasta and metadata files.
#'
#' @param links_df A data.frame comprising a subset of tophits links, if provided, links are not read from the files. NOTE: Each links must (at a minimum) contain the fields: <pos1> & <pos2>.
#' If alleles from the <ntop_links> are to be visualised, ensure that the links_df is sorted such that the most important link is at the top.
#' @param lr_tophits_path Path to the lr_tophits file (typically Tophits/lr_tophits.tsv).
#' @param lr_annotated_links_path Path to the lr_annotated_links file (typically Annotated_links/lr_links_annotated.tsv).
#' @param sr_tophits_path Path to the lr_tophits file (typically Tophits/sr_tophits.tsv).
#' @param sr_annotated_links_path Path to the sr_annotated_links file (typically Annotated_links/sr_links_annotated.tsv).
#'
#' @param fasta_path fasta file comprising the sites associated with links
#' @param pos_file_path position file containing fasta sites
#'
#' @param perform_midpoint_rooting Specify wheter the tree should be midpoint rooted (default = T)
#'
#' @param metadata_df data.frame with metadata. Must contain a column labelled 'id' that exactly matches the sequence names in the fasta file and tip.labels of the tree.
#' @param ntop_links Specify the number of links to plot from tophits file(s). ntop_links will be selected from each file if both LR and SR links should be visualised (default = 10).
#' @param from draw the SNPs that are within the region starting <from> and <to> below (default = NULL). SNPs that are also in links with SNPs from the chosen region will be plotted.
#' @param to If <from> is provided, <to> muts also be provided and <ntop_links> will be ignored (default = NULL).
#'
#' @param offset_metadata specify the separation between tree and metadata layer. These parameters should be modified to tune the plot dimensions as required. (default = NULL)
#' @param offset_alleles specify the speciation between tree and alleles layer (default = NULL)
#' @param width_metadata specify the width of the metadata layer as a proportion of the tree size (default = NULL)
#' @param width_alleles specify the width of the alleles layer as a proportion of the tree size (default = NULL)
#'
#' @param plot_save_path path to save the tree plot (default = NULL, show plot instead of saving)
#' @param plot_height if saving, specify the plot height in inches (default = 20)
#' @param plot_width if saving, specify the plot height in inches (default = 15). Modify these numbers to tune the plot dimensions as required.
#'
#' @examples
#' \dontrun{
#' genomewide_LDMap(lr_links_path, sr_links_path, plot_save_path)
#' }
#'
#' @export
view_tree = function(tree_path, perform_midpoint_rooting = T, metadata_df = NULL,
                     fasta_path, pos_file_path, links_df = NULL,
                     lr_tophits_path = NULL, lr_annotated_links_path = NULL, #lr_fasta_path = NULL, lr_pos_path = NULL,
                     sr_tophits_path = NULL, sr_annotated_links_path = NULL, #sr_fasta_path = NULL, sr_pos_path = NULL,
                     ntop_links = 10, from = NULL, to = NULL, offset_metadata= NULL, offset_alleles = NULL,
                     width_metadata = NULL, width_alleles = NULL, plot_save_path = NULL, plot_height = 20,
                     plot_width = 15){

  # Change the file to give an optional fasta + pos file genereated from LDWeaver::snpdat_to_fa()


  # sanity checks and IO

  ### tree_file
  tree_path = normalizePath(tree_path)
  tree = ape::read.tree(tree_path); if(perform_midpoint_rooting) tree = phytools::midpoint.root(tree)

  if(!is.null(links_df)){
    top_hits = links_df; rm(links_df)
  } else {
    # links file(s)
    if(!is.null(lr_tophits_path)) lrt = LDWeaver::read_TopHits(normalizePath(lr_tophits_path)) else lrt = NULL
    if(!is.null(lr_annotated_links_path)) lra = LDWeaver::read_TopHits(normalizePath(lr_annotated_links_path)) else lra = NULL
    if(!is.null(sr_tophits_path)) {
      srt = LDWeaver::read_TopHits(normalizePath(sr_tophits_path))
      rm_srp_idx = which(tolower(colnames(srt)) == "srp")
      if(length(rm_srp_idx) == 1) srt = srt[,-rm_srp_idx] else stop("sr_tophits file does not contain the srp column!")
    } else srt = NULL
    if(!is.null(sr_annotated_links_path)) {
      sra = LDWeaver::read_TopHits(normalizePath(sr_annotated_links_path))
      rm_srp_idx = which(tolower(colnames(sra)) == "srp")
      if(length(rm_srp_idx) == 1) sra = sra[,-rm_srp_idx] else stop("sr_annotated_links file does not contain the srp column!")
    } else sra = NULL

    srl = rbind(srt, sra); rm(sra, srt); if(!is.null(srl)) {srl = srl[!duplicated(srl), ]; srl = cbind(srl, link = "sr")}
    lrl = rbind(lrt, lra); rm(lra, lrt); if(!is.null(lrl)) {lrl = lrl[!duplicated(lrl), ]; lrl = cbind(lrl, link = "lr")}
    top_hits = rbind(srl, lrl); rm(srl, lrl)
  }
  # if(!is.null(lr_tophits_path) & !is.null(lr_annotated_links_path)) stop("Only one of the tophits or annotated links must be provided")
  # if(!is.null(sr_tophits_path) & !is.null(sr_annotated_links_path)) stop("Only one of the tophits or annotated links must be provided")

  # fasta/pos file
  pos = as.numeric(readLines(pos_file_path))
  fasta = read_fasta(fasta_path, pos, tree$tip.label)

  if(!is.null(metadata_df)){
    md_id_col = which("id"  == tolower(colnames(metadata_df)))
    if(length(md_id_col) != 1) stop("Metadata file must contain an ID column")
  }

  # if(!is.null(lr_tophits_path) | !is.null(lr_annotated_links_path)) {
  #   if(!is.null(lr_tophits_path)) lr_tophits_path = normalizePath(lr_tophits_path)
  #   if(!is.null(lr_annotated_links_path)) lr_annotated_links_path = normalizePath(lr_annotated_links_path)
  #   if(is.null(lr_fasta_path) | is.null(lr_pos_path)) stop("All 3 LR files must be provided to visualise")
  #   lr_fasta_path = normalizePath(lr_fasta_path)
  #   lr_pos_path = normalizePath(lr_pos_path)
  #   lr_vis = T
  # } else {
  #   lr_vis = F
  # }
  #
  # if(!is.null(sr_tophits_path) | !is.null(sr_annotated_links_path)) {
  #   if(!is.null(sr_tophits_path)) sr_tophits_path = normalizePath(sr_tophits_path)
  #   if(!is.null(sr_annotated_links_path)) sr_annotated_links_path = normalizePath(sr_annotated_links_path)
  #   if(is.null(sr_fasta_path) | is.null(sr_pos_path)) stop("All 3 SR files must be provided to visualise")
  #   sr_fasta_path = normalizePath(sr_fasta_path)
  #   sr_pos_path = normalizePath(sr_pos_path)
  #   sr_vis = T
  # } else {
  #   sr_vis = F
  # }

  # if(lr_vis & sr_vis) vis_lrsr = T else vis_lrsr = F
  #
  # if(lr_vis == F & sr_vis == F) stop("At least one set of link files (LR or SR) must be provided")

  # Use from-to positions
  if(!is.null(from)) if(is.null(to)) stop("<to> must also be provided")
  if(!is.null(to)) if(is.null(from)) stop("<from> must also be provided")
  if(!is.null(from) & !is.null(to)){
    if(to < from) stop("<from> must be less than <to>")
    if(from < 0) stop("<from> must be positive")
    if(to < 0) stop("<to> must be positive")
    from = round(from)
    to = round(to)
    ntop_links = NULL
  }

  # Use ntop_links
  if(!is.null(ntop_links)){
    from = to = NULL # this must already be the case
    if(ntop_links < 0) stop("<ntop_links> must be positive")
    if(ntop_links > 10) warning("Plot may be cluttered due to large <ntop_links> value")
  }


  # Read the tophits
  # if(vis_lrsr){ # both must be read
  #
  #   ## sr_links
  #   if(!is.null(sr_tophits_path)){
  #     sr_tophits = cbind(LDWeaver::read_TopHits(sr_tophits_path), link = "sr") # tophits file
  #   } else {
  #     sr_tophits = cbind(LDWeaver::read_AnnotatedLinks(sr_annotated_links_path), link = "sr") # annotated links file
  #   }
  #   rm_srp_idx = which(tolower(colnames(sr_tophits)) == "srp")
  #   if(length(rm_srp_idx) == 1) sr_tophits = sr_tophits[,-rm_srp_idx] else stop("sr_tophits file does not contain the srp column!")
  #
  #   ## lr_links
  #   if(!is.null(lr_tophits_path)){
  #     top_hits = rbind(sr_tophits,
  #                      cbind(LDWeaver::read_TopHits(lr_tophits_path), link = "lr"))
  #   } else {
  #     top_hits = rbind(sr_tophits,
  #                      cbind(LDWeaver::read_AnnotatedLinks(lr_annotated_links_path), link = "lr"))
  #   }
  #
  # } else if(lr_vis){
  #   if(!is.null(lr_tophits_path)){
  #     top_hits = cbind(LDWeaver::read_TopHits(lr_tophits_path), link = "lr")
  #   } else {
  #     top_hits = cbind(LDWeaver::read_AnnotatedLinks(lr_annotated_links_path), link = "lr")
  #   }
  #
  # } else if(sr_vis){
  #   if(!is.null(sr_tophits_path)){
  #     top_hits = cbind(LDWeaver::read_TopHits(sr_tophits_path), link = "lr")
  #   } else {
  #     top_hits = cbind(LDWeaver::read_AnnotatedLinks(sr_annotated_links_path), link = "lr")
  #   }
  # }

  # Read the fasta/pos files
  # if(vis_lrsr){
  #   pos_lr = as.numeric(readLines(lr_pos_path))
  #   pos_sr = as.numeric(readLines(sr_pos_path))
  #
  #   fasta_lr = read_fasta(lr_fasta_path, pos_lr, tree$tip.label)
  #   fasta_sr = read_fasta(sr_fasta_path, pos_sr, tree$tip.label)
  #
  #   dup = c()
  #   for(i in 1:length(pos_lr)) dup[i] = pos_lr[i] %in% pos_sr # duplicated columns in the lr dataset
  #   if(all(dup)){ # Everything in lr is already in sr
  #     pos = pos_sr
  #     fasta = fasta_sr
  #   } else if(any(dup)){ # some LR values are duplicated
  #     pos_lr = pos_lr[!dup]
  #     fasta_lr = fasta_lr[,!dup]
  #
  #     pos = c(pos_sr, pos_lr)
  #     fasta = cbind(fasta_sr, fasta_lr); fasta = fasta[, order(pos)]
  #     pos = sort(pos)
  #
  #     if(! all(as.numeric(colnames(fasta)) == pos) ) stop("Mismatch between fasta and position files")
  #
  #   } else { # no duplicates
  #
  #     pos = c(pos_sr, pos_lr)
  #     fasta = cbind(fasta_sr, fasta_lr); fasta = fasta[, order(pos)]
  #     pos = sort(pos)
  #
  #     if(! all(as.numeric(colnames(fasta)) == pos) ) stop("Mismatch between fasta and position files")
  #
  #   }
  # } else if(lr_vis){
  #   pos = as.numeric(readLines(lr_pos_path))
  #   fasta = read_fasta(lr_fasta_path, pos, tree$tip.label)
  # } else if(sr_vis){
  #   pos = as.numeric(readLines(sr_pos_path))
  #   fasta = read_fasta(sr_fasta_path, pos, tree$tip.label)
  # }

  # add SNP data from tophits to the plot
  snp_chosen_nums = c()
  if(!is.null(ntop_links)){
    if(!is.null(links_df)){ # links_df is provided, just take the top set of links as the most important
      snp_chosen_nums = c(snp_chosen_nums, 1:ntop_links)
    } else { # else we have the lr/sr info and can pick from both pools
      lr_l = which(top_hits$link == "lr")
      if(length(lr_l) > 0) snp_chosen_nums = c(snp_chosen_nums, lr_l[1:ntop_links])
      sr_l = which(top_hits$link == "sr")
      if(length(sr_l) > 0) snp_chosen_nums = c(snp_chosen_nums, sr_l[1:ntop_links])
    }
  }


  # print(snp_chosen_nums)

  if(!is.null(from) & !is.null(to)){
    snp_chosen_nums = unique(c(which(top_hits$pos1 >= from & top_hits$pos1 <= to),
                               which(top_hits$pos2 >= from & top_hits$pos2 <= to)))

  }

  # print(snp_chosen_nums)

  df7 = data.frame()
  if(length(snp_chosen_nums) > 0){
    pos_plot = sort(unique(c(matrix(c(top_hits$pos1[snp_chosen_nums], top_hits$pos2[snp_chosen_nums]), byrow = T, nrow = 2))))
    mask = rep(T, length(pos_plot))
    snp_pos2 = c()
    for(i in 1:length(pos_plot)){
      idx = which(pos_plot[i] == pos)
      if(length(idx) != 1) {
        mask[i] = F
        warning(paste(pos_plot[i], "not available in the provided fasta file(s)"))
        next
      }
      snp_pos2 = c(snp_pos2, idx)
    }

    df7 = data.frame(fasta[,snp_pos2])
    rownames(df7) = tree$tip.label
    if(any(mask)){
      pos_plot = pos_plot[mask]
      colnames(df7) = as.character(pos_plot)
    } else {
      df7 = data.frame()
    }


  }


  # add data from metadata to the plot
  df5 = data.frame()
  if(!is.null(metadata_df)){
    md_id = metadata_df[,md_id_col]
    md_ord_idx = c()
    for(i in 1:length(tree$tip.label)){
      tmp = which(tree$tip.label[i] == md_id)
      if(length(tmp) == 0) stop("Entry in tree$tip.label missing in <metadata_df> ids")
      md_ord_idx = c(md_ord_idx, tmp[1]) # multiple matches can be quietly allowed
    }

    metadata_df = metadata_df[md_ord_idx, ]
    if(nrow(metadata_df) > 0){
      df5 = data.frame(metadata_df[, -md_id_col])
      colnames(df5) = colnames(metadata_df)[-md_id_col]
      rownames(df5) = tree$tip.label
    }
  }
  p = ggtree::ggtree(tree, ladderize = T, layout = "dendrogram")

  if(is.null(width_metadata)) width_metadata = ncol(df7)/200
  if(is.null(offset_metadata)) offset_metadata = 0
  if(is.null(offset_alleles)) offset_alleles = ncol(df7)/100
  if(is.null(width_alleles)) width_alleles = 5

  cat(paste("Current_Values:\nwidth_metadata=",width_metadata, "\noffset_metadata=", offset_metadata, "\nwidth_alleles=", width_alleles, "\noffset_alleles=", offset_alleles,
            "\n"))

  if(nrow(df5) > 0){
    p = ggtree::gheatmap(p, df5, offset = offset_metadata, width = width_metadata, colnames_angle = 0, colnames_offset_y =0, hjust = 1) +
      ggplot2::scale_fill_viridis_d(option="C", name="Metadata")
    p = p + ggnewscale::new_scale_fill()
  }


  if(nrow(df7) > 0){
    p = ggtree::gheatmap(p, df7, offset = offset_alleles, width = width_alleles, colnames_angle = 0, colnames_offset_y =0, hjust = 1) +
      ggplot2::scale_fill_viridis_d(option="H", name="Alleles")
  }



  # p

  if(!is.null(plot_save_path)){
    ggplot2::ggsave(filename = plot_save_path, plot= p, width = plot_width, height = plot_height, units = 'in')
    p
  } else {
    p
  }


}

read_fasta = function(fasta_path, pos, order){
  fasta_out = .readFasta(fasta_path, length(pos))
  seqs = matrix(NA, ncol = fasta_out$seq.length, nrow = fasta_out$num.seqs)

  for(i in 1:fasta_out$num.seqs){
    seqs[i,] = unlist(strsplit(fasta_out$seq.s[i],""))
  }


  order_idx = c()
  for(i in order){
    idx = which(i == fasta_out$seq.names)
    if(length(idx) != 1) {
      stop("Sequence names mismatch between provided tree file and fasta file")
    }
    order_idx = c(order_idx, idx)
  }

  seqs = seqs[order_idx, ]
  rownames(seqs) = order
  colnames(seqs) = as.character(pos)
  return(seqs)
}

