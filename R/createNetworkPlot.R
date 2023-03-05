#' create_network
#'
#' Function to generate a network view of top GWES hits
#'
#' @importFrom plyr . ddply
#' @importFrom stats quantile
#' @importFrom igraph graph_from_edgelist set.edge.attribute
#' @importFrom ggraph ggraph geom_edge_arc2 scale_edge_colour_discrete geom_node_label
#' @importFrom ggplot2 theme_void theme ggsave
#'
#' @param tophits data frame with top short or long range GWES links, returned from BacGWES::perform_snpEff_annotations()
#' @param netplot_path folder to save the network (default = NULL, return ggplot object)
#' @param plot_title title for the network plot (default = NULL)
#' @param separator break pattern in annotation (default = ":", this is the default used by SnpEff)
#' @param max_plot_nodes specify the maximum number of nodes to visualise in the plot (default = NULL, use all). Large values will result in a cluttered output, consider increasing plot size.
#' @param plot_w specify plot width in pixels (default = 6000)
#' @param plot_h specify plot height in pixels (default = 4000)
#' @param min_links_to_include specify the minimum number of links required to retain the node in the plot (default = 2)
#'
#' @return Network plots (ggplot object)
#'
#' @examples
#' \dontrun{
#' create_network(tophits, netplot_path)
#' }
#'
#' @export
create_network = function(tophits, netplot_path = NULL, plot_title = NULL, separator = ":",
                          max_plot_nodes = NULL, plot_w = 6000, plot_h = 4000,  min_links_to_include = 2){ #,links_type = "SR",
  # ){
  cat("Preparing Network Plot ... ")
  Num_Links = name = weights = NULL # avoid ggplot NSE issue (https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/)

  t0 = Sys.time()

  if(is.null(plot_title)) plot_title = ""

  # Let's ensure 50 nodes are present here
  c50 = T
  if(is.null(max_plot_nodes)){
    avail_links = nrow(tophits)
    max_plot_nodes = nrow(tophits)
  } else {
    avail_links = max_plot_nodes
  }



  while (c50){

    p1a = unname(sapply(tophits$pos1_ann[1:avail_links], function(x) unlist(strsplit(x, separator))[1]))
    p2a = unname(sapply(tophits$pos2_ann[1:avail_links], function(x) unlist(strsplit(x, separator))[1]))


    # if(links_type == "SR"){
    #   df = data.frame(p1a, p2a, w = tophits$srp[1:avail_links])
    # } else if(links_type == "LR") {
    df = data.frame(p1a, p2a, w = tophits$MI[1:avail_links])
    # }
    df_uq =  plyr::ddply(df, .(p1a,p2a), nrow)

    if(nrow(df_uq) >= max_plot_nodes) c50=F
    avail_links = avail_links + 1
    if(avail_links > nrow(tophits)) c50=F
  }

  df_uq = df_uq[df_uq$V1 >= min_links_to_include, ]

  # get the max srp/MI val. for each link (this is now redundant!)
  for(i in 1:nrow(df_uq)){
    df_uq$w[i] = max(df$w[which(df_uq$p1a[i] == df$p1a & df_uq$p2a[i] == df$p2a)])
  }

  # swapping positions might improve the way the plot looks
  swp_lr = which(sapply(df_uq$p1a, function(x) x %in% df_uq$p2a))
  if(length(swp_lr) > 0){
    tmp = df_uq$p1a[swp_lr]
    df_uq$p1a[swp_lr] = df_uq$p2a[swp_lr]
    df_uq$p2a[swp_lr] = tmp
  }

  swp_rl = which(sapply(df_uq$p2a, function(x) x %in% df_uq$p1a))
  if(length(swp_rl) > 0){
    tmp = df_uq$p1a[swp_rl]
    df_uq$p1a[swp_rl] = df_uq$p2a[swp_rl]
    df_uq$p2a[swp_rl] = tmp
  }


  # some links might be the same, merge them
  pst1 = paste(df_uq$p1a, df_uq$p2a)
  pst2 = paste(df_uq$p2a, df_uq$p1a)

  lnks_to_merge = which(sapply(pst1, function(x) x%in% pst2))
  if(length(lnks_to_merge) > 0){
    for(i in lnks_to_merge){
      mrg = unlist(strsplit(pst1[i], ' '))
      l1 = which(mrg[1] == df_uq$p1a & mrg[2] == df_uq$p2a)
      l2 = which(mrg[2] == df_uq$p1a & mrg[1] == df_uq$p2a)
      if(length(l1) == 1 & length(l2) == 1){ # the link to merge exists!
        df_uq$V1[l1] = sum(df_uq$V1[c(l1,l2)])
        df_uq$w[l1] = max(df_uq$w[c(l1,l2)])
        df_uq = df_uq[-l2, ]
      }
    }
  }

  # also need to drop loops
  kps = which(df_uq$p1a != df_uq$p2a)
  if(length(kps) > 0){
    p1a_d = df_uq$p1a[kps]
    p2a_d = df_uq$p2a[kps]
    s_weights = df_uq$w[kps]
    s_labs = df_uq$V1[kps]
  } else {
    stop("Everything is a loop!")
  }
  s_weights = (s_weights)/max(s_weights)

  g1 = igraph::graph_from_edgelist(matrix(c(p1a_d, p2a_d), ncol = 2))
  g1 = igraph::set.edge.attribute(g1, name = "weights", value = s_weights)
  g1 = igraph::set.edge.attribute(g1, name = "Num_Links", value = as.factor(s_labs))

  p1 = ggraph::ggraph(g1, layout = "nicely") +
    # ggraph::geom_node_point(size = 0) +
    ggraph::geom_edge_arc2(ggplot2::aes(label = Num_Links, colour = Num_Links, width = weights, alpha = weights),
                           label_colour = NA, angle_calc = "along", force_flip = T,check_overlap = T, position = "jitter", label_alpha = 1) +
    ggraph::scale_edge_colour_discrete(guide = "legend") +
    # ggraph::geom_edge_arc2(aes(colour = labs)) +
    ggraph::geom_node_label(ggplot2::aes(label = name)) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(legend.position = "bottom")

  if(is.null(netplot_path)){
    return(p1)
  } else {
    ggplot2::ggsave(plot = p1, filename = netplot_path, width = plot_w, height = plot_h, units = "px", device = "png", dpi = 300, bg = "white")
    return(p1)
  }

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))

}


#' create_network_for_gene
#'
#' Function to generate a network view of GWES hits for a given gene
#'
#' @param gene_name case sensitive name of genome region (must match the SnpEff annotated name)
#' @param sr_annotated_path path to short-range annotated links tsv file (default = NULL)
#' @param lr_annotated_path path to long-range annotated links tsv file (default = NULL)
#' @param drop_syXsy specify whether to drop synonymous-synonymous links (default = T)
#' @param drop_indirect specify whether to drop ARACNE indirect links (default = T)
#' @param level set 1 for only nearest neighbours, set 2 for an additional layer (default = 1)
#' @param separator break pattern in annotation (default = ":", this is the default used by SnpEff)
#' @param min_links_to_include specify the minimum number of links required to retain the node in the plot (default = 2)
#'
#'
#' @return data_frame to use as input for BacGWES::create_network()
#'
#' @examples
#' \dontrun{
#' create_network_for_gene(tophits, netplot_path)
#' }
#'
#' @export
create_network_for_gene = function(gene_name, sr_annotated_path = NULL, lr_annotated_path = NULL,
                                   drop_syXsy = T, drop_indirect = T, level = 1, separator = ":",
                                   min_links_to_include = 3) {


  if(is.null(sr_annotated_path) & is.null(sr_annotated_path)) stop("<sr> or <lr> annotated_link tsv file path must be provided!")
  if(level != 1) if(level != 2) stop('Level must be 1 or 2')

  tophits_df_sr = c()
  if(!is.null(sr_annotated_path)){
    sr_ann = BacGWES::read_AnnotatedLinks(sr_annotated_path)
    sr_idx = sort(unique(c(grep(gene_name, sr_ann$pos1_ann), grep(gene_name, sr_ann$pos2_ann))))
    if(length(sr_idx) > 0) {
      sr_links = sr_ann[sr_idx, ]
      tophits_df_sr = data.frame(pos1 = sr_links$pos1,
                                 pos2 = sr_links$pos2,
                                 pos1_ann = sr_links$pos1_ann,
                                 pos2_ann = sr_links$pos2_ann,
                                 MI = sr_links$MI,
                                 links = sr_links$links,
                                 ARACNE = sr_links$ARACNE)
    }
  }

  tophits_df_lr = c()
  if(!is.null(lr_annotated_path)){
    lr_ann = BacGWES::read_AnnotatedLinks(lr_annotated_path)
    lr_idx = sort(unique(c(grep(gene_name, lr_ann$pos1_ann), grep(gene_name, lr_ann$pos2_ann))))
    if(length(lr_idx) > 0) {
      lr_links = lr_ann[lr_idx, ]
      tophits_df_lr = data.frame(pos1 = lr_links$pos1,
                                 pos2 = lr_links$pos2,
                                 pos1_ann = lr_links$pos1_ann,
                                 pos2_ann = lr_links$pos2_ann,
                                 MI = lr_links$MI,
                                 links = lr_links$links,
                                 ARACNE = lr_links$ARACNE)
    }
  }

  tophits_df = rbind(tophits_df_sr, tophits_df_lr)

  # drop syXsy
  if(drop_syXsy) tophits_df = tophits_df[-which(tophits_df$links == "syXsy"), ]
  if(drop_indirect) tophits_df = tophits_df[tophits_df$ARACNE == 1, ]


  if(level == 2){ # we need to go through the pipeline again for each new gene
    c50 = T
    avail_links = nrow(tophits_df)
    max_plot_nodes = nrow(tophits_df)

    while (c50){

      p1a = unname(sapply(tophits_df$pos1_ann[1:avail_links], function(x) unlist(strsplit(x, separator))[1]))
      p2a = unname(sapply(tophits_df$pos2_ann[1:avail_links], function(x) unlist(strsplit(x, separator))[1]))
      df = data.frame(p1a, p2a, w = tophits_df$MI[1:avail_links])
      df_uq =  plyr::ddply(df, .(p1a,p2a), nrow)

      if(nrow(df_uq) >= max_plot_nodes) c50=F
      avail_links = avail_links + 1
      if(avail_links > nrow(tophits_df)) c50=F
    }

    df_uq = df_uq[df_uq$V1 >= min_links_to_include, ] # maybe this should be retained here (can be noisy for some cases, decide later)

    genes = unique(c(df_uq$p1a, df_uq$p2a))
    genes = genes[-which(genes == gene_name)]

    for(gene in genes){
      tophits_df_sr = c()
      if(!is.null(sr_annotated_path)){
        sr_idx = sort(unique(c(grep(gene, sr_ann$pos1_ann), grep(gene, sr_ann$pos2_ann))))
        if(length(sr_idx) > 0) {
          sr_links = sr_ann[sr_idx, ]
          tophits_df_sr = data.frame(pos1 = sr_links$pos1,
                                     pos2 = sr_links$pos2,
                                     pos1_ann = sr_links$pos1_ann,
                                     pos2_ann = sr_links$pos2_ann,
                                     MI = sr_links$MI,
                                     links = sr_links$links,
                                     ARACNE = sr_links$ARACNE)
        }
      }

      tophits_df_lr = c()
      if(!is.null(lr_annotated_path)){
        lr_idx = sort(unique(c(grep(gene, lr_ann$pos1_ann), grep(gene, lr_ann$pos2_ann))))
        if(length(lr_idx) > 0) {
          lr_links = lr_ann[lr_idx, ]
          tophits_df_lr = data.frame(pos1 = lr_links$pos1,
                                     pos2 = lr_links$pos2,
                                     pos1_ann = lr_links$pos1_ann,
                                     pos2_ann = lr_links$pos2_ann,
                                     MI = lr_links$MI,
                                     links = lr_links$links,
                                     ARACNE = lr_links$ARACNE)
        }
      }

      tophits_df = rbind(tophits_df, tophits_df_sr, tophits_df_lr)
    }
    if(drop_syXsy) tophits_df = tophits_df[-which(tophits_df$links == "syXsy"), ]
    if(drop_indirect) tophits_df = tophits_df[tophits_df$ARACNE == 1, ]


  }
  dups = duplicated(tophits_df)
  if(any(dups)) tophits_df = tophits_df[!dups, ]

  return(tophits_df)
  # BacGWES::create_network(tophits_df, links_type = "LR", min_links_to_include = 5)

}
