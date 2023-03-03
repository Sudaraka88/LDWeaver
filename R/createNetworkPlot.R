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
#' @param netplot_path folder to save the network
#' @param plot_title title for the network plot (default = NULL)
#' @param separator break pattern in annotation (default = ":", this is the default used by SnpEff)
#' @param max_plot_nodes specify the maximum number of nodes to visualise in the plot (default = NULL, use all). Large values will result in a cluttered output, consider increasing plot size.
#' @param plot_w specify plot width in pixels (default = 6000)
#' @param plot_h specify plot height in pixels (default = 4000)
#' @param links_type specify the links type long-range "LR" or short-range "SR" (default = "SR")
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' create_network(tophits, netplot_path)
#' }
#'
#' @export
create_network = function(tophits, netplot_path, plot_title = NULL, separator = ":",
                          max_plot_nodes = NULL, plot_w = 6000, plot_h = 4000, links_type = "SR"){
  cat("Preparing Network Plot ... ")
  Num_Links = name = NULL # avoid ggplot NSE issue (https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/)

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


    if(links_type == "SR"){
      df = data.frame(p1a, p2a, w = tophits$srp[1:avail_links])
    } else if(links_type == "LR") {
      df = data.frame(p1a, p2a, w = tophits$MI[1:avail_links])
    }
    df_uq =  plyr::ddply(df, .(p1a,p2a), nrow)

    if(nrow(df_uq) >= max_plot_nodes) c50=F
    avail_links = avail_links + 1
    if(avail_links > nrow(tophits)) c50=F
  }

  df_uq = df_uq[df_uq$V1 > 1, ]

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
    # s_weights = df_uq$w[kps]
    s_labs = df_uq$V1[kps]
  } else {
    stop("Everything is a loop!")
  }
  # s_weights = (s_weights)/max(s_weights)

  g1 = igraph::graph_from_edgelist(matrix(c(p1a_d, p2a_d), ncol = 2))
  # g1 = igraph::set.edge.attribute(g1, name = "weights", value = s_weights)
  g1 = igraph::set.edge.attribute(g1, name = "Num_Links", value = as.factor(s_labs))

  p1 = ggraph::ggraph(g1, layout = "nicely") +
    # ggraph::geom_node_point(size = 0) +
    ggraph::geom_edge_arc2(aes(label = Num_Links, colour = Num_Links), label_colour = NA, angle_calc = "along", force_flip = T,check_overlap = T) +
    ggraph::scale_edge_colour_discrete(guide = "legend") +
    # ggraph::geom_edge_arc2(aes(colour = labs)) +
    ggraph::geom_node_label(aes(label = name)) +
    ggplot2::theme_void() +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(legend.position = "bottom")

  # p1




    # g1 = intergraph::asNetwork(g1)
    # g1 = network::set.edge.attribute(x = g1, attrname = "weights", value = s_weights)
    # qs = stats::quantile(s_weights)
    # g1 = network::set.edge.attribute(x = g1, attrname = "color", value = ifelse(g1 %e% "weights" < qs[2], "grey75",
    #                                                                             ifelse(g1 %e% "weights" < qs[3], "grey50",
    #                                                                                    ifelse (g1 %e% "weights" < qs[4], "black", "red"))))
    #
    # g1 = network::set.edge.attribute(x = g1, attrname = "labs", value = s_labs)





  # net = ggnetwork::fortify(g1)
  #
  #
  # ggplot(net, aes(x, y, xend = xend, yend = yend)) +
  #   ggnetwork::geom_edges() +
  #   ggnetwork::geom_nodes(size = 1, color = "black")+
  #   ggnetwork::geom_edges(aes(color = weights)) +
  #   scale_color_gradient(low = "grey50", high = "tomato") +
  #   ggnetwork::geom_edgelabel(aes(label = labs)) +
  #   ggnetwork::geom_nodetext_repel(aes(label = vertex.names),show.legend = F)
  #
  # +
  #   ggnetwork::theme_blank()
  # ggplot(ggnetwork::ggnetwork(g1))

  # add weights


  # add colours



  # GGally::ggnet2(g1, label = T, edge.size = "weights", edge.color = "color", label.color = "blue4",
  #             node.color = "aliceblue", node.size = 6, shape = 13, edge.label.color = "darkorchid1", edge.label = "labs",
  #             label.size = 4, mode = "fruchtermanreingold")
  #
  # + ggrepel::geom_text_repel(aes(label = g1$vertex.names))
  #
  #
  # ?GGally::ggnet2()
  #
  # p1
  #
  # + ggrepel::geom_text_repel()
  #
  # p1

  ggplot2::ggsave(plot = p1, filename = netplot_path, width = plot_w, height = plot_h, units = "px", device = "png", dpi = 300, bg = "white")
  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))

}
