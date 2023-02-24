#' create_network
#'
#' Function to generate a network view of top GWES hits
#'
#' @importFrom htmlwidgets saveWidget
#' @importFrom stats cutree hclust dist
#' @importFrom genbankr cds
#' @importFrom plyr . ddply
#' @importFrom chromoMap chromoMap
#'
#' @param srlinks_tophits data frame with top short range GWES links, returned from BacGWES::perform_snpEff_annotations()
#' @param netplot_path folder to save the network
#' @param separator break pattern for annotation (default = ":", which is the default used by SnpEff)
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' create_network(srlinks_tophits, netplot_path)
#' }
#'
#' @export
create_network = function(srlinks_tophits, netplot_path, separator = ":", max_plot_nodes = 50){
  # Let's ensure 50 nodes are present here
  c50 = T
  avail_links = max_plot_nodes

  while (c50){

    p1a = unname(sapply(srlinks_tophits$pos1_ann[1:avail_links], function(x) unlist(strsplit(x, separator))[1]))
    p2a = unname(sapply(srlinks_tophits$pos2_ann[1:avail_links], function(x) unlist(strsplit(x, separator))[1]))


    df = data.frame(p1a, p2a, w = srlinks_tophits$srp[1:avail_links])
    df_uq =  plyr::ddply(df, .(p1a,p2a), nrow)

    if(nrow(df_uq) >= max_plot_nodes) c50=F
    avail_links = avail_links + 1
  }
  # get the max srp for each link
  for(i in 1:nrow(df_uq)){
    df_uq$w[i] = max(df$w[which(df_uq$p1a[i] == df$p1a & df_uq$p2a[i] == df$p2a)])
  }

  # swapping positions might improve the way the plot looks
  swp_lr = which(sapply(df_uq$p1a, function(x) x %in% df_uq$p2a))
  tmp = df_uq$p1a[swp_lr]
  df_uq$p1a[swp_lr] = df_uq$p2a[swp_lr]
  df_uq$p2a[swp_lr] = tmp

  swp_rl = which(sapply(df_uq$p2a, function(x) x %in% df_uq$p1a))
  tmp = df_uq$p1a[swp_rl]
  df_uq$p1a[swp_rl] = df_uq$p2a[swp_rl]
  df_uq$p2a[swp_rl] = tmp

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
    print("Everything is a loop!")
  }
  s_weights = (s_weights)/max(s_weights)

  g1 = igraph::graph_from_edgelist(matrix(c(p1a_d, p2a_d), ncol = 2))
  g1 = intergraph::asNetwork(g1)
  # add weights
  g1 = network::set.edge.attribute(x = g1, attrname = "weights", value = s_weights)

  # add colours
  qs = quantile(s_weights)
  g1 = network::set.edge.attribute(x = g1, attrname = "color", value = ifelse(g1 %e% "weights" < qs[2], "grey75",
                                                                              ifelse(g1 %e% "weights" < qs[3], "grey50",
                                                                                     ifelse (g1 %e% "weights" < qs[4], "black", "red"))))

  g1 = network::set.edge.attribute(x = g1, attrname = "labs", value = s_labs)

  p1 = ggnet2(g1, label = T, edge.size = "weights", edge.color = "color", label.color = "blue4",
              node.color = "aliceblue", node.size = 6, shape = 13, edge.label = "labs", edge.label.color = "darkorchid1",
              label.size = 4, mode = "fruchtermanreingold")

  # p1
  ggsave(plot = p1, filename = netplot_path, width = 6000, height = 4000, units = "px", device = "png", dpi = 300)

}
