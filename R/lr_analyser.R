#' analyse_long_range_links
#'
#' Function to analyse long range links. These links can be computed using BacGWES or Spydrpick (faster). Spydrpick is recommended for larger datasets.
#'
#' @importFrom Rfast2 Quantile
#'
#' @param dset name of the dataset, all outputs will be saved to the folder <dset>, possible to use the same dset name used for the short-range analysis
#' @param lr_links_path path to saved LR links, generally called <lr_links.tsv> or the SpydrPick links file
#' @param sr_links_path path to saved SR links, generally called <sr_links.tsv>
#' @param are_lrlinks_ordered are the links sorted in descending MI order (default = F), leave at default if unsure
#' @param SnpEff_Annotate specify whether to perform annotations using SnpEff (default = F)
#' @param snpeff_jar_path path to <snpEff.jar> (default = NULL). If annotations are not required, set SnpEff_Annotate = F
#' @param gbk_path path to genbank file, not needed if SnpEff_Annotate = F
#' @param snp.dat output from BacGWES::parse_fasta_alignment(), not needed if SnpEff_Annotate = F
#' @param cds_var output from BacGWES::estimate_variation_in_CDS(), not needed if SnpEff_Annotate = F
#' @param max_tophits specify the maximum number of long range links to save as <lr_tophits.tsv>. Note: all short-range links will be annotated (and saved separately),
#' but only the top <max_tophits> will be used for visualisation (default = 500)
#'
#' @examples
#' \dontrun{
#' analyse_long_range_links = function("dset", "dset/lr_links.tsv", "dset/sr_links.tsv")
#' }
#' @export
analyse_long_range_links = function(dset, lr_links_path, sr_links_path, are_lrlinks_ordered = F, SnpEff_Annotate = F,
                                    snpeff_jar_path = NULL, gbk_path = NULL, snp.dat = NULL, cds_var = NULL, max_tophits = 500){
                                    # tanglegram_break_segments = 5){

  # it makes sense to have a larger max_tophits for long range links - there will be a lot more of long-range links compared to short
  t_global = Sys.time()

  if(!file.exists(dset)) dir.create(dset) # folder to save

  len = MI = NULL # avoid ggplot NSE issue (https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/)

  # NOTE: if links are from spydrpick, ARACNE might already be there (DO NOT RE-RUN!)
  # NOTE: spydrpick does not add clusters, add them from paint (requires cds_var)

  if(SnpEff_Annotate == T) {
    if(is.null(snpeff_jar_path)) stop("You must specify <snpeff_jar_path> for annotations. To run without annotations, re-run analyse_long_range_links() with SnpEff_Annotate = F")
    if(!file.exists(snpeff_jar_path)) stop(paste("<SnpEff.jar> not found at:", snpeff_jar_path, "please check the path provided"))
    if(is.null(snpeff_jar_path)) stop("You must specify <snpeff_jar_path> for annotations. To run without annotations, re-run analyse_long_range_links() with SnpEff_Annotate = F")
    if(is.null(gbk_path)) stop("You must specify <gbk_path> for annotations. To run without annotations, re-run analyse_long_range_links() with SnpEff_Annotate = F")
    if(is.null(snp.dat)) stop("You must specify <gbk_path> for annotations. To run without annotations, re-run analyse_long_range_links() with SnpEff_Annotate = F")
    if(is.null(cds_var)) stop("You must specify <gbk_path> for annotations. To run without annotations, re-run analyse_long_range_links() with SnpEff_Annotate = F")
  }

  # lr_links_path = "~/Desktop/BacGWES_RUN/maela/lr_links.tsv"
  lr_links = BacGWES::readLongRangeLinks(lr_links_path)

  # sr_links_path = "~/Desktop/BacGWES_RUN/maela/sr_links.tsv"
  sr_links = readShortRangeLinks(sr_links_path)
  # sr_links = BacGWES::readShortRangeLinks(sr_links_path)

  q13 = Rfast2::Quantile(lr_links$MI, probs = c(0.25, 0.75)) # Global threshold

  # Can we make cluster based plots (doesn't seem to make much sense)
  clusts = sort(unique(c(lr_links$c1, lr_links$c2))) # these are numeric
  lr_links$clust = NA
  for(i in 1:length(clusts)){
    for(j in i:length(clusts)){
      lnks = which((lr_links$c1 == i & lr_links$c2 == j) | (lr_links$c2 == i & lr_links$c1 == j))
      lr_links$clust[lnks] = paste(i,j,sep="")
      q13 = rbind(q13,
                  Rfast2::Quantile(lr_links$MI[lnks], probs = c(0.25, 0.75))) # Global threshold
    }
  }
  thresholds = t(apply(q13, 1, function(x) x[2] + (x[2] - x[1])*c(1.5, 3)))

  lr_links_ARACNE_check = rbind(data.frame(pos1 = lr_links$pos1, pos2 = lr_links$pos2, MI = lr_links$MI),
                                data.frame(pos1 = sr_links$pos1, pos2 = sr_links$pos2, MI = sr_links$MI))

  lr_links_ARACNE_check = lr_links_ARACNE_check[lr_links_ARACNE_check$MI > min(thresholds), ]

  lr_links_red = lr_links[lr_links$MI > min(thresholds), ] # we only care about outliers

  lr_links_red$ARACNE = runARACNE(lr_links_red, lr_links_ARACNE_check)

  # sr_links are sorted (largest to smallest)
  if(!are_lrlinks_ordered){ # if they are not sorted, sort below
    lr_links_red = lr_links_red[order(lr_links_red$MI, decreasing = T), ]
    rownames(lr_links_red) = NULL
  }

  lr_plt_path = file.path(dset, "lr_gwes.png")

  lr_links_red = lr_links_red[rev(1:nrow(lr_links_red)), ]
  p1 = ggplot2::ggplot() + ggplot2::geom_point(data = lr_links_red[which(lr_links_red$ARACNE == F), ], ggplot2::aes(x = len, y = MI), col = "#C0C0C0") +
    ggplot2::geom_point(data = lr_links_red[which(lr_links_red$ARACNE == 1), ], ggplot2::aes(x = len, y = MI), col = "#0868ac") +
    ggplot2::geom_hline(yintercept = max(thresholds), col = "#db4325") +
    ggplot2::theme_light() #+
  # ggplot2::facet_wrap('~clust')

  ggplot2::ggsave(plot = p1, filename = lr_plt_path, width = 2200, height = 1200, units = "px")
  if(SnpEff_Annotate == F){
    cat(paste("\nDone in", round(difftime(Sys.time(), t_global, units = "mins"), 3), "m ** \n"))
  } else {
    tophits_path = file.path(dset, "lr_tophits.tsv")
    parsed_gbk_path = file.path(dset, "parsed_gbk.rds")
    if(file.exists(parsed_gbk_path)){
      gbk = readRDS(parsed_gbk_path )
    } else {
      cat("Reading the GBK file \n")
      gbk = BacGWES::parse_genbank_file(gbk_path = gbk_path, g = snp.dat$g, length_check = T) # will return 1 if fails
    }

    tophits = BacGWES::perform_snpEff_annotations(dset_name = dset, annotation_folder = file.path(getwd(), dset),
                                                  snpeff_jar = snpeff_jar_path, gbk = gbk, gbk_path = gbk_path,
                                                  cds_var = cds_var, links_df = lr_links_red, snp.dat = snp.dat,
                                                  tophits_path = tophits_path, max_tophits = max_tophits, links_type = "LR")

    # Tanglegram is difficult to read when plotted like this, best to avoid!
    # tanglegram_path = file.path(dset, "LR_Tanglegram")
    # if(!file.exists(tanglegram_path)) dir.create(tanglegram_path)
    # BacGWES::create_tanglegram(tophits = tophits, gbk = gbk, tanglegram_folder = tanglegram_path,
    #                            break_segments = tanglegram_break_segments, links_type = "LR")

    cat("\n")
    gwesexplorer_path = file.path(dset, "LR_GWESExplorer")
    if(!file.exists(gwesexplorer_path)) dir.create(gwesexplorer_path)
    BacGWES::write_output_for_gwes_explorer(snp.dat = snp.dat, tophits = tophits,
                                            gwes_explorer_folder = gwesexplorer_path, links_type = "LR")

    cat("\n")

    netplot_path = file.path(dset, "lr_network_plot.png")
    BacGWES::create_network(tophits = tophits, netplot_path = netplot_path, links_type = "LR")

  }
  cat(paste("\nDone in", round(difftime(Sys.time(), t_global, units = "mins"), 3), "m ** \n"))
  # return(p1)

}
