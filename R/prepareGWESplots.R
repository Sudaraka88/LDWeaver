#' make_gwes_plots
#'
#' Function to generate long and short range GWES plots
#'
#' @importFrom ggplot2 ggplot aes geom_point ggtitle ggsave scale_colour_gradientn facet_wrap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils read.table
#'
#' @param lr_links data.frame containing lr_links or the path to saved lr_links TSV file.
#' Usually output from perform_MI_computation(), (default = NULL). Instead of using this function to generate the lr-gwes plot, use the dedicated
#' BacGWES::analyse_long_range_links() function for long range link analysis (and plotting).
#' @param sr_links data.frame containing sr_links or the path to saved sr_links TSV file.
#' Usually output from perform_MI_computation(), (default = NULL)
#' @param plt_folder specify the folder to save generated plots (default = NULL, will be saved to a folder called PLOTS in getwd())
#' @param are_srlinks_ordered if False, links will be ordered before plotting, which will plot higher srp links on top of lower ones (default = F).
#' Set True only if the links are already in the preferred order. Note: This function will reverse the link order to match the ggplot's top to bottom plotting order.
#'
#' @examples
#' \dontrun{
#' sr_links_red = perform_MI_computation(snp.dat = snp.dat, hdw = hdw, cds_var = cds_var, ncores = 10,
#' lr_save_path = "lr.tsv", plt_folder = "plots")
#' }
#'
#' @export
make_gwes_plots = function(lr_links=NULL, sr_links = NULL, plt_folder = NULL, are_srlinks_ordered = F){
  cat("Generating short and long range GWES plots ... ")
  t0 = Sys.time()
  len = MI = srp_max = NULL # avoid ggplot NSE issue (https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/)

  if(is.null(plt_folder)) plt_folder = file.path(getwd(), "PLOTS");

  if(!file.exists(plt_folder)) dir.create(plt_folder)

  # sanity checks
  lr_links_ready = T
  if(is.null(lr_links)) lr_links_ready = F

  if(lr_links_ready){
    if(!is.data.frame(lr_links)){ #lr_links can be a path to the TSV or a data.frame
      # lr_links must be a path then
      if(file.exists(lr_links)){
        lr_links = read.table(lr_links, sep = '\t', header = F, quote = "", comment.char = "")
      } else { # nothing at the given path
        lr_links_ready = F
      }
    }
    # print(head(lr_links))

    if(ncol(lr_links) == 6){ # force this order?
      colnames(lr_links) = c("pos1", "pos2", "c1", "c2", "len", "MI")
    } else {
      lr_links_ready = F
    }

    if(!lr_links_ready){
      stop("lr_links must either be (1) a data.frame with lr_links or (2) the path to the saved tsv file from perform_MI_computation()")
      # return(-1)
    }

    plr = ggplot2::ggplot(data = lr_links, ggplot2::aes(x = len, y = MI)) +
      ggplot2::geom_point()

    # save path
    lr_plt_path = file.path(plt_folder, "lr_gwes.png")
    ggplot2::ggsave(plot = plr, filename = lr_plt_path, width = 4800, height = 1200, units = "px")
  }


  # sanity checks
  sr_links_ready= T
  if(is.null(sr_links)) sr_links_ready = F

  if(sr_links_ready){
    if(!is.data.frame(sr_links)){ #sr_links can be a path to the TSV or a data.frame
      # sr_links must be a path then
      if(file.exists(sr_links)){
        sr_links = read.table(sr_links, sep = '\t', header = F, quote = "", comment.char = "")
      } else { # nothing at the given path
        sr_links_ready = F
      }
    }

    if(ncol(sr_links) == 9){ # force this order?
      colnames(sr_links) = c("clust_c", "pos1", "pos2", "clust1", "clust2", "len", "MI", "srp_max", "ARACNE")
    } else {
      sr_links_ready = F
    }

    if(!sr_links_ready){
      stop("sr_links must either be (1) a data.frame with sr_links or (2) the path to the saved tsv file from perform_MI_computation()")
      # return(-1)
    }


    # sr_links are sorted (largest to smallest)
    if(!are_srlinks_ordered){ # if they are not sorted, sort below
      sr_links = sr_links[order(sr_links$srp_max, decreasing = T), ]
      rownames(sr_links) = NULL
    }

    # Note: ggplot plots from top to bottom of DF, we need to reverse order the links so that the top links appear on top!
    sr_links = sr_links[rev(1:nrow(sr_links)), ]
    p1 = ggplot2::ggplot() +
      ggplot2::geom_point(data = sr_links[which(sr_links$ARACNE == 0), ], ggplot2::aes(x = len, y = MI), col = "#C0C0C0") +
      ggplot2::geom_point(data = sr_links[which(sr_links$ARACNE == 1), ], ggplot2::aes(x = len, y = MI, col = srp_max)) +
      ggplot2::scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(6, "RdYlBu"))) +
      ggplot2::facet_wrap('~clust_c') +
      ggplot2::theme_light()

    sr_segregated_plt_path = file.path(plt_folder, "sr_gwes_clust.png")
    ggplot2::ggsave(plot = p1, filename = sr_segregated_plt_path, width = 2200, height = 1200, units = "px")

    p2 = ggplot2::ggplot() +
      ggplot2::geom_point(data = sr_links[which(sr_links$ARACNE == 0), ], ggplot2::aes(x = len, y = MI), col = "#C0C0C0") +
      ggplot2::geom_point(data = sr_links[which(sr_links$ARACNE == 1), ], ggplot2::aes(x = len, y = MI, col = srp_max)) +
      ggplot2::scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(6, "RdYlBu"))) +
      ggplot2::theme_light()

    sr_combined_plt_path = file.path(plt_folder, "sr_gwes_combi.png")
    ggplot2::ggsave(plot = p2, filename = sr_combined_plt_path, width = 2200, height = 1200, units = "px")
  }
  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
}
