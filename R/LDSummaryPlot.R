#' genomewide_LDMap
#'
#' Function to generate the bird's-eye view of the genomewide LD structure (measured as weighted-MI)
#'
#' @importFrom MatrixExtra crossprod
#' @importFrom Matrix sparseMatrix t
#' @importFrom grDevices colorRampPalette png
#' @importFrom heatmap3 heatmap3
#'
#' @param lr_links_path path to saved <lr_links.tsv> or spydrpick edges (links) file
#' @param sr_links_path path to saved SR links, generally called <sr_links.tsv>
#' @param plot_save_path path to save the LD Map
#' @param reducer positive value used to reduce the LDMap size (default = NULL, will use length(positions)/1e3).
#' The default paramter should work well, reducing the parameter will increase resolution at the cost of resources.
#' @param plot_title title for the LD Map (default = NULL)
#' @param links_from_spydrpick are the links computed using spydrpick (default = F)
#' @param from draw the plot within the region starting <from> (default = NULL)
#' @param to draw the plot within the region ending in <to> (default = NULL)
#'
#' @examples
#' \dontrun{
#' genomewide_LDMap(lr_links_path, sr_links_path, plot_save_path)
#' }
#'
#' @export
genomewide_LDMap = function(lr_links_path, sr_links_path, plot_save_path, reducer = NULL, plot_title = NULL,
                            links_from_spydrpick = F, from = NULL, to = NULL){
  t0 = Sys.time()

  # sanity checks
  if(!is.null(reducer)){
    if(reducer < 0){
      warning("<reducer> for genomewide_LDMap should be >0, set to default")
      reducer = NULL
    }
  }
  # from & to
  if(!is.null(from) & is.null(to)) stop("If <from> is provided, <to> must be provided as well!")
  if(!is.null(to) & is.null(from)) stop("If <to> is provided, <from> must be provided as well!")

  if(is.null(from) & is.null(to)) {
    gwplot = T
    } else {
      gwplot = F
      if(to <= from) stop("<to> must be greater than <from>!")
      if(from < 0 | to < 0) stop("<from> and <to> must be positive values")
      from = round(from)
      to = round(to)
    }

  cat("Reading links... \n")
  lr_links = LDWeaver::read_LongRangeLinks(lr_links_path, links_from_spydrpick = links_from_spydrpick)
  sr_links = LDWeaver::read_ShortRangeLinks(sr_links_path)

  cat("Converting to SparseMatrix... \n")

  if(gwplot){ # genomewide plot
    pos_vec = sort(unique(c(lr_links$pos1, lr_links$pos2, sr_links$pos1, sr_links$pos2)))
    if(is.null(plot_title)) plot_title = "Genomewide LD plot"
  } else { # we need to filter
    pos_vec = sort(unique(c(lr_links$pos1, lr_links$pos2, sr_links$pos1, sr_links$pos2)))

    tot_pos = round(length(pos_vec)/100)
    pos_vec = pos_vec[pos_vec < to & pos_vec > from]
    if(length(pos_vec) < tot_pos) stop("Not enough SNPs in the <from>-<to> range to generate the plot, please choose different values")
    if(is.null(plot_title)) plot_title = "LD plot"

    lr_links = lr_links[ (lr_links$pos1 >= from & lr_links$pos1 <= to) & (lr_links$pos2 >= from & lr_links$pos2 <= to), ]
    sr_links = sr_links[ (sr_links$pos1 >= from & sr_links$pos1 <= to) & (sr_links$pos2 >= from & sr_links$pos2 <= to), ]
  }

  if(nrow(lr_links) > 0){
    lr_links$pos1_i = as.numeric(factor(lr_links$pos1, levels = pos_vec))
    lr_links$pos2_i = as.numeric(factor(lr_links$pos2, levels = pos_vec))
  }
  if(nrow(sr_links) > 0){
    sr_links$pos1_i = as.numeric(factor(sr_links$pos1, levels = pos_vec))
    sr_links$pos2_i = as.numeric(factor(sr_links$pos2, levels = pos_vec))
  }

  i = c(lr_links$pos1_i, lr_links$pos2_i, sr_links$pos1_i, sr_links$pos2_i)
  j = c(lr_links$pos2_i, lr_links$pos1_i, sr_links$pos2_i, sr_links$pos1_i)
  x = c(lr_links$MI, lr_links$MI, sr_links$MI, sr_links$MI)
  #
  sprs = Matrix::sparseMatrix(i = i, j = j, x = x,
                              dims = c(length(pos_vec), length(pos_vec)),
                              dimnames = list( as.character(pos_vec),  as.character(pos_vec) ))

  if(is.null(reducer)){
    reducer = round(length(pos_vec)/1e3)
  } else {
    reducer = round(reducer)
  }

  if(reducer > 1){
    # Matrix::isSymmetric(sprs)
    cat(paste("Reducing dimensionality from", length(pos_vec) ,"x", length(pos_vec), "-> "))
    x = .mat(length(pos_vec), reducer)
    cat(paste(ncol(x) ,"x", ncol(x), "...\n"))
    reduced_mat = MatrixExtra::crossprod(x, MatrixExtra::crossprod(sprs, x))
    htm = as(reduced_mat, 'matrix')/reducer^2
    nms = as.character(pos_vec[seq(from = 1, to = length(pos_vec), by = reducer-1)])
    colnames(htm) = rownames(htm) = nms[1:nrow(htm)]
  } else {
    cat(paste("No reduction performed, chosen segment might be too sparse!"))
    htm = as(sprs, 'matrix')
    colnames(htm) = rownames(htm) = as.character(pos_vec)
  }



  # spam::image(ldmx)

  # diag(htm) = NA

  htm = log10(htm + 1e-5)
  htm = .rescale01(htm)

  cat("Preparing plot...")

  col = grDevices::colorRampPalette(c("white", "#E1B9B4","#AE452C","#802418"))(2056)

  grDevices::png(plot_save_path, width = 5000, height = 5250, res = 600) # are we leaving this static?
  heatmap3::heatmap3(htm, symm = T, Rowv = NA, Colv = NA,
                     balanceColor = F, col = col, cexRow = 1, cexCol = 1,
                     legendfun = function() ._legendfun2(legend=plot_title,col=c("white"),
                                                       lwd = 0, cex = 1) )
  grDevices::dev.off()

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 3), "s \n"))
}

#' ._legendfun2
#' @param legend legend
#' @param lwd lwd
#' @param cex cex
#' @param col col
#' @param ... ...
._legendfun2 <- function (legend = c("Group A", "Group B"), lwd = 3, cex = 1.1,
                        col = c("red", "blue"), ...) {
  plot(0, xaxt = "n", bty = "n", yaxt = "n", type = "n", xlab = "",
       ylab = "")
  legend("bottomright", legend = legend, lwd = lwd, col = col,
         bty = "n", cex = cex, ...)
}

#' .rescale01
#'
#' @param v vector to rescale between 0-1
#'
#' @return rescaled vector
#'
#' @examples
#' \dontrun{
#' .rescale01(sample(10,25,10))
#' }
.rescale01 = function(v){
  mv = min(v)
  rn = max(v) - min(v)

  rsc = (v - mv)/rn
  return(rsc)
}

#' .mat
#'
#' @param n number of dimensions
#' @param r reduction factor
#'
#' @return kernel reduction
#'
#' @examples
#' \dontrun{
#' .mat(10, 5)
#' }
.mat <- function(n, r) {
  suppressWarnings(matrix(c(rep(1, r), rep(0, n)), n, n/r))
}

