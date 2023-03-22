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
#' @param plot_title title for the LD Map (default = NULL)
#' @param reducer positive value used to reduce the LDMap size (default = NULL, will use length(positions)/1e3).
#' The default paramter should work well, reducing the parameter will increase resolution at the cost of resources.
#' @param links_from_spydrpick are the links computed using spydrpick (default = F)
#'
#' @examples
#' \dontrun{
#' genomewide_LDMap(lr_links_path, sr_links_path, plot_save_path)
#' }
#'
#' @export
genomewide_LDMap = function(lr_links_path, sr_links_path, plot_save_path, plot_title = NULL,
                            reducer = NULL, links_from_spydrpick = F){
  t0 = Sys.time()

  if(!is.null(reducer)){
    if(reducer < 0){
      warning("<reducer> for genomewide_LDMap should be >0, set to default")
      reducer = NULL
    }
  }

  if(is.null(plot_title)) plot_title = "Genomewide LD plot"

  cat("Reading links... \n")
  lr_links = LDWeaver::read_LongRangeLinks(lr_links_path, links_from_spydrpick = links_from_spydrpick)
  sr_links = LDWeaver::read_ShortRangeLinks(sr_links_path)

  cat("Converting to SparseMatrix... \n")
  pos_vec = sort(unique(c(lr_links$pos1, lr_links$pos2, sr_links$pos1, sr_links$pos2)))
  lr_links$pos1_i = as.numeric(factor(lr_links$pos1, levels = pos_vec))
  lr_links$pos2_i = as.numeric(factor(lr_links$pos2, levels = pos_vec))

  sr_links$pos1_i = as.numeric(factor(sr_links$pos1, levels = pos_vec))
  sr_links$pos2_i = as.numeric(factor(sr_links$pos2, levels = pos_vec))

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

  # Matrix::isSymmetric(sprs)
  cat(paste("Reducing dimensionality from", length(pos_vec) ,"x", length(pos_vec), "-> "))
  x = .mat(length(pos_vec), reducer)
  cat(paste(ncol(x) ,"x", ncol(x), "...\n"))
  reduced_mat = MatrixExtra::crossprod(x, MatrixExtra::crossprod(sprs, x))

  # spam::image(ldmx)
  htm = as(reduced_mat, 'matrix')/reducer^2
  # diag(htm) = NA

  nms = as.character(pos_vec[seq(from = 1, to = length(pos_vec), by = reducer-1)])

  colnames(htm) = rownames(htm) = nms[1:nrow(htm)]

  htm = log10(htm + 1e-5)
  htm = .rescale01(htm)

  cat("Preparing plot...")

  col = grDevices::colorRampPalette(c("white", "#E1B9B4","#AE452C","#802418"))(2056)

  grDevices::png(plot_save_path, width = 5000, height = 5250, res = 600) # are we leaving this static?
  heatmap3::heatmap3(htm, symm = T, Rowv = NA, Colv = NA,
                     balanceColor = F, col = col, cexRow = 1, cexCol = 1,
                     legendfun = function() legendfun2(legend=plot_title,col=c("white"),
                                                       lwd = 0, cex = 1) )
  grDevices::dev.off()

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 3), "s \n"))
}

legendfun2 <- function (legend = c("Group A", "Group B"), lwd = 3, cex = 1.1,
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

