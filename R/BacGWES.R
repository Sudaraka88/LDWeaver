#' BacGWES
#'
#' Function to compute pairwise MI (long range links saved to text file), model
#' the background LD levels within each cluster (plots saved) and
#' compute short-range p-values for each short-range links. If requested, this
#' function can also run the ARACNE algorithm
#'
#' @importFrom MatrixExtra tcrossprod t
#' @importFrom Matrix rowSums
#' @importFrom fitdistrplus fitdist
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot aes geom_point ggtitle ggsave
#' @importFrom dplyr `%>%` summarise
#' @importFrom stats quantile coef fitted pbeta
#' @importFrom utils write.table
#' @importFrom RcppArmadillo fastLm
#' @importFrom Rfast2 Intersect
#' @importFrom utils txtProgressBar
#'
#' @param snp.dat output from parsing the multi fasta alignment using BacGWES::parse_fasta_alignment()
#' @param hdw vector of Hamming distance weights, output from BacGWES::estimate_Hamming_distance_weights()
#' @param cds_var output from BacGWES::estimate_variation_in_CDS()
#' @param ncores specify the number of cores to use
#' @param lr_save_path specify the location to save long range MI links as a tsv file (default = NULL, will be auto set)
#' @param sr_save_path specify the location to save short range MI links as a tsv file (default = NULL, will be auto set), links below srp_cutoff will not be saved
#' @param plt_folder specify the folder to save generated plots (default = NULL, will be saved to a folder called PLOTS in getwd())
#' @param sr_dist specify the short-range basepair separation (default = 20000)
#' @param discard_MI_threshold_lr specify minimum MI value to retain long range links (default = 0.25)
#' @param max_blk_sz specify maximum block size for MI computation (default = 10000), larger sizes require more RAM
#' @param srp_cutoff specify the short-range -log10(p) cut-off value to discard short-range links before returning the data.frame. This setting has no impact on the
#' modelling since all links are used. However, setting a threshold > 2 will generally reduce the memory usage, plotting time (default = 3, i.e. corresponding to p = 0.001),
#' and run time for ARACNE. If all links are required to be returned, set to 0 (i.e. corresponding to p = 1)
#' @param runARACNE specify whether to run ARACNE (default = TRUE), if set to FAULT, all links will be marked as ARACNE=1
#'
#' @return R data frame with short range GWES links (plots and long range links will be saved)
#'
#' @examples
#' \dontrun{
#' sr_links_red = perform_MI_computation(snp.dat = snp.dat, hdw = hdw, cds_var = cds_var, ncores = 10,
#' lr_save_path = "lr.tsv", plt_folder = "plots")
#' }
#'
#' @export
