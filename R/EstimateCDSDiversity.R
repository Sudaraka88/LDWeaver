#' estimate_variation_in_CDS
#'
#' Function to estimate the variation within each coding region, the output from this function
#' can be used to segment the genome into diversity-based clusters.
#'
#' @importFrom genbankr cds getSeq
#' @importFrom GenomicRanges start width end
#' @importFrom Matrix rowSums colSums
#' @importFrom data.table data.table setattr %between% .I
#' @importFrom stats kmeans
#' @importFrom ggplot2 ggplot geom_point aes ggsave
#'
#' @param gbk output from parsing the genbank file using BacGWES::parse_genbank_file()
#' @param snp.dat output from parsing the multi fasta alignment using BacGWES::parse_fasta_alignment()
#' @param num_clusts_CDS parition to genome into num_clusts_CDS regions using k-means (default = 3)
#' @param ncores specify the number of cores to use
#' @param clust_plt_path specify path to save CDS variation plot
#'
#' @return R list with CDS variation and allele distribution details
#'
#' @examples
#' \dontrun{
#' cds_var <- estimate_variation_in_CDS(gbk, snp.dat, ncores = 10)
#' }
#' @export
estimate_variation_in_CDS = function(gbk, snp.dat, ncores, num_clusts_CDS = 3, clust_plt_path = NULL){
  # This method is only approximate, but much MUCH faster and easier on resources
  # TODO: Include the higher accuracy function
  t0 = Sys.time()
  cds_reg = genbankr::cds(gbk)
  starts = GenomicRanges::start(cds_reg)
  widths = GenomicRanges::width(cds_reg)
  ends = GenomicRanges::end(cds_reg)
  ncds = length(starts)
  var_estimate = rep(NA, ncds)
  # convert ref to a CharacterVector
  ref = unlist(unname(strsplit(as.character(genbankr::getSeq(gbk)), '')))[snp.dat$POS]

  variation = matrix(c(Matrix::rowSums(snp.dat$snp.matrix_A),
                       Matrix::rowSums(snp.dat$snp.matrix_C),
                       Matrix::rowSums(snp.dat$snp.matrix_G),
                       Matrix::rowSums(snp.dat$snp.matrix_T),
                       Matrix::rowSums(snp.dat$snp.matrix_N)),
                     ncol = snp.dat$nsnp, byrow = T)

  # Generate a reference masking mx with 0 at reference allele
  reference = matrix(rep(1, 5*snp.dat$nsnp), nrow = 5); .ACGTN2num(reference, ref, ncores)

  variation_wo_ref = variation*reference # Remove ref. allele sums

  # side-job - prepare ALT allele for VCF format (required for snpEff)
  alpha = c("A", "C", "G", "T", "*")
  alt = apply(variation_wo_ref > 0, 2, function(x) paste(alpha[which(x)], sep = ',', collapse = ','))

  snp_var = Matrix::colSums(variation_wo_ref) # Variation at each SNP

  POS_dt = data.table::data.table(w = snp.dat$POS)
  data.table::setattr(POS_dt, "sorted", "w")
  for(cds in 1:ncds){
    # pos_idx = POS_dt[, .I[POS_dt$w %between% c(starts[cds], ends[cds])]]
    pos_idx = POS_dt[, .I[POS_dt$w %between% c(starts[cds], ends[cds])]]
    if(length(pos_idx) > 0) var_estimate[cds] = sum(snp_var[pos_idx])/widths[cds]
    # Divide by width to normalise between CDS (i.e. long CDS with fewer SNPs have smaller variability)
  }
  cds_idx = !is.na(var_estimate)

  rownames(variation) = c("A", "C", "G", "T", "N")

  var_estimate = var_estimate[cds_idx]
  cds_start = starts[cds_idx]
  cds_end = ends[cds_idx]

  clusts = perform_clustering(var_estimate, nclust = num_clusts_CDS) # perform clustering
  paint = painter(snp.dat$POS, clusts, cds_start, cds_end) # apply paint to intergenic regions

  cds_var = list(var_estimate = var_estimate, cds_start = cds_start, cds_end = cds_end,
                 clusts = clusts, paint = paint, ref = ref, alt = alt, allele_table = variation,
                 nclust = num_clusts_CDS)

  if(is.null(clust_plt_path)) clust_plt_path = "clust_plt.png"

  prepare_cluster_plot(cds_var, clust_plt_path)

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  return(cds_var)
}

# wrapper around k-means to perform clustering
perform_clustering = function(var_estimate, nclust = 3){
  km = stats::kmeans(var_estimate, centers = nclust, nstart = 10)
  # Reshuffle in to descending order
  km_ord = order(table(km$cluster), decreasing = T)
  km_clst_ord = km$cluster
  for(i in 1:length(km_ord)){
    if(i != km_ord[i]){
      km_clst_ord[which(km$cluster == km_ord[i])] = i
    }
  }

  km_ord = order(table(km$cluster), decreasing = T)
  km_clst_ord = km$cluster
  for(i in 1:length(km_ord)){
    if(i != km_ord[i]){
      km_clst_ord[which(km$cluster == km_ord[i])] = i
    }
  }

  cutoff = max(var_estimate[km_clst_ord == 1])
  return(list(km_clst_ord = km_clst_ord, cutoff = cutoff))
}

# Add clustering to regions that missed out (i.e. intergenic regions)
painter = function(POS, clusts, cds_start, cds_end){
  paint = rep(0, length(POS))
  for(i in 1:length(table(clusts$km_clst_ord))){
    c1 = cbind(cds_start[which(clusts$km_clst_ord==i)], cds_end[which(clusts$km_clst_ord==i)])
    for(j in 1:nrow(c1)){
      paint[which(c1[j,1] < POS & POS < c1[j,2])] = i
    }
  }

  ##### let's update paint using actual regions
  begin = 1
  curr_val = paint[1]
  prev_val = paint[1]
  region_mat = c()
  update = FALSE
  for(i in 2:length(paint)){
    if(paint[i] != prev_val){ # new value is different!
      end = i-1 # ended in the previous iteration
      region_mat = cbind(region_mat, c(prev_val, begin, end))
      begin = i
      prev_val = paint[i] # update as the new val
      update = TRUE
    }
    # edge case
    if(i == length(paint)){ # kast entry
      if(update == TRUE) break # there was a change in the last step, shouldn't update again
      region_mat = cbind(region_mat, c(prev_val, begin, i))
    }
    update = FALSE
  }

  # two edge cases:
  # starting SNP(s) not labelled - append label from the right
  if(region_mat[1,1] == 0){
    r = region_mat[,1]
    paint[r[2]:r[3]] = region_mat[1,2] # update paint
    region_mat[1,1] = region_mat[1,2] # update region_mat
  }

  # ending SNP(s) not labelled - append label from the left
  if(region_mat[1,ncol(region_mat)] == 0){
    r = region_mat[,ncol(region_mat)]
    paint[r[2]:r[3]] = region_mat[1,(ncol(region_mat)-1)] # update paint
    region_mat[1,ncol(region_mat)] = region_mat[1,(ncol(region_mat)-1)] # update region_mat
  }

  # Let's take the 0 regions and check on either side
  rm0s = which(region_mat[1,] == 0)
  for(i in 1:length(rm0s)){ # these are the 0 regions
    r = region_mat[,rm0s[i]]
    if(r[2] == r[3]){ # isolated 0, add the value from the left
      paint[r[2]] = region_mat[1,(rm0s[i]-1)]
    } else { # longer region
      ss = round((r[3] - r[2])/2) # mid point of region
      paint[r[2]:(r[2]+ss)] = region_mat[1,(rm0s[i]-1)] # region to the left
      paint[(r[2]+ss+1):r[3]] = region_mat[1,(rm0s[i]+1)] # region to the right
    }
  }
  return(paint)
}

prepare_cluster_plot = function(cds_var, clust_plt_path){
  POS = Variation = Cluster = NULL
  clust_df = data.frame(POS = cds_var$cds_start, Variation = cds_var$var_estimate,
                        Cluster = factor(cds_var$clusts$km_clst_ord))
  p_clust = ggplot2::ggplot(data = clust_df) + ggplot2::geom_point(ggplot2::aes(x = POS, y = Variation, col = Cluster))
  ggplot2::ggsave(plot = p_clust, filename = clust_plt_path, width = 2200, height = 1200, units = "px")

}
