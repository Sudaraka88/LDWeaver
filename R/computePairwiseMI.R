#' perform_MI_computation
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
#' @param lr_retain_level specify the long-range MI retaining percentile (default = 0.99) - in each block, only the top 1\% of lr MI links will be retained
#' @param max_blk_sz specify maximum block size for MI computation (default = 10000), larger sizes require more RAM
#' @param srp_cutoff specify the short-range -log10(p) cut-off value to discard short-range links before returning the data.frame. This setting has no impact on the
#' modelling since all links are used. However, setting a threshold > 2 will generally reduce the memory usage, plotting time (default = 3, i.e. corresponding to p = 0.001),
#' and run time for ARACNE. If all links are required to be returned, set to 0 (i.e. corresponding to p = 1)
#' @param runARACNE specify whether to run ARACNE (default = TRUE), if set to FAULT, all links will be marked as ARACNE=1
#' @param perform_SR_analysis_only skip the long range link analysis (default = FALSE)
#' @param order_links return and save links after sorting in short-range p-value order, most to least important (default = T)
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
perform_MI_computation = function(snp.dat, hdw, cds_var, ncores, lr_save_path = NULL, sr_save_path = NULL, plt_folder = NULL,
                                  sr_dist = 20000, lr_retain_level = 0.99, max_blk_sz = 10000, srp_cutoff = 3, runARACNE = TRUE,
                                  perform_SR_analysis_only = FALSE, order_links = T){
  t000 = Sys.time()
  # TODO: if no paths are given, we need a way to stop overwriting (use timestamp()?)
  if(is.null(lr_save_path)) lr_save_path = file.path(getwd(), "lr_links.tsv")
  if(is.null(sr_save_path)) sr_save_path = file.path(getwd(), "sr_links.tsv")
  if(is.null(plt_folder)) plt_folder = file.path(getwd(), "PLOTS");

  if(!file.exists(plt_folder)) dir.create(plt_folder)

  cat("Begin MI computation... \n")
  # Break down the computation into blocks
  max_blk_sz = round(max_blk_sz, -3) # rounding to 1000s
  MI_cmp_blks = make_blocks(snp.dat$nsnp, max_blk_sz)
  nblcks = nrow(MI_cmp_blks)

  sr_links = list() # list to hold short-range links (fast computation and avoid indexing)
  for(i in 1:cds_var$nclust) sr_links[[i]] = data.frame()

  # pre-computed parameters
  neff = sum(hdw)
  hsq = diag(sqrt(hdw))


  for(i in 1:nblcks){
    t0 = Sys.time()
    cat(paste("Block", i, "of", nblcks, "..."))
    from_ = MI_cmp_blks$from_s[i]:MI_cmp_blks$from_e[i]
    to_ = MI_cmp_blks$to_s[i]:MI_cmp_blks$to_e[i]
    # sr_links = perform_MI_computation_ACGTN(snp.dat = snp.dat, hdw = hdw, from = MI_cmp_blks$from_s[i]:MI_cmp_blks$from_e[i],
    #                                         to = MI_cmp_blks$to_s[i]:MI_cmp_blks$to_e[i], paint = cds_var$paint,
    #                                         nclust = cds_var$nclust, sr_dist = sr_dist, lr_retain_level = lr_retain_level,
    #                                         lr_save_path = lr_save_path, ncores = ncores, sr_links = sr_links)
    sr_links = perform_MI_computation_ACGTN(snp.dat = snp.dat, neff = neff, hsq = hsq, cds_var = cds_var,
                                            lr_save_path = lr_save_path, from = from_, sr_dist = 20000,
                                            lr_retain_level = 0.99, to = to_, ncores = ncores, sr_links = sr_links,
                                            perform_SR_analysis_only = perform_SR_analysis_only)

    cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))
  }
  # In case a link gets through with len > sr_dist, we should filter that out in <mergeNsort_sr_links()>
  sr_links_red = mergeNsort_sr_links(cds_var = cds_var, sr_links = sr_links, sr_dist = sr_dist, plt_path = plt_folder, srp_cutoff = srp_cutoff)
  # sr_links_red = sr_links_df[sr_links_df$srp_max > srp_cutoff, ] # This is an arbitrary filter for a nice plot and quicker processing

  if(runARACNE){
    cat(paste("Running ARACNE on", nrow(sr_links_red), "links... \n"))
    ARACNE = runAracne(sr_links_red)
    sr_links_red$ARACNE = as.numeric(ARACNE)
  } else {
    warning('ARACNE not run, all values will be set to 1')
    sr_links_red$ARACNE = 1
  }

  # sorting after ARACNE
  if(order_links){ # not ordering can speed up annotations
    sr_links_red = sr_links_red[order(sr_links_red$srp_max, decreasing = T), ]
    rownames(sr_links_red) = NULL
  }

  # can we omit this save? sr_links_annotated is the annotated version of this
  write.table(x = sr_links_red, file = sr_save_path, append = T, quote = F, row.names = F, col.names = F, sep = '\t')
  cat(paste("All done in", round(difftime(Sys.time(), t000, units = "mins"), 2), "mins \n"))

  return(sr_links_red) # will be returned as a data.frame for easy continuation of pipeline

}

make_blocks = function(nsnp, max_blk_sz = 10000){ # create the blocks (from_s, from_e, to_s, to_e) for <perform_MI_computation_ACGTN>
  # TODO: we should probably choose the block size based on RAM availability?
  fromtodf = data.frame()
  part1 = ceiling(nsnp/max_blk_sz)
  from_s = c()
  from_e = c()
  for(i in 1:part1){
    from_s = c(from_s, (i-1)*max_blk_sz + 1)
    from_e = c(from_e, min(i*max_blk_sz, nsnp))
  }

  for(i in 1:part1){
    for(j in i:part1){
      fromtodf = rbind(fromtodf, c(from_s[i], from_e[i], from_s[j], from_e[j]))
    }
  }
  colnames(fromtodf) = c("from_s", "from_e", "to_s", "to_e")
  return(fromtodf)
}

perform_MI_computation_ACGTN = function(snp.dat, from, to, neff, hsq, cds_var, lr_save_path, ncores, sr_links, sr_dist = 20000, lr_retain_level = 0.99, perform_SR_analysis_only = F){
  # from, to are vectors

  # These are static, best passed in here <potential inputs>
  # neff = sum(hdw)
  # hsq = diag(sqrt(hdw))
  POS_f = as.numeric(snp.dat$POS[from])
  POS_t = as.numeric(snp.dat$POS[to])

  if(perform_SR_analysis_only){ # We can save (a lot!) of time if only SR analysis is required

    # perform a length check on <POS_f, POS_t> with all <POS_t, POS_f> - drop sites that don't form links < sr_dist
    kp_f = sapply(POS_f, function(x) any( abs(0.5*snp.dat$g - abs((POS_t - x)%%snp.dat$g  - 0.5*snp.dat$g)) < sr_dist) )
    kp_t = sapply(POS_t, function(x) any( abs(0.5*snp.dat$g - abs((POS_f - x)%%snp.dat$g  - 0.5*snp.dat$g)) < sr_dist) )
    from = from[kp_f]
    to = to[kp_t]
    POS_f = as.numeric(snp.dat$POS[from])
    POS_t = as.numeric(snp.dat$POS[to])

  }



  # paint_f = paint[from]
  # paint_t = paint[to]
  paint_f = cds_var$paint[from]
  paint_t = cds_var$paint[to]


  if(length(from) == length(to)){
    fromISto = all(from == to) # center square mx, no need to repeat certain calculations
  } else {
    fromISto = FALSE
  }

  rf = snp.dat$r[from]; if(fromISto) rt = rf else rt = snp.dat$r[to]
  uqf = snp.dat$uqe[from, ]; if(fromISto) uqt = uqf else uqt = snp.dat$uqe[to, ]

  # from
  {
    # tAf = as(snp.dat$snp.matrix_A[from, ], 'unpackedMatrix'); tAfh = MatrixExtra::tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2) # crashes in linuxMint
    tAf = as(snp.dat$snp.matrix_A[from, ], 'lgeMatrix'); tAfh = MatrixExtra::tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2)
    tCf = as(snp.dat$snp.matrix_C[from, ], 'lgeMatrix'); tCfh = MatrixExtra::tcrossprod(tCf, hsq); pCf = Matrix::rowSums(tCfh^2)
    tGf = as(snp.dat$snp.matrix_G[from, ], 'lgeMatrix'); tGfh = MatrixExtra::tcrossprod(tGf, hsq); pGf = Matrix::rowSums(tGfh^2)
    tTf = as(snp.dat$snp.matrix_T[from, ], 'lgeMatrix'); tTfh = MatrixExtra::tcrossprod(tTf, hsq); pTf = Matrix::rowSums(tTfh^2)
    tNf = as(snp.dat$snp.matrix_N[from, ], 'lgeMatrix'); tNfh = MatrixExtra::tcrossprod(tNf, hsq); pNf = Matrix::rowSums(tNfh^2)
  }

  if(fromISto){
    tAt = tAf; tAth = tAfh; pAt = pAf
    tCt = tCf; tCth = tCfh; pCt = pCf
    tGt = tGf; tGth = tGfh; pGt = pGf
    tTt = tTf; tTth = tTfh; pTt = pTf
    tNt = tNf; tNth = tNfh; pNt = pNf
  } else {
    # to
    tAt = as(snp.dat$snp.matrix_A[to, ], 'lgeMatrix'); tAth = MatrixExtra::tcrossprod(tAt, hsq); pAt = Matrix::rowSums(tAth^2)#+rt*0.5
    tCt = as(snp.dat$snp.matrix_C[to, ], 'lgeMatrix'); tCth = MatrixExtra::tcrossprod(tCt, hsq); pCt = Matrix::rowSums(tCth^2)#+rt*0.5
    tGt = as(snp.dat$snp.matrix_G[to, ], 'lgeMatrix'); tGth = MatrixExtra::tcrossprod(tGt, hsq); pGt = Matrix::rowSums(tGth^2)#+rt*0.5
    tTt = as(snp.dat$snp.matrix_T[to, ], 'lgeMatrix'); tTth = MatrixExtra::tcrossprod(tTt, hsq); pTt = Matrix::rowSums(tTth^2)#+rt*0.5
    tNt = as(snp.dat$snp.matrix_N[to, ], 'lgeMatrix'); tNth = MatrixExtra::tcrossprod(tNt, hsq); pNt = Matrix::rowSums(tNth^2)#+rt*0.5
  }

  den = neff + MatrixExtra::tcrossprod(snp.dat$r[from], snp.dat$r[to]) * 0.5
  rft = MatrixExtra::t(MatrixExtra::tcrossprod(rf, rt))*0.25
  rf = 0.5*rf
  rt = 0.5*rt
  # snp_i = 3
  # snp_j = 41
  {
    t0 = Sys.time()
    MI = matrix(rep(0, length(from)*length(to)), nrow = length(from))

    computeMI_Sprase(MI, tAfh, tAth, pAf, pAt, rf, rt, rft, uqf[,1], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tAfh, tCth, pAf, pCt, rf, rt, rft, uqf[,1], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tAfh, tGth, pAf, pGt, rf, rt, rft, uqf[,1], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tAfh, tTth, pAf, pTt, rf, rt, rft, uqf[,1], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tAfh, tNth, pAf, pNt, rf, rt, rft, uqf[,1], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])

    computeMI_Sprase(MI, tCfh, tAth, pCf, pAt, rf, rt, rft, uqf[,2], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tCfh, tCth, pCf, pCt, rf, rt, rft, uqf[,2], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tCfh, tGth, pCf, pGt, rf, rt, rft, uqf[,2], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tCfh, tTth, pCf, pTt, rf, rt, rft, uqf[,2], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tCfh, tNth, pCf, pNt, rf, rt, rft, uqf[,2], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])

    computeMI_Sprase(MI, tGfh, tAth, pGf, pAt, rf, rt, rft, uqf[,3], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tGfh, tCth, pGf, pCt, rf, rt, rft, uqf[,3], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tGfh, tGth, pGf, pGt, rf, rt, rft, uqf[,3], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tGfh, tTth, pGf, pTt, rf, rt, rft, uqf[,3], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tGfh, tNth, pGf, pNt, rf, rt, rft, uqf[,3], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])

    computeMI_Sprase(MI, tTfh, tAth, pTf, pAt, rf, rt, rft, uqf[,4], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tTfh, tCth, pTf, pCt, rf, rt, rft, uqf[,4], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tTfh, tGth, pTf, pGt, rf, rt, rft, uqf[,4], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tTfh, tTth, pTf, pTt, rf, rt, rft, uqf[,4], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tTfh, tNth, pTf, pNt, rf, rt, rft, uqf[,4], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])

    computeMI_Sprase(MI, tNfh, tAth, pNf, pAt, rf, rt, rft, uqf[,5], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tNfh, tCth, pNf, pCt, rf, rt, rft, uqf[,5], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tNfh, tGth, pNf, pGt, rf, rt, rft, uqf[,5], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tNfh, tTth, pNf, pTt, rf, rt, rft, uqf[,5], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
    computeMI_Sprase(MI, tNfh, tNth, pNf, pNt, rf, rt, rft, uqf[,5], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])


  }

  # cat("\n")
  # cat(str(MI))
  # cat("\n")
  # Once MI is computed, we need to save it to a text file
  # Diagnoal blocks of the big matrix should only use upper.tri entries
  # Can we also get the clustering done?
  if(fromISto){
    ind <- which( lower.tri(t(MI),diag=FALSE) , arr.ind = TRUE )
  } else {
    ind <- rbind(which( upper.tri(MI, diag=FALSE) , arr.ind = TRUE), which( lower.tri(MI, diag=FALSE) , arr.ind = TRUE))
  }

  # cat("\n")
  # cat(paste("r:", min(ind[,1]), ",R:", max(ind[,1]), ",c:", min(ind[,2]), ",C:", max(ind[,2]), sep = ""))
  # cat("\n")


  lr_present = sr_present = F

  pos2 = POS_f[ind[,1]]
  pos1 = POS_t[ind[,2]]

  clust2 = paint_f[ind[,1]]
  clust1 = paint_t[ind[,2]]


  MI_df = data.frame(pos1 = pos1,
                     pos2 = pos2,
                     clust1 = clust1,
                     clust2 = clust2,
                     len = 0.5*snp.dat$g - abs((pos1 - pos2)%%snp.dat$g  - 0.5*snp.dat$g),
                     MI = MI[ind])

  sr_lr_switch = (MI_df$len <= sr_dist)
  if(all(sr_lr_switch==TRUE)) {
    MI_df_sr = MI_df
    sr_present = T
  } else if(all(sr_lr_switch==FALSE)){
    MI_df_lr = MI_df
    lr_present = T
  } else {
    MI_df_sr = MI_df[sr_lr_switch, ]
    MI_df_lr = MI_df[!sr_lr_switch, ]
    sr_present = lr_present = T
  }

  # Discard MI values w len > sr_dist using discard_MI_threshold
  if(lr_present & !perform_SR_analysis_only){ # if only the SR analysis is requested, discard this section
    # WARNING! setting a low quantile will create a HUGE lr_links.tsv file!
    cat("... Adding LR links to file ...") # debug
    disc_thresh = stats::quantile(MI_df_lr$MI, lr_retain_level)
    len_filt = MI_df_lr$MI >= disc_thresh
    if(!all(len_filt == FALSE)){
      MI_df_lr = MI_df_lr[len_filt, ]
      write.table(x = MI_df_lr, file = lr_save_path, append = T, quote = F, row.names = F, col.names = F, sep = '\t')
    }
  }

  if(sr_present){
    cat("... Updating sr_links list ...") # debug
    # write.table(x = MI_df_sr, file = sr_save_path, append = T, quote = F, row.names = F, col.names = F)
    clust_mat = matrix(c(MI_df_sr$clust1, MI_df_sr$clust2), nrow = nrow(MI_df_sr))

    # quickly loop through across nclust and append link info to the relevant object
    for(i in 1:cds_var$nclust){
      clust_link_idx = which(.compareToRow(clust_mat, i))
      if(length(clust_link_idx) > 0){
        sr_links[[i]] <- rbind(sr_links[[i]], MI_df_sr[clust_link_idx, ])
      }

    }

  }
  return(sr_links)
  # cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))


}

computeMI_Sprase = function(MI_t, tX, tY, pX, pY, rX, rY, RXY, uqX, uqY, den, ncores){
  # t0 = Sys.time();
  pxy_t = MatrixExtra::tcrossprod(tX, tY) + 0.5; pxy_t = as(pxy_t, "matrix") # convert to dense mx
  uq_t = MatrixExtra::tcrossprod(uqX, uqY)
  # rX = 0.5*rX
  # rY = 0.5*rY
  # RXY = t(MatrixExtra::tcrossprod(rX, rY))
  pXrX = MatrixExtra::tcrossprod(pX*rX, rep(1,length(pY)))
  # pYrY = MatrixExtra::t(MatrixExtra::tcrossprod(pY*rY, rep(1,length(pY))))
  pYrY = MatrixExtra::tcrossprod(rep(1,length(pX)), pY*rY) # same as above, already transposed
  pxpy_tt = MatrixExtra::tcrossprod(pX, pY) # + RXY + pXrX + pYrY


  # difftime(Sys.time(), t0)

  # t0 = Sys.time();
  # MI_t = matrix(rep(0, nrow(pxy_t)*ncol(pxy_t)), nrow = nrow(pxy_t))
  # MI_t = uq_t*pxy_t/den*log(pxy_t/pxpy_t*den)
  .fastHadamard(MI_t, den, uq_t, pxy_t, pxpy_tt, RXY, pXrX, pYrY, ncores)
  # difftime(Sys.time(), t0)

  # return(MI_t)
}

mergeNsort_sr_links = function(cds_var, sr_links, sr_dist, plt_path, srp_cutoff){
  sr_links_df = data.frame()
  duplink_df = data.frame()

  if(cds_var$nclust != length(sr_links)){
    warning("Cluster mismatch detected, stopping!")
    return(-1)
  }
  t00 = Sys.time()
  for(i in 1:cds_var$nclust){
    len = MI = fit = srp_max = NULL # avoid dplyr NSE issue (https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/)

    cat(paste("Clust", i,"of", cds_var$nclust, "\nStep 1/3 - Extract links ..."))
    t0 = Sys.time()

    sr_links_t = sr_links[[i]]
    # Add a filter in case there is a long range link here
    sr_links_t = sr_links_t[!is.na(sr_links_t$len), ];
    sr_links_t = sr_links_t[(sr_links_t$len < sr_dist), ]
    sr_links_t = sr_links_t[(sr_links_t$len > 0), ]

    # colnames(sr_links_t) = c("len", "MI")
    maxvls = sr_links_t %>% dplyr::group_by(len) %>% dplyr::summarise(max = quantile(MI, 0.95))

    cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))

    cat("Step 2/3 - Model decay ...")
    t0 = Sys.time()
    fit_l = RcppArmadillo::fastLm(X = cbind(log(maxvls$len), 1), y = log(maxvls$max))
    maxvls$fit = mean_dist = exp(fitted(fit_l))
    p_cf = ggplot2::ggplot(data = maxvls, ggplot2::aes(x = len)) +
      ggplot2::geom_point(ggplot2::aes(y = max)) +
      ggplot2::geom_line(ggplot2::aes(y = fit), col = "red") +
      ggplot2::ggtitle(paste("Clust", i))

    # print(summary(maxvls$max))
    # plot(maxvls)
    ggplot2::ggsave(filename = file.path(plt_path, paste("c", i, "_fit.png", sep= "")), plot =  p_cf, width = 2200, height = 1200, units = "px")

    cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))

    cat("Step 3/3 - Fit links ...")
    t0 = Sys.time()

    # Fitting all links together
    diff_dat = sr_links_t$MI - mean_dist[sr_links_t$len] # val - mean_val (mean_val = fit)
    diff_dat_idx = which(diff_dat > 0) # vals > 0 (i.e. vals > mean_val)
    diff_dat = cbind(diff_dat_idx, diff_dat[diff_dat_idx])

    beta_fit = coef(fitdistrplus::fitdist(diff_dat[,2], "beta")) # beta distributed (Assumption)
    srp_est = -pbeta(diff_dat[,2], shape1 = beta_fit[1], shape2 = beta_fit[2], lower.tail = F, log.p =  T)
    srp_links = cbind(diff_dat[,1], srp_est)
    sr_links_t$srp_max = NA

    sr_links_t$srp_max[srp_links[,1]] = srp_links[,2] # update the links in sr_links_t
    sr_links_t = sr_links_t[!is.na(sr_links_t$srp_max), ]
    # let's have a separate DF for duplicated links
    dup_links = sr_links_t$clust1 != sr_links_t$clust2
    if(any(dup_links))any_dup_links = T else any_dup_links = F

    if(any_dup_links){
      sr_links_df = rbind(sr_links_df, data.frame(clust_c = i, sr_links_t[!dup_links, ]))
      duplink_df = rbind(duplink_df, data.frame(clust_c = i, sr_links_t[dup_links, ]))
    } else {
      sr_links_df = rbind(sr_links_df, data.frame(clust_c = i, sr_links_t))
    }

    cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))

  }
  cat("Cleaning up links ... ")

  if(nrow(duplink_df) > 0){


    duplink_dt = data.table::as.data.table(duplink_df)
    keys = colnames(duplink_dt)[!grepl('srp_max', colnames(duplink_dt))]
    keys = keys[!grepl('clust_c', keys)]

    # duplink_dt_red = duplink_dt[, list(srp_max = max(srp_max), idx = .I), keys]

    duplink_dt_red = duplink_dt[, .I[which.max(srp_max)], by = keys]$V1


    # duplink_dt_red = duplink_dt[, lapply(.SD, max), by = .(srp_max)]
    # duplink_df_red = duplink_df %>% group_by(keys) %>% summarize(max_srp = max(max_srp))

    sr_links_df = rbind(sr_links_df,
                        duplink_df[duplink_dt_red, ])
  }

  # if we do the sorting after ARACNE, it could be faster!
  # sr_links_df = sr_links_df[order(sr_links_df$srp_max, decreasing = T), ]
  # rownames(sr_links_df) = NULL

  # nlinks = nrow(sr_links_df)
  # srp_thresh = -log10(0.01/nlinks) # Bonferroni (until we find a better cut-off)
  # sr_links_red = sr_links_df[!is.na(sr_links_df$srp_max), ]
  # sum(sr_links_df$srp_max > srp_thresh)

  # sr_links_red = sr_links_df[sr_links_df$srp_max > 3, ] # This is an arbitrary filter for a nice plot
  # sr_links_red = sr_links_red[order(sr_links_red$srp_max, decreasing = T), ] # we need to drop as many links as possible before coming here
  # row.names(sr_links_red) = NULL

  sr_links_red = sr_links_df[sr_links_df$srp_max > srp_cutoff, ] # This is an arbitrary filter for a nice plot and quicker processing

  cat(paste("All done in", round(difftime(Sys.time(), t00, units = "secs"), 2), "s \n"))
  return(sr_links_red)

}

runAracne = function(sr_links_red){
  t0 = Sys.time()
  # TODO: make this function faster using openMP
  # links red is the reduced set of links
  # sr_links_df is useless here, no data is taken from it! sr_links_red is a subset of it!
  nlinks = nrow(sr_links_red)
  pos_mat = matrix(c(sr_links_red$pos1, sr_links_red$pos2), nrow = nlinks) # for the reduced link set
  MIs = matrix(sr_links_red$MI)

  # POS = matrix(POS, nrow = length(POS)) # convert to mx for fast searching

  ARACNE = rep(T, nlinks)
  t0 = Sys.time()
  pb = utils::txtProgressBar(min = 1, max = nlinks, initial = 1)

  pX_ = 0
  # pZ_ = 0 # links are not in pos2 order, probably not worth checking pos2 for order

  for(i in 1:nlinks){
    # microbenchmark::microbenchmark(
    #   {
    utils::setTxtProgressBar(pb,i)
    redo_comXZ = F
    pX = pos_mat[i,1]; # X = which(.compareToRow(POS, pX)) #which(POS %in% pX)
    pZ = pos_mat[i,2]; #Z = which(.compareToRow(POS, pZ)) #which(POS %in% pZ)

    if(pX != pX_){ # skip this check if pX is unchanged
      idX = which(.compareToRow(pos_mat, pX)) #decodeIndex(index[[X]])
      matX = cbind(pos_mat[idX,1], pos_mat[idX,2]); matX = matX[matX != pX]
      pX_ = pX
    }
    # if(pZ != pZ_){
    idZ = which(.compareToRow(pos_mat, pZ)) #decodeIndex(index[[Z]])
    matZ = cbind(pos_mat[idZ,1], pos_mat[idZ,2]); matZ = matZ[matZ != pZ]
    # pZ_ = pZ
    # }

    comXZ = Rfast2::Intersect(matX, matZ) # either of matX or matZ should change for a new link, cannot skip!

    if(length(comXZ) > 0){ # This is the only sr_links
      # MI0 = sr_links_df$MI[Rfast2::Intersect(idX, idZ)]
      MI0X = MIs[idX[which(.compareToRow(matrix(matX), comXZ))], 1]
      MI0Z = MIs[idZ[which(.compareToRow(matrix(matZ), comXZ))], 1]
      ARACNE[i] = .compareTriplet(MI0X, MI0Z, MIs[i,1])
    }
  }
  close(pb)
  cat(paste("\nDone in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  return(ARACNE)
}


