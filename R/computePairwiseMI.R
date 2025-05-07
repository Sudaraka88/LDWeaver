#' perform_MI_computation
#'
#' Function to compute pairwise MI (long range links saved to text file), model
#' the background LD levels within each cluster (plots saved) and
#' compute short-range p-values for each short-range links. If requested, this
#' function can also run the ARACNE algorithm
#'
#' @import MatrixExtra
#' @importFrom Matrix rowSums
#' @importFrom fitdistrplus fitdist
#' @importFrom data.table as.data.table
#' @importFrom ggplot2 ggplot aes geom_point ggtitle ggsave
#' @importFrom dplyr `%>%` summarise
#' @importFrom stats quantile coef fitted pbeta
#' @importFrom utils write.table
#' @importFrom RcppArmadillo fastLm
#' @importFrom utils txtProgressBar
#'
#' @param snp.dat output from parsing the multi fasta alignment using LDWeaver::parse_fasta_alignment()
#' @param hdw vector of Hamming distance weights, output from LDWeaver::estimate_Hamming_distance_weights()
#' @param cds_var output from LDWeaver::estimate_variation_in_CDS()
#' @param ncores specify the number of cores to use
#' @param lr_save_path specify the location to save long range MI links as a tsv file (default = NULL, will be auto set)
#' @param sr_save_path specify the location to save short range MI links as a tsv file (default = NULL, will be auto set), links below srp_cutoff will not be saved
#' @param plt_folder specify the folder to save generated plots (default = NULL, will be saved to a folder called PLOTS in getwd())
#' @param sr_dist specify the short-range basepair separation (default = 20000)
#' @param lr_retain_links specify the maximum number of long-range MI links to retain (default = 1000000) - in each block, only a top subset of links will be saved
#' @param max_blk_sz specify maximum block size for MI computation (default = 10000), larger sizes require more RAM
#' @param srp_cutoff specify the short-range -log10(p) cut-off value to discard short-range links before returning the data.frame. This setting has no impact on the
#' modelling since all links are used. However, setting a threshold > 2 will generally reduce the memory usage, plotting time (default = 3, i.e. corresponding to p = 0.001),
#' and run time for ARACNE. If all links are required to be returned, set to 0 (i.e. corresponding to p = 1)
#' @param runARACNE specify whether to run ARACNE on short-range links (default = TRUE), if set to FAULT, all links will be marked as ARACNE=1
#' @param perform_SR_analysis_only skip the long range link analysis (default = FALSE)
#' @param order_links return and save links after sorting in short-range p-value order, most to least important (default = T)
#' @param mega_dset set TRUE for mega scale datasets (default = F)
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
                                  sr_dist = 20000, lr_retain_links = 1e6, max_blk_sz = 10000, srp_cutoff = 3, runARACNE = TRUE,
                                  perform_SR_analysis_only = FALSE, order_links = T, mega_dset = F){

  ## DEBUG LINES - DO NOT DELETE and REMEMBER TO COMMENT
  # lr_save_path = "testscript_op/lr_links_spam.tsv"
  # sr_save_path = "testscript_op/sr_links_spam.tsv"
  # plt_folder = "testscript_op"
  # sr_dist = 20000; lr_retain_links = 1e6; max_blk_sz = 10000; srp_cutoff = 3
  # Rcpp::sourceCpp("src/computeMI.cpp"); Rcpp::sourceCpp("src/fintersect.cpp")

  ## 20240507: MatrixExtra 0.1.15 only accepts a limited number of signatures for parallel crossprod/tcrossprod see ?MatrixExtra::matmult - Changed all signatures to maintain multicore functionality

  t000 = Sys.time()
  # TODO: if no paths are given, we need a way to stop overwriting (use timestamp()?)
  if(is.null(lr_save_path)) lr_save_path = file.path(getwd(), "lr_links.tsv")
  if(is.null(sr_save_path)) sr_save_path = file.path(getwd(), "sr_links.tsv")
  if(is.null(plt_folder)) plt_folder = file.path(getwd(), "PLOTS")

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
  # hsq = as(diag(sqrt(hdw)), 'RsparseMatrix') # Try making this an RsparseMatrix
  ### mega datasets
  if(mega_dset){ # Using SPAM
    if(!requireNamespace("spam") & !requireNamespace("spam64")){
      message("This feature requires spam and spam64 packages.")
      return(invisible())
    } else {
      # hsq = spam::as.spam(hsq)
    }
    hsq = as(diag(sqrt(hdw)), 'RsparseMatrix')# This is likely not optimal for spam, maintain previous behaviour
  } else {
    hsq = as(diag(sqrt(hdw)), 'RsparseMatrix') # Try making this an RsparseMatrix
  }

  if(!perform_SR_analysis_only){
    # estimate the number of lr_links to decide on a discard threshold using  all or 1% of SNPs
    snp_subset = min(snp.dat$nsnp, round(snp.dat$nsnp*0.1))
    set.seed(1988) # set a seed, so that the same subset of positions get picked below (not absolutely needed, but good for reproducibility)
    lr_link_count = sapply(snp.dat$POS[sample(snp.dat$nsnp, snp_subset)], function(x) sum((0.5*snp.dat$g - abs((x - snp.dat$POS)%%snp.dat$g  - 0.5*snp.dat$g))>sr_dist))
    lr_links_approx = sum(lr_link_count)/snp_subset*snp.dat$nsnp/2

  } else {
    lr_links_approx = NULL
  }

  for(i in 1:nblcks){
    t0 = Sys.time()
    cat(paste("Block", i, "of", nblcks, "..."))
    from = MI_cmp_blks$from_s[i]:MI_cmp_blks$from_e[i]
    to = MI_cmp_blks$to_s[i]:MI_cmp_blks$to_e[i]
    # cat("**debug** snp_subset =",snp_subset, " lr_links_approx=",lr_links_approx," lr_retain_links=",lr_retain_links, " **debug**\n")
    sr_links = perform_MI_computation_ACGTN(snp.dat = snp.dat, neff = neff, hsq = hsq, cds_var = cds_var,
                                            lr_save_path = lr_save_path, from = from, sr_dist = sr_dist,
                                            lr_retain_links = lr_retain_links, to = to, ncores = ncores, sr_links = sr_links,
                                            perform_SR_analysis_only = perform_SR_analysis_only, lr_links_approx = lr_links_approx,
                                            mega_dset = mega_dset)

    cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))
  }
  # In case a link gets through with len > sr_dist, we should filter that out in <mergeNsort_sr_links()>
  sr_links_all = mergeNsort_sr_links(cds_var = cds_var, sr_links = sr_links, sr_dist = sr_dist, plt_path = plt_folder, srp_cutoff = srp_cutoff)
  sr_links_red = sr_links_all$sr_links_red
  sr_links_ARACNE_check = sr_links_all$sr_links_ARACNE_check # do not pass outside this function
  rm(sr_links_all)
  # sr_links_red = sr_links_df[sr_links_df$srp_max > srp_cutoff, ] # This is an arbitrary filter for a nice plot and quicker processing

  if(runARACNE){
    # cat(paste("Running ARACNE on", nrow(sr_links_red), "links... \n"))
    ARACNE = runARACNE(sr_links_red, sr_links_ARACNE_check)
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

perform_MI_computation_ACGTN = function(snp.dat, from, to, neff, hsq, cds_var,
                                        lr_save_path, ncores, sr_links, sr_dist = 20000,
                                        lr_retain_links = lr_retain_links, perform_SR_analysis_only = F,
                                        lr_links_approx = NULL, mega_dset = F){
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


  ### mega datasets
  if(mega_dset){ # Using SPAM

    tAf = spam::as.matrix(snp.dat$snp.matrix_A[from, ]); tAfh = tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2)
    tCf = spam::as.matrix(snp.dat$snp.matrix_C[from, ]); tCfh = tcrossprod(tCf, hsq); pCf = Matrix::rowSums(tCfh^2)
    tGf = spam::as.matrix(snp.dat$snp.matrix_G[from, ]); tGfh = tcrossprod(tGf, hsq); pGf = Matrix::rowSums(tGfh^2)
    tTf = spam::as.matrix(snp.dat$snp.matrix_T[from, ]); tTfh = tcrossprod(tTf, hsq); pTf = Matrix::rowSums(tTfh^2)
    tNf = spam::as.matrix(snp.dat$snp.matrix_N[from, ]); tNfh = tcrossprod(tNf, hsq); pNf = Matrix::rowSums(tNfh^2)

    if(fromISto){
      tAt = tAf; tAth = tAfh; pAt = pAf
      tCt = tCf; tCth = tCfh; pCt = pCf
      tGt = tGf; tGth = tGfh; pGt = pGf
      tTt = tTf; tTth = tTfh; pTt = pTf
      tNt = tNf; tNth = tNfh; pNt = pNf
    } else {
      # to
      tAt = spam::as.matrix(snp.dat$snp.matrix_A[to, ]); tAth = tcrossprod(tAt, hsq); pAt = Matrix::rowSums(tAth^2)
      tCt = spam::as.matrix(snp.dat$snp.matrix_C[to, ]); tCth = tcrossprod(tCt, hsq); pCt = Matrix::rowSums(tCth^2)
      tGt = spam::as.matrix(snp.dat$snp.matrix_G[to, ]); tGth = tcrossprod(tGt, hsq); pGt = Matrix::rowSums(tGth^2)
      tTt = spam::as.matrix(snp.dat$snp.matrix_T[to, ]); tTth = tcrossprod(tTt, hsq); pTt = Matrix::rowSums(tTth^2)
      tNt = spam::as.matrix(snp.dat$snp.matrix_N[to, ]); tNth = tcrossprod(tNt, hsq); pNt = Matrix::rowSums(tNth^2)
      # }

    }

  } else {
    ### MATRIX MODE CODE DO NOT CHANGE BELOW
    # from
    # tAf = as(snp.dat$snp.matrix_A[from, ], 'unpackedMatrix'); tAfh = MatrixExtra::tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2) # crashes in linuxMint
    tAf = as(snp.dat$snp.matrix_A[from, ], 'matrix'); tAfh = tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2)
    tCf = as(snp.dat$snp.matrix_C[from, ], 'matrix'); tCfh = tcrossprod(tCf, hsq); pCf = Matrix::rowSums(tCfh^2)
    tGf = as(snp.dat$snp.matrix_G[from, ], 'matrix'); tGfh = tcrossprod(tGf, hsq); pGf = Matrix::rowSums(tGfh^2)
    tTf = as(snp.dat$snp.matrix_T[from, ], 'matrix'); tTfh = tcrossprod(tTf, hsq); pTf = Matrix::rowSums(tTfh^2)
    tNf = as(snp.dat$snp.matrix_N[from, ], 'matrix'); tNfh = tcrossprod(tNf, hsq); pNf = Matrix::rowSums(tNfh^2)

    if(fromISto){
      tAt = tAf; tAth = tAfh; pAt = pAf
      tCt = tCf; tCth = tCfh; pCt = pCf
      tGt = tGf; tGth = tGfh; pGt = pGf
      tTt = tTf; tTth = tTfh; pTt = pTf
      tNt = tNf; tNth = tNfh; pNt = pNf
    } else {
      # to
      tAt = as(snp.dat$snp.matrix_A[to, ], 'matrix'); tAth = tcrossprod(tAt, hsq); pAt = Matrix::rowSums(tAth^2)#+rt*0.5
      tCt = as(snp.dat$snp.matrix_C[to, ], 'matrix'); tCth = tcrossprod(tCt, hsq); pCt = Matrix::rowSums(tCth^2)#+rt*0.5
      tGt = as(snp.dat$snp.matrix_G[to, ], 'matrix'); tGth = tcrossprod(tGt, hsq); pGt = Matrix::rowSums(tGth^2)#+rt*0.5
      tTt = as(snp.dat$snp.matrix_T[to, ], 'matrix'); tTth = tcrossprod(tTt, hsq); pTt = Matrix::rowSums(tTth^2)#+rt*0.5
      tNt = as(snp.dat$snp.matrix_N[to, ], 'matrix'); tNth = tcrossprod(tNt, hsq); pNt = Matrix::rowSums(tNth^2)#+rt*0.5
    }
  }

  den = neff + tcrossprod(snp.dat$r[from], snp.dat$r[to]) * 0.5
  rft = t(tcrossprod(rf, rt))*0.25
  rf = 0.5*rf
  rt = 0.5*rt
  # snp_i = 3
  # snp_j = 41
  # {
  # t0 = Sys.time()
  MI = matrix(rep(0, length(from)*length(to)), nrow = length(from))
  # tX = tAfh; tY = tAth; pX = pAf; pY = pAt; rX = rf; rY = rt; RXY = rft; uqX = uqf[,1]; uqY = uqt[,1]
  computeMI_Sprase(MI, tAfh, tAth, pAf, pAt, rf, rt, rft, uqf[,1], uqt[,1], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tAfh, tCth, pAf, pCt, rf, rt, rft, uqf[,1], uqt[,2], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tAfh, tGth, pAf, pGt, rf, rt, rft, uqf[,1], uqt[,3], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tAfh, tTth, pAf, pTt, rf, rt, rft, uqf[,1], uqt[,4], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tAfh, tNth, pAf, pNt, rf, rt, rft, uqf[,1], uqt[,5], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])

  computeMI_Sprase(MI, tCfh, tAth, pCf, pAt, rf, rt, rft, uqf[,2], uqt[,1], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tCfh, tCth, pCf, pCt, rf, rt, rft, uqf[,2], uqt[,2], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tCfh, tGth, pCf, pGt, rf, rt, rft, uqf[,2], uqt[,3], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tCfh, tTth, pCf, pTt, rf, rt, rft, uqf[,2], uqt[,4], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tCfh, tNth, pCf, pNt, rf, rt, rft, uqf[,2], uqt[,5], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])

  computeMI_Sprase(MI, tGfh, tAth, pGf, pAt, rf, rt, rft, uqf[,3], uqt[,1], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tGfh, tCth, pGf, pCt, rf, rt, rft, uqf[,3], uqt[,2], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tGfh, tGth, pGf, pGt, rf, rt, rft, uqf[,3], uqt[,3], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tGfh, tTth, pGf, pTt, rf, rt, rft, uqf[,3], uqt[,4], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tGfh, tNth, pGf, pNt, rf, rt, rft, uqf[,3], uqt[,5], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])

  computeMI_Sprase(MI, tTfh, tAth, pTf, pAt, rf, rt, rft, uqf[,4], uqt[,1], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tTfh, tCth, pTf, pCt, rf, rt, rft, uqf[,4], uqt[,2], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tTfh, tGth, pTf, pGt, rf, rt, rft, uqf[,4], uqt[,3], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tTfh, tTth, pTf, pTt, rf, rt, rft, uqf[,4], uqt[,4], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tTfh, tNth, pTf, pNt, rf, rt, rft, uqf[,4], uqt[,5], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])

  computeMI_Sprase(MI, tNfh, tAth, pNf, pAt, rf, rt, rft, uqf[,5], uqt[,1], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tNfh, tCth, pNf, pCt, rf, rt, rft, uqf[,5], uqt[,2], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tNfh, tGth, pNf, pGt, rf, rt, rft, uqf[,5], uqt[,3], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tNfh, tTth, pNf, pTt, rf, rt, rft, uqf[,5], uqt[,4], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])
  computeMI_Sprase(MI, tNfh, tNth, pNf, pNt, rf, rt, rft, uqf[,5], uqt[,5], den, ncores) #, mega_dset = mega_dset); #print(MI[snp_i,snp_j])

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
    # disc_thresh = stats::quantile(MI_df_lr$MI, lr_retain_level)
    n_lr_links = nrow(MI_df_lr)
    # Following prob works well for most cases where n_lr_links_total >> lr_retain_links (i.e. save 1M out of 1B),
    # set the maximum to 0 (i.e. no negative values), should we set the max to 1 if prob > 1
    prob = max(c(0, (1 - ((lr_retain_links * (n_lr_links / lr_links_approx)) / n_lr_links))))
    # disc_thresh = Rfast2::Quantile(MI_df_lr$MI, probs = prob)
    disc_thresh = stats::quantile(MI_df_lr$MI, probs = prob)

    # cat("**debug** n_lr_links =",n_lr_links," prob=",prob, " disc_thresh=",disc_thresh, " **debug**\n")

    len_filt = MI_df_lr$MI >= disc_thresh
    cat(paste("... Adding ", sum(len_filt) ," LR links with MI>",round(disc_thresh, 3) ," to file ...", sep = "")) # debug
    if(!all(len_filt == FALSE)){
      MI_df_lr = MI_df_lr[len_filt, ]
      write.table(x = MI_df_lr, file = lr_save_path, append = T, quote = F, row.names = F, col.names = F, sep = '\t')
    }
  }

  if(sr_present){
    # write.table(x = MI_df_sr, file = sr_save_path, append = T, quote = F, row.names = F, col.names = F)
    n_sr_links = nrow(MI_df_sr)
    clust_mat = matrix(c(MI_df_sr$clust1, MI_df_sr$clust2), nrow = n_sr_links)

    # quickly loop through across nclust and append link info to the relevant object
    for(i in 1:cds_var$nclust){
      clust_link_idx = which(.compareToRow(clust_mat, i))
      if(length(clust_link_idx) > 0){
        sr_links[[i]] <- rbind(sr_links[[i]], MI_df_sr[clust_link_idx, ])
      }

    }
    cat(paste("... Adding", n_sr_links, "SR links to list ...")) # debug

  }
  return(sr_links)
  # cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))


}

# tAfh, tAth, pAf, pAt, rf, rt, rft, uqf[,1], uqt[,1]
# computeMI_Sprase = function(MI_t, tX, tY, pX, pY, rX, rY, RXY, uqX, uqY, den, ncores, mega_dset = F){ # for deletion
computeMI_Sprase = function(MI_t, tX, tY, pX, pY, rX, rY, RXY, uqX, uqY, den, ncores){
  pxy_t = tcrossprod(tX, as(tY, 'RsparseMatrix')) + 0.5; pxy_t = as(pxy_t, "matrix") # convert to dense mx
  uq_t = tcrossprod(uqX, uqY)
  pXrX = as(tcrossprod(as(pX*rX, 'RsparseMatrix'), rep(1,length(pY))), 'matrix')
  pYrY = as(tcrossprod(rep(1,length(pX)), as(pY*rY, 'RsparseMatrix')), 'matrix') # same as above, already transposed
  pxpy_tt = tcrossprod(pX, pY) # + RXY + pXrX + pYrY
  .fastHadamard(MI_t, den, uq_t, pxy_t, pxpy_tt, RXY, pXrX, pYrY, ncores)

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
    maxvls = sr_links_t %>% dplyr::group_by(len) %>% dplyr::summarise(max = stats::quantile(MI, 0.95))

    cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))

    cat("Step 2/3 - Model decay ...")
    t0 = Sys.time()
    fit_l = RcppArmadillo::fastLm(X = cbind(log(maxvls$len), 1), y = log(maxvls$max))
    maxvls$fit = mean_dist = exp(fitted(fit_l))
    p_cf = ggplot2::ggplot(data = maxvls, ggplot2::aes(x = len)) +
      ggplot2::geom_point(ggplot2::aes(y = max)) +
      ggplot2::geom_line(ggplot2::aes(y = fit), col = "red") +
      ggplot2::ggtitle(paste("Clust", i)) +
      ggplot2::xlab("Basepair separation") +
      ggplot2::ylab("MI (95th percentile)")

    # print(summary(maxvls$max))
    # plot(maxvls)
    saveRDS(object = maxvls, file = file.path(plt_path, paste("c", i, "_fit_data.rds", sep= "")))
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

    duplink_dt_red = duplink_dt[, .I[which.max(srp_max)], by = keys]$V1

    sr_links_df = rbind(sr_links_df,
                        duplink_df[duplink_dt_red, ])
  }

  # if we do the sorting after ARACNE, it could be faster!
  sr_links_red = sr_links_df[sr_links_df$srp_max > srp_cutoff, ] # This is an arbitrary filter for a nice plot and quicker processing
  sr_links_ARACNE_check = sr_links_df[sr_links_df$MI >= min(sr_links_red$MI), ]

  cat(paste("All done in", round(difftime(Sys.time(), t00, units = "secs"), 2), "s \n"))
  return(list(sr_links_red = sr_links_red, sr_links_ARACNE_check = sr_links_ARACNE_check))

}

