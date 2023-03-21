### For now, let's recompute MI and extract only the mean summary, should eventually combine with the rest of the pipeline...

# lr_links_path = "~/cloudstor/Research/Bacteria/EpistasisLDBreaker/bac_kpneumo/NTUH-K2044_human_kpneumo.filtered_ge001maf_le015gf_gt01-lt05states.L150134n1029.spydrpick_couplings.1-based.14204744edges"
# sr_links_path = "~/Desktop/BacGWES_RUN/kpneumo/sr_links.tsv"
# links_from_spydrpick = T
#
# lr_links_path = "~/Desktop/BacGWES_RUN/maela/lr_links.tsv"
# sr_links_path = "~/Desktop/BacGWES_RUN/maela/sr_links.tsv"
# links_from_spydrpick = F
#
# lr_links_path = "~/Desktop/BacGWES_RUN/msch/lr_links.tsv"
# sr_links_path = "~/Desktop/BacGWES_RUN/msch/sr_links.tsv"
# links_from_spydrpick = F
#
# lr_links_path = "~/Desktop/BacGWES_RUN/cjejuni/lr_links.tsv"
# sr_links_path = "~/Desktop/BacGWES_RUN/cjejuni/sr_links.tsv"
# links_from_spydrpick = F
#
#
# genomewide_LDMap(lr_links_path = lr_links_path, sr_links_path = sr_links_path, links_from_spydrpick = links_from_spydrpick)


genomewide_LDMap = function(lr_links_path, sr_links_path, reducer = NULL, links_from_spydrpick = F){
  # quick way without requiring snp_data
  # lr_links = LDWeaver::read_LongRangeLinks("~/Desktop/BacGWES_RUN/msch/lr_links.tsv")
  # sr_links = LDWeaver::read_ShortRangeLinks("~/Desktop/BacGWES_RUN/msch/sr_links.tsv")
  # reducer = 25

  t0 = Sys.time()
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
  }

  # Matrix::isSymmetric(sprs)
  cat(paste("Reducing dimensionality from", length(pos_vec) ,"x", length(pos_vec), "-> "))
  x = mat(length(pos_vec), reducer)
  cat(paste(ncol(x) ,"x", ncol(x), "...\n"))
  reduced_mat = MatrixExtra::crossprod(x, MatrixExtra::crossprod(sprs, x))

  # spam::image(ldmx)
  htm = as(reduced_mat, 'matrix')/reducer^2
  # diag(htm) = NA

  nms = as.character(pos_vec[seq(from = 1, to = length(pos_vec), by = reducer-1)])

  colnames(htm) = rownames(htm) = nms[1:nrow(htm)]

  htm = log10(htm + 1e-5)
  htm = rescale01(htm)

  cat("Preparing plot...")
  heatmap3::heatmap3(htm, symm = T, Rowv = NA, Colv = NA,
                     balanceColor = F, col = colorRampPalette(c("white", "#E1B9B4","#AE452C","#802418"))(2056))


  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 3), "s \n"))
}



rescale01 = function(v){
  mv = min(v)
  rn = max(v) - min(v)

  rsc = (v - mv)/rn
  return(rsc)
}

mat <- function(n, r) {
  suppressWarnings(matrix(c(rep(1, r), rep(0, n)), n, n/r))
}

# showLegend <- function(
#   legend = c("Group A", "Group B"),
#   lwd = 3,
#   cex = 1.1,
#   col = c("red", "blue"),
#   ...
# )
# # heatmap(as(sprs, 'matrix'), symm = T, Rowv = NA, Colv = NA)
# colnames(ldmx) = rownames(ldmx) = as.character(pos_vec)
# heatmap(ldmx, symm = T, Rowv = NA, Colv = NA)
# # heatmap(log10(ldmx), symm = T, Rowv = NA, Colv = NA)

# ComplexHeatmap::Heatmap(ldmx, row_dend_reorder = F
# )

#
#
# spam::image(lr_sprs) # This does not look that good!
#   t000 = Sys.time()
#   cat("Begin MI computation... \n")
#   snp.dat = readRDS("~/Desktop/BacGWES_RUN/spn_gps/snp_ACGTN.rds")
#   hdw = readRDS("~/Desktop/BacGWES_RUN/spn_gps/hdw.rds")
#   cds_var = readRDS("~/Desktop/BacGWES_RUN/spn_gps/cds_var.rds")
#   ncores = parallel::detectCores()
#   max_blk_sz = 100 # think of a clever way to find this
#
#   MI_cmp_blks = make_blocks(snp.dat$nsnp, max_blk_sz)
#   nblcks = nrow(MI_cmp_blks)
#
#   sqrt(2*nblcks)
#
#   # MI_cmp_blks$MMI = 0
#   stats = data.frame()
#
#   # pre-computed parameters
#   neff = sum(hdw)
#   hsq = diag(sqrt(hdw))
#
#   pb = utils::txtProgressBar(min = 1, max = nblcks, initial = 1)
#
#   for(i in 1:nblcks){
#     utils::setTxtProgressBar(pb,i)
#     t0 = Sys.time()
#     # cat(paste("Block", i, "of", nblcks, "..."))
#     from_ = MI_cmp_blks$from_s[i]:MI_cmp_blks$from_e[i]
#     to_ = MI_cmp_blks$to_s[i]:MI_cmp_blks$to_e[i]
#
#     stats = rbind(stats,
#                   find_MI_block_stats(snp.dat = snp.dat, neff = neff, hsq = hsq, cds_var = cds_var,
#                                             from = from_, to = to_, ncores = ncores))
#
#     # cat(paste(" Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))
#   }
#
#
#   MI_cmp_blks$i = as.numeric(as.factor(MI_cmp_blks$from_s))
#   MI_cmp_blks$j = as.numeric(as.factor(MI_cmp_blks$to_s))
#
#   ldmx = Matrix::sparseMatrix(i = c(MI_cmp_blks$i, MI_cmp_blks$j),
#                               j = c(MI_cmp_blks$j, MI_cmp_blks$i),
#                               x = c(stats$mean, stats$mean))
#
#   Matrix::isSymmetric(ldmx)
#
#   # spam::image(ldmx)
#   ldmx = as(ldmx, 'matrix')
#
#
#   diag(ldmx) = NA
#
#   colnames(ldmx) = rownames(ldmx) = as.character(snp.dat$POS[unique(MI_cmp_blks$to_s)])
#   heatmap(ldmx, symm = T, Rowv = NA, Colv = NA)
#   # heatmap(log10(ldmx), symm = T, Rowv = NA, Colv = NA)
#
#   ComplexHeatmap::Heatmap(ldmx, row_dend_reorder = F
#                           )
# }
#
# make_blocks = function(nsnp, max_blk_sz = 10000){ # create the blocks (from_s, from_e, to_s, to_e) for <perform_MI_computation_ACGTN>
#   # TODO: we should probably choose the block size based on RAM availability?
#   fromtodf = data.frame()
#   part1 = ceiling(nsnp/max_blk_sz)
#   from_s = c()
#   from_e = c()
#   for(i in 1:part1){
#     from_s = c(from_s, (i-1)*max_blk_sz + 1)
#     from_e = c(from_e, min(i*max_blk_sz, nsnp))
#   }
#
#   for(i in 1:part1){
#     for(j in i:part1){
#       fromtodf = rbind(fromtodf, c(from_s[i], from_e[i], from_s[j], from_e[j]))
#     }
#   }
#   colnames(fromtodf) = c("from_s", "from_e", "to_s", "to_e")
#   return(fromtodf)
# }
#
#
# # Send in much smaller chunks
# find_MI_block_stats = function(snp.dat, from, to, neff, hsq, cds_var, ncores){
#   # from, to are vectors
#   POS_f = as.numeric(snp.dat$POS[from])
#   POS_t = as.numeric(snp.dat$POS[to])
#
#   paint_f = cds_var$paint[from]
#   paint_t = cds_var$paint[to]
#
#
#   if(length(from) == length(to)){
#     fromISto = all(from == to) # center square mx, no need to repeat certain calculations
#   } else {
#     fromISto = FALSE
#   }
#
#   rf = snp.dat$r[from]; if(fromISto) rt = rf else rt = snp.dat$r[to]
#   uqf = snp.dat$uqe[from, ]; if(fromISto) uqt = uqf else uqt = snp.dat$uqe[to, ]
#
#   # from
#   {
#     # tAf = as(snp.dat$snp.matrix_A[from, ], 'unpackedMatrix'); tAfh = MatrixExtra::tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2) # crashes in linuxMint
#     tAf = as(snp.dat$snp.matrix_A[from, ], 'lgeMatrix'); tAfh = MatrixExtra::tcrossprod(tAf, hsq); pAf = Matrix::rowSums(tAfh^2)
#     tCf = as(snp.dat$snp.matrix_C[from, ], 'lgeMatrix'); tCfh = MatrixExtra::tcrossprod(tCf, hsq); pCf = Matrix::rowSums(tCfh^2)
#     tGf = as(snp.dat$snp.matrix_G[from, ], 'lgeMatrix'); tGfh = MatrixExtra::tcrossprod(tGf, hsq); pGf = Matrix::rowSums(tGfh^2)
#     tTf = as(snp.dat$snp.matrix_T[from, ], 'lgeMatrix'); tTfh = MatrixExtra::tcrossprod(tTf, hsq); pTf = Matrix::rowSums(tTfh^2)
#     tNf = as(snp.dat$snp.matrix_N[from, ], 'lgeMatrix'); tNfh = MatrixExtra::tcrossprod(tNf, hsq); pNf = Matrix::rowSums(tNfh^2)
#   }
#
#   if(fromISto){
#     tAt = tAf; tAth = tAfh; pAt = pAf
#     tCt = tCf; tCth = tCfh; pCt = pCf
#     tGt = tGf; tGth = tGfh; pGt = pGf
#     tTt = tTf; tTth = tTfh; pTt = pTf
#     tNt = tNf; tNth = tNfh; pNt = pNf
#   } else {
#     # to
#     tAt = as(snp.dat$snp.matrix_A[to, ], 'lgeMatrix'); tAth = MatrixExtra::tcrossprod(tAt, hsq); pAt = Matrix::rowSums(tAth^2)#+rt*0.5
#     tCt = as(snp.dat$snp.matrix_C[to, ], 'lgeMatrix'); tCth = MatrixExtra::tcrossprod(tCt, hsq); pCt = Matrix::rowSums(tCth^2)#+rt*0.5
#     tGt = as(snp.dat$snp.matrix_G[to, ], 'lgeMatrix'); tGth = MatrixExtra::tcrossprod(tGt, hsq); pGt = Matrix::rowSums(tGth^2)#+rt*0.5
#     tTt = as(snp.dat$snp.matrix_T[to, ], 'lgeMatrix'); tTth = MatrixExtra::tcrossprod(tTt, hsq); pTt = Matrix::rowSums(tTth^2)#+rt*0.5
#     tNt = as(snp.dat$snp.matrix_N[to, ], 'lgeMatrix'); tNth = MatrixExtra::tcrossprod(tNt, hsq); pNt = Matrix::rowSums(tNth^2)#+rt*0.5
#   }
#
#   den = neff + MatrixExtra::tcrossprod(snp.dat$r[from], snp.dat$r[to]) * 0.5
#   rft = MatrixExtra::t(MatrixExtra::tcrossprod(rf, rt))*0.25
#   rf = 0.5*rf
#   rt = 0.5*rt
#   # snp_i = 3
#   # snp_j = 41
#   {
#     # t0 = Sys.time()
#     MI = matrix(rep(0, length(from)*length(to)), nrow = length(from))
#
#     computeMI_Sprase(MI, tAfh, tAth, pAf, pAt, rf, rt, rft, uqf[,1], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tAfh, tCth, pAf, pCt, rf, rt, rft, uqf[,1], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tAfh, tGth, pAf, pGt, rf, rt, rft, uqf[,1], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tAfh, tTth, pAf, pTt, rf, rt, rft, uqf[,1], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tAfh, tNth, pAf, pNt, rf, rt, rft, uqf[,1], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])
#
#     computeMI_Sprase(MI, tCfh, tAth, pCf, pAt, rf, rt, rft, uqf[,2], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tCfh, tCth, pCf, pCt, rf, rt, rft, uqf[,2], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tCfh, tGth, pCf, pGt, rf, rt, rft, uqf[,2], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tCfh, tTth, pCf, pTt, rf, rt, rft, uqf[,2], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tCfh, tNth, pCf, pNt, rf, rt, rft, uqf[,2], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])
#
#     computeMI_Sprase(MI, tGfh, tAth, pGf, pAt, rf, rt, rft, uqf[,3], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tGfh, tCth, pGf, pCt, rf, rt, rft, uqf[,3], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tGfh, tGth, pGf, pGt, rf, rt, rft, uqf[,3], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tGfh, tTth, pGf, pTt, rf, rt, rft, uqf[,3], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tGfh, tNth, pGf, pNt, rf, rt, rft, uqf[,3], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])
#
#     computeMI_Sprase(MI, tTfh, tAth, pTf, pAt, rf, rt, rft, uqf[,4], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tTfh, tCth, pTf, pCt, rf, rt, rft, uqf[,4], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tTfh, tGth, pTf, pGt, rf, rt, rft, uqf[,4], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tTfh, tTth, pTf, pTt, rf, rt, rft, uqf[,4], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tTfh, tNth, pTf, pNt, rf, rt, rft, uqf[,4], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])
#
#     computeMI_Sprase(MI, tNfh, tAth, pNf, pAt, rf, rt, rft, uqf[,5], uqt[,1], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tNfh, tCth, pNf, pCt, rf, rt, rft, uqf[,5], uqt[,2], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tNfh, tGth, pNf, pGt, rf, rt, rft, uqf[,5], uqt[,3], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tNfh, tTth, pNf, pTt, rf, rt, rft, uqf[,5], uqt[,4], den, ncores); #print(MI[snp_i,snp_j])
#     computeMI_Sprase(MI, tNfh, tNth, pNf, pNt, rf, rt, rft, uqf[,5], uqt[,5], den, ncores); #print(MI[snp_i,snp_j])
#
#
#   }
#
#   return(data.frame(q25 = quantile(MI, 0.25),
#                     q50 = quantile(MI, 0.5),
#                     q75 = quantile(MI, 0.75),
#                     mean = mean(MI),
#                     max = max(MI)))
#
# }
#
#
# computeMI_Sprase = function(MI_t, tX, tY, pX, pY, rX, rY, RXY, uqX, uqY, den, ncores){
#   # t0 = Sys.time();
#   pxy_t = MatrixExtra::tcrossprod(tX, tY) + 0.5; pxy_t = as(pxy_t, "matrix") # convert to dense mx
#   uq_t = MatrixExtra::tcrossprod(uqX, uqY)
#   # rX = 0.5*rX
#   # rY = 0.5*rY
#   # RXY = t(MatrixExtra::tcrossprod(rX, rY))
#   pXrX = MatrixExtra::tcrossprod(pX*rX, rep(1,length(pY)))
#   # pYrY = MatrixExtra::t(MatrixExtra::tcrossprod(pY*rY, rep(1,length(pY))))
#   pYrY = MatrixExtra::tcrossprod(rep(1,length(pX)), pY*rY) # same as above, already transposed
#   pxpy_tt = MatrixExtra::tcrossprod(pX, pY) # + RXY + pXrX + pYrY
#
#
#   # difftime(Sys.time(), t0)
#
#   # t0 = Sys.time();
#   # MI_t = matrix(rep(0, nrow(pxy_t)*ncol(pxy_t)), nrow = nrow(pxy_t))
#   # MI_t = uq_t*pxy_t/den*log(pxy_t/pxpy_t*den)
#   .fastHadamard(MI_t, den, uq_t, pxy_t, pxpy_tt, RXY, pXrX, pYrY, ncores)
#   # difftime(Sys.time(), t0)
#
#   # return(MI_t)
# }
