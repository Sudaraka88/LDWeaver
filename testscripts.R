snp.dat.spam = LDWeaver::parse_fasta_alignment("inst/extdata/sample.aln.gz", mega_dset = T)
snp.dat.mx = LDWeaver::parse_fasta_alignment("inst/extdata/sample.aln.gz")

gbk = LDWeaver::parse_genbank_file("inst/extdata/sample.gbk", length_check = F)

cds_var.spam = LDWeaver::estimate_variation_in_CDS(gbk = gbk$gbk, snp.dat = snp.dat.spam, ncores = 3, mega_dset = T, clust_plt_path = "testscript_op/cds_plt_spam.png" )
cds_var.mx = LDWeaver::estimate_variation_in_CDS(gbk = gbk$gbk, snp.dat = snp.dat.mx, ncores = 3, clust_plt_path = "testscript_op/cds_plt_mx.png")

all.equal(cds_var.mx, cds_var.spam)

hdw.mx = LDWeaver::estimate_Hamming_distance_weights(snp.dat = snp.dat.mx)
hdw.spam = LDWeaver::estimate_Hamming_distance_weights(snp.dat = snp.dat.spam, mega_dset = T)

sr_links.mx = LDWeaver::perform_MI_computation(snp.dat = snp.dat.mx, hdw = hdw.mx, cds_var = cds_var.mx, ncores = 10,
                                            lr_save_path = "testscript_op/lr_links_mx.tsv", sr_save_path = "testscript_op/sr_links_mx.tsv",
                                            plt_folder = "testscript_op", sr_dist = 20000, lr_retain_links = 1e6,
                                            max_blk_sz = 1e4, srp_cutoff = 3, runARACNE = T,
                                            perform_SR_analysis_only = F, order_links = T)

sr_links.spam = LDWeaver::perform_MI_computation(snp.dat = snp.dat.spam, hdw = hdw.spam, cds_var = cds_var.spam, ncores = 10,
                                               lr_save_path = "testscript_op/lr_links_spam.tsv", sr_save_path = "testscript_op/sr_links_spam.tsv",
                                               plt_folder = "testscript_op", sr_dist = 20000, lr_retain_links = 1e6,
                                               max_blk_sz = 1e4, srp_cutoff = 3, runARACNE = T,
                                               perform_SR_analysis_only = F, order_links = T, mega_dset = T)




all.equal(sr_links.mx, sr_links.spam)
sr_links.mx[which(sr_links.mx$pos1 != sr_links.spam$pos1),]

## Are the matrices the same?
# tmp = spam::as.dgCMatrix.spam(snp.dat.spam$snp.matrix_N)
# tmp2 = as(snp.dat$snp.matrix_N, 'dgCMatrix'). ## This is not possible for large matrices
#
# all(tmp@i == tmp2@i)
# all(tmp@p == tmp2@p)
# all(tmp@x == tmp2@x)
# all(tmp@Dim == tmp2@Dim)

## Check the crossprod
# tmpc_mx = Matrix::crossprod(snp.dat$snp.matrix_A)
# tmpc_spam = spam::crossprod.spam(snp.dat.spam$snp.matrix_A)
# tmpc_mx[1:10, 1:10]
# tmpc_spam[1:10, 1:10]
