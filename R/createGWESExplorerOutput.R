#' write_output_for_gwes_explorer
#'
#' Function to generate long and short range GWES plots
#'
#' @importFrom Matrix which
#' @importFrom utils write.table
#'
#'
#' @param snp.dat output from parsing the multi fasta alignment using BacGWES::parse_fasta_alignment()
#' @param srlinks_tophits data frame with top short range GWES links, returned from BacGWES::perform_snpEff_annotations()
#' @param gwes_explorer_folder folder to save the inputs for GWESExplorer
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' write_output_for_gwes_explorer(snp.dat, srlinks_tophits, gwes_explorer_folder)
#' }
#'
#' @export
write_output_for_gwes_explorer = function(snp.dat, srlinks_tophits, gwes_explorer_folder){
  t0 = Sys.time()
  cat("Preparing Outputs for GWESExplorer...")
  if(!file.exists(gwes_explorer_folder)) dir.create(gwes_explorer_folder)
  gwes_explorer_loci_pth = file.path(gwes_explorer_folder, "snps.loci")
  gwes_explorer_aln_pth = file.path(gwes_explorer_folder, "snps.aln")
  gwes_explorer_outliers_pth = file.path(gwes_explorer_folder, "snps.outliers")

  gwex_snps = sort(unique(c(srlinks_tophits$pos1, srlinks_tophits$pos2)))
  gwes_snps_idx = sapply(gwex_snps, function(x) which(snp.dat$POS %in% x))

  # loci file
  if(file.exists(gwes_explorer_loci_pth)) unlink(gwes_explorer_loci_pth)
  write.table(x = gwex_snps, file = gwes_explorer_loci_pth, quote = F, col.names = F, row.names = F)

  # aln file
  ## Using the snp_pos_idx, we need to save a chosen SNP FASTA for gwes_explorer
  fasta = matrix(rep(NA, snp.dat$nseq*length(gwes_snps_idx)), nrow = length(gwes_snps_idx))
  tidx = Matrix::which(snp.dat$snp.matrix_A[gwes_snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "A"
  tidx = Matrix::which(snp.dat$snp.matrix_C[gwes_snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "C"
  tidx = Matrix::which(snp.dat$snp.matrix_G[gwes_snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "G"
  tidx = Matrix::which(snp.dat$snp.matrix_T[gwes_snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "T"
  tidx = Matrix::which(snp.dat$snp.matrix_N[gwes_snps_idx,] == TRUE); if(length(tidx) > 0) fasta[tidx] = "N"
  fasta = t(fasta)

  if(file.exists(gwes_explorer_aln_pth))  unlink(gwes_explorer_aln_pth)
  for(i in 1:snp.dat$nseq){
    write.table(paste(">", snp.dat$seq.names[i], sep = ""), gwes_explorer_aln_pth, quote = F, col.names = F, row.names = F, append = T)
    write.table(paste(fasta[i,], collapse = ""), gwes_explorer_aln_pth, quote = F, col.names = F, row.names = F, append = T)
    # write.table("\n", gwes_explorer_aln_pth, quote = F, col.names = F, row.names = F, append = T)
  }

  # outliers file
  outliers = data.frame(Pos_1 = as.numeric(srlinks_tophits$pos1),
                        Pos_2 = as.numeric(srlinks_tophits$pos2),
                        Distance = as.numeric(srlinks_tophits$len),
                        Direct = as.numeric(srlinks_tophits$ARACNE),
                        MI = as.numeric(srlinks_tophits$srp),
                        MI_wogaps = as.numeric(srlinks_tophits$MI))

  if(file.exists(gwes_explorer_outliers_pth))  unlink(gwes_explorer_outliers_pth)
  write.table(x = outliers, file = gwes_explorer_outliers_pth, quote = F, col.names = T, row.names = F)

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
}
