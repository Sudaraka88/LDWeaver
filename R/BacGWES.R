#' BacGWES
#'
#' Function to run the BacGWES pipeline
#'
#' @importFrom parallel detectCores
#' @importFrom utils timestamp
#'
#' @param dset name of the dataset, all outputs will be saved to the folder <dset>
#' @param aln_path path to the multi fasta alignment
#' @param gbk_path path to genbank file
#' @param snp_filt_method specify the filtering method for SNP extraction: 'relaxed' or 'default' (default = 'default')
#' @param snpeff_jar_path path to <snpEff.jar>. If unavailable or not required, set SnpEff_Annotate = F
#' @param SnpEff_Annotate specify whether to perform annotations using SnpEff
#' @param sr_dist links less than <sr_dist> apart are considered 'short range' (default = 20000), range 1000 - 25000 bp.
#' @param discard_MI_threshold_lr specify minimum MI value to retain long range links (default = 0.25), range 0-1 (excluding both). Values closer to 0 will retain more values
#' @param max_tophits specify the maximum number of short range links to save as <tophits.tsv>. Note: all short-range links will be annotated (and saved separately),
#' but only the top <max_tophits> will be used for visualisation (default = 250), range 50 - 1000
#' @param num_clusts_CDS parition to genome into num_clusts_CDS regions using k-means (default = 3) range 1 - 10
#' @param srp_cutoff specify the short-range -log10(p) cut-off value to discard short-range links before returning the data.frame. This setting has no impact on the
#' modelling since all links are used. However, setting a threshold > 2 will generally reduce the memory usage, plotting time (default = 3, i.e. corresponding to p = 0.001),
#' and run time for ARACNE. If all links are required to be returned, set to 0 (i.e. corresponding to p = 1), range 0 - 5
#' @param tanglegram_break_segments specify the number of genome segments to prepare - one tanglegram per segment (default = 5), range 1 - 10
#' @param multicore specify whether to use parallel processing (default = T)
#' @param ncores specify the number of cores to use for parallel processing (default = NULL), will auto detect if NULL
#' @param max_blk_sz specify maximum block size for MI computation (default = 10000), larger sizes require more RAM, range 1000 - 100000
#'
#' @return Numeric Value 0 if successful (all generated outputs will be saved)
#'
#' @examples
#' \dontrun{
#' sr_links_red = BacGWES(dset = "efcm", aln_path = "<efcm_aln>", gbk_path = "<efcm.gbk>")
#' }
#'
#' @export
BacGWES = function(dset, aln_path, gbk_path, snp_filt_method = "default", snpeff_jar_path = NULL, SnpEff_Annotate = T,
                   sr_dist = 20000, discard_MI_threshold_lr = 0.25, max_tophits = 250, num_clusts_CDS = 3, srp_cutoff = 3,
                   tanglegram_break_segments = 5, multicore = T, max_blk_sz = 10000, ncores = NULL){
  # Build blocks
  # BLK1: Extract SNPs and create sparse Mx from MSA (fasta)
  # BLK2: Parse GBK
  # BLK3: Estimate diversity within each CDS, cluster and paint < # possible inputs on methods>
  # BLK4: Compute Hamming Distance weights
  # BLK5: Compute MI between all links, sr_links model fitter, ARACNE
  # BLK6: GWES_plots
  # BLK7: Snpeff annotation pipeline, dtermine tophits
  # BLK8: Tanglegram (depends: chromoMap)
  # BLK9: GWESExplorer (depends: GWESExplorer)

  # Welcome message
  timestamp()
  cat(paste("\n\n Performing GWES analysis on:", dset, "\n\n"))

  # Sanity checks
  # annotations
  if(SnpEff_Annotate == T) {
    if(is.null(snpeff_jar_path)) stop("You must specify <snpeff_jar_path> for annotations. To run without annotations, set SnpEff_Annotate = F")
  }

  # multicore
  if(multicore == T) { # multicore requested (default)
    max_n_cores = parallel::detectCores()
    if(!is.null(ncores)){ # ncores provided as well
      if( (ncores <= 1) || (ncores > max_n_cores) ){ # negative ncores or ncores > max_n_cores
        warning(paste("Specified ncores = ", ncores, "unsupported, using", max_n_cores, "instead."))
        ncores = max_n_cores # set to ncores = max_n_cores
      }
    } else { # ncores not provided
      ncores = max_n_cores # set to ncores = max_ncores
    }

  } else ncores = 1  # multicore not requested

  # parameters
  if(sr_dist < 1000 | sr_dist > 25000) {
    warning(paste("Unable to use the provided value for <sr_dist>, using", 20000))
    sr_dist = 20000
  }

  if(discard_MI_threshold_lr <= 0 | discard_MI_threshold_lr >= 1) {
    warning(paste("Unable to use the provided value for <num_clusts_CDS>, using", 0.25))
    discard_MI_threshold_lr = 0.25
  }

  if(max_tophits < 50 | max_tophits > 1000) {
    warning(paste("Unable to use the provided value for <max_tophits>, using", 250))
    sr_dist = 250
  }

  if(num_clusts_CDS < 1 | num_clusts_CDS > 10) {
    warning(paste("Unable to use the provided value for <num_clusts_CDS>, using", 3))
    num_clusts_CDS = 3
  }

  if(srp_cutoff < 0 | srp_cutoff > 5) {
    warning(paste("Unable to use the provided value for <srp_cutoff>, using", 3))
    srp_cutoff = 3
  }

  if(tanglegram_break_segments < 0 | tanglegram_break_segments > 10) {
    warning(paste("Unable to use the provided value for <tanglegram_break_segments>, using", 5))
    tanglegram_break_segments = 5
  }

  if(max_blk_sz < 1000 | max_blk_sz > 100000) {
    warning(paste("Unable to use the provided value for <max_blk_sz>, using", 10000, "...!If this value is causing the function to crash, consider reducing!..."))
    max_blk_sz = 10000
  }




  # setup paths
  if(!file.exists(dset)) dir.create(dset) # save everything in here
  ACGTN_snp_path = file.path(dset, "snp_ACGTN.rds")
  parsed_gbk_path = file.path(dset, "parsed_gbk.rds")
  cds_var_path = file.path(dset, "cds_var.rds")
  hdw_path = file.path(dset, "hdw.rds")
  lr_save_path = file.path(dset, "lr_links.tsv")
  sr_save_path = file.path(dset, "sr_links.tsv")
  tophits_path = file.path(dset, "tophits.tsv")
  clust_plt_path = file.path(dset, "CDS_clustering.png")

  # BLK1
  cat("\n\n #################### BLOCK 1 #################### \n\n")
  if(!file.exists(ACGTN_snp_path)) {
    cat("Performing snp extraction from the alignment \n")
    snp.dat = BacGWES::parse_fasta_alignment(aln_path = aln_path, method = snp_filt_method)
    cat("Savings snp.dat...")
    saveRDS(snp.dat, ACGTN_snp_path)
    cat("Done!\n")
  }else{
    cat("Loading previous snp matrix \n")
    snp.dat = readRDS(ACGTN_snp_path)
  }

  # BLK2
  cat("\n\n #################### BLOCK 2 #################### \n\n")
  if(!file.exists(parsed_gbk_path)) {
    cat("Reading the GBK file \n")
    gbk = BacGWES::parse_genbank_file(gbk_path = gbk_path, g = snp.dat$g) # will return 1 if fails
    saveRDS(gbk, parsed_gbk_path)
  } else {
    cat("Loading parsed gbk file \n")
    gbk = readRDS(parsed_gbk_path)
  }

  # BLK3
  cat("\n\n #################### BLOCK 3 #################### \n\n")
  if(!file.exists(cds_var_path)) {
    cat("Estimating the variation in CDS \n")
    cds_var = BacGWES::estimate_variation_in_CDS(gbk = gbk, snp.dat = snp.dat, ncores = ncores, num_clusts_CDS = num_clusts_CDS, clust_plt_path = clust_plt_path)
    saveRDS(cds_var, cds_var_path)
  } else {
    cat("Loading previous CDS variation estimates \n")
    cds_var = readRDS(cds_var_path)
  }

  # BLK4
  cat("\n\n #################### BLOCK 4 #################### \n\n")
  if(!file.exists(hdw_path)) {
    cat("Estimating per sequence Hamming distance \n")
    # hdw = perform_pop_struct_correction_sparse(snp.matrix = snp.dat$snp.matrix, nsnp = snp.dat$nsnp)
    hdw = BacGWES::estimate_Hamming_distance_weights(snp.dat = snp.dat)
    saveRDS(hdw, hdw_path)
  } else {
    cat("Loading previous Hamming distance estimates \n")
    hdw = readRDS(hdw_path)
  }

  # BLK5
  cat("\n\n #################### BLOCK 5 #################### \n\n")
  if(file.exists(lr_save_path) & file.exists(sr_save_path)) {
    cat("Loading previous MI computation \n")
    sr_links = read.table(sr_save_path)
    colnames(sr_links) = c("clust_c", "pos1", "pos2", "clust1", "clust2", "len", "MI", "srp_max", "ARACNE")
  } else {

    cat("Commencing MI computation \n")
    sr_links = BacGWES::perform_MI_computation(snp.dat = snp.dat, hdw = hdw, cds_var = cds_var, ncores = ncores,
                                               lr_save_path = lr_save_path, sr_save_path = sr_save_path,
                                               plt_folder = dset, sr_dist = sr_dist, discard_MI_threshold_lr = discard_MI_threshold_lr,
                                               max_blk_sz = max_blk_sz, srp_cutoff = srp_cutoff, runARACNE = T)
  }

  lr_links = read.table(lr_save_path) # This is written as a tsv file, need to load for plotting
  colnames(lr_links) = c("pos1", "pos2", "c1", "c2", "len", "MI")

  # BLK6
  cat("\n\n #################### BLOCK 6 #################### \n\n")
  BacGWES::make_gwes_plots(lr_links = lr_links, sr_links = sr_links, plt_folder = dset)

  # BLK7
  cat("\n\n #################### BLOCK 7 #################### \n\n")
  if(SnpEff_Annotate == F){
    cat("\n\n ** All done ** \n")
    return(0)
  }

  # Additional paths if annotations are requested
  # tanglegram
  tanglegram_path = file.path(dset, "Tanglegram")
  if(!file.exists(tanglegram_path)) dir.create(tanglegram_path)
  # GWESExplorer
  gwesexplorer_path = file.path(dset, "GWESExplorer")
  if(!file.exists(gwesexplorer_path)) dir.create(gwesexplorer_path)

  if(!file.exists(tophits_path)){
    tophits = BacGWES::perform_snpEff_annotations(dset_name = dset, annotation_folder = file.path(getwd(), dset),
                                                  snpeff_jar = snpeff_jar_path, gbk = gbk, gbk_path = gbk_path,
                                                  cds_var = cds_var, sr_links = sr_links, snp.dat = snp.dat,
                                                  tophits_path = tophits_path, max_tophits = max_tophits)
  } else {
    cat("Loading previous top hits \n")
    tophits = read.table(tophits_path, sep = '\t', header = T)
  }

  # BLK8
  cat("\n\n #################### BLOCK 7 #################### \n\n")
  BacGWES::create_tanglegram(srlinks_tophits = tophits, gbk = gbk, tanglegram_folder = tanglegram_path, break_segments = tanglegram_break_segments)

  # BLK9
  cat("\n\n #################### BLOCK 7 #################### \n\n")
  BacGWES::write_output_for_gwes_explorer(snp.dat = snp.dat, srlinks_tophits = tophits, gwes_explorer_folder = gwesexplorer_path)

  cat("\n\n ** All done ** \n")

}
