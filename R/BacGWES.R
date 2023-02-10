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
#' @param check_gbk_fasta_lengths check if the gbk reference sequence length matches the with fasta alignment (default = T)
#' @param snp_filt_method specify the filtering method for SNP extraction: 'relaxed' or 'default' (default = 'default')
#' @param gap_freq sites with a gap frequency >gap_greq will be dropped (default = 0.15)
#' @param maf_freq sites with a minor allele frequency <maf_freq will be dropped (default = 0.01)
#' @param snpeff_jar_path path to <snpEff.jar>. If unavailable or not required, set SnpEff_Annotate = F
#' @param hdw_threshold Hamming distance similarity threshold (default = 0.1, i.e. 10\%) - add more?
#' @param perform_SR_analysis_only specify whether to only perform the short range link analysis (default = FALSE)
#' @param SnpEff_Annotate specify whether to perform annotations using SnpEff
#' @param sr_dist links less than <sr_dist> apart are considered 'short range' (default = 20000), range 1000 - 25000 bp.
#' @param lr_retain_level specify the long-range MI retaining percentile (default = 0.99) - in each block, only the top 1\% of lr MI links will be retained, range 0.1-0.999. WARNING! Low values will results in a HUGE lr_links.tsv file!
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
#'
#' @return Numeric Value 0 if successful (all generated outputs will be saved)
#'
#' @examples
#' \dontrun{
#' sr_links_red = BacGWES(dset = "efcm", aln_path = "<efcm_aln>", gbk_path = "<efcm.gbk>")
#' }
#' @export
BacGWES = function(dset, aln_path, gbk_path, check_gbk_fasta_lengths = T, snp_filt_method = "default",
                   gap_freq = 0.15, maf_freq = 0.01, snpeff_jar_path = NULL, hdw_threshold = 0.1,
                   perform_SR_analysis_only = F, SnpEff_Annotate = T, sr_dist = 20000, lr_retain_level = 0.99,
                   max_tophits = 250, num_clusts_CDS = 3, srp_cutoff = 3, tanglegram_break_segments = 5,
                   multicore = T, max_blk_sz = 10000, ncores = NULL){
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

<<<<<<< HEAD
  #TODO: Change lr_retain_level to save #n links
=======
  # Welcome message
  timestamp()
  cat(paste("\n\n Performing GWES analysis on:", dset, "\n\n"))
  #TODO: show a set of parameters before starting the run (save paths, file sizes, etc.)
>>>>>>> 20b1541ce6952880e14a18cfb521005767805046

  # Sanity checks
  # annotations
  if(SnpEff_Annotate == T) {
    if(is.null(snpeff_jar_path)) stop("You must specify <snpeff_jar_path> for annotations. To run without annotations, set SnpEff_Annotate = F")
    if(!file.exists(snpeff_jar_path)) stop(paste("<SnpEff.jar> not found at:", snpeff_jar_path, "please check the path provided"))
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

  if(lr_retain_level <= 0.1 | lr_retain_level >= 0.9999) {
    warning(paste("Unable to use the provided value for <lr_retain_level>, using", 0.99))
    lr_retain_level = 0.99
  }

  if(lr_retain_level < 0.8) warning("The given lr_retain_level can results in a large lr_links.tsv file!")

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

  # normalise_input_paths
  aln_path = normalizePath(aln_path)
  gbk_path = normalizePath(gbk_path)
  if(!is.null(snpeff_jar_path)) snpeff_jar_path = normalizePath(snpeff_jar_path)


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


  ######## Welcome message ########
  {
    timestamp()
    if(ncores > 1) cat(paste("\n\n Performing GWES analysis on:", dset, " - using", ncores, "cores\n\n"))
    if(ncores == 1) cat(paste("\n\n Performing GWES analysis on:", dset, "\n\n"))
    if(perform_SR_analysis_only) cat("Only short-range analysis requested. \n")
    cat(paste("All outputs will be saved to:", normalizePath(dset), "\n"))
    cat(paste("\n *** Input paths *** \n\n"))
    cat(paste("* Alignment:", aln_path, "\n"))
    cat(paste("* GenBank Annotation:", gbk_path, "\n"))
    if(!is.null(snpeff_jar_path)) cat(paste("* SnpEff Annotations will be performed on short-range links. SnpEff path:", snpeff_jar_path, "\n"))

    cat(paste("\n *** Parameters *** \n\n"))

    if(snp_filt_method == "default") {
      cat(paste("Default SNP filtering: sites with gap_freq <", gap_freq, "and non-gap minor allele freq >", maf_freq, "will be retained. \n"))
    } else {
      cat(paste("Relaxed SNP filtering: sites with gap_freq <", gap_freq, "and minor allele freq >", maf_freq, "will be retained. \n"))
    }
    cat(paste("Hamming distance calculation weight:", hdw_threshold, "\n"))
    cat(paste("Links <=", sr_dist, "bp-apart will be classified as short-range (sr-links) \n"))
    if(!perform_SR_analysis_only) cat(paste("Approx. top", lr_retain_level, "long range links will be saved \n"))
    cat(paste("Top sr-links with -log10(p) >", srp_cutoff, "will be saved \n"))
    cat(paste("Tanglegram/GWESExplorer outputs will illustrate upto:", max_tophits, "top sr-links \n"))
    cat(paste("MI Computation will use a max block size of:", max_blk_sz, "x", max_blk_sz, "SNPs! Reduce <max_blk_sz> if RAM is scarce!\n\n"))
    cat(paste("~~~~~ https://github.com/Sudaraka88/BacGWES/ ~~~~~"))
  }
  ######## <Welcome message> ########
  t_global = Sys.time() # global timer

  # BLK1
  cat("\n\n #################### BLOCK 1 #################### \n\n")
  if(!file.exists(ACGTN_snp_path)) {
    t0 = Sys.time()
    cat(paste("Parsing Alignment:", aln_path ,"\n"))
    snp.dat = BacGWES::parse_fasta_alignment(aln_path = aln_path, method = snp_filt_method, gap_freq = gap_freq, maf_freq = maf_freq)
    cat("Step 5: Savings snp.dat...")
    saveRDS(snp.dat, ACGTN_snp_path)
    cat(paste("BLOCK 1 complete in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))
  }else{
    cat("Loading previous snp matrix \n")
    snp.dat = readRDS(ACGTN_snp_path)
  }

  # BLK2
  cat("\n\n #################### BLOCK 2 #################### \n\n")
  if(!file.exists(parsed_gbk_path)) {
    cat("Reading the GBK file \n")
    gbk = BacGWES::parse_genbank_file(gbk_path = gbk_path, g = snp.dat$g, length_check = check_gbk_fasta_lengths) # will return 1 if fails
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
    hdw = BacGWES::estimate_Hamming_distance_weights(snp.dat = snp.dat, threshold = hdw_threshold)
    saveRDS(hdw, hdw_path)
  } else {
    cat("Loading previous Hamming distance estimates \n")
    hdw = readRDS(hdw_path)
  }

  # BLK5
  cat("\n\n #################### BLOCK 5 #################### \n\n")
  if( (perform_SR_analysis_only &  file.exists(sr_save_path)) |   (file.exists(lr_save_path) &  file.exists(sr_save_path)) ) {
    cat("Loading previous MI computation \n")
    sr_links = read.table(sr_save_path)
    colnames(sr_links) = c("clust_c", "pos1", "pos2", "clust1", "clust2", "len", "MI", "srp_max", "ARACNE")
  } else {

    cat("Commencing MI computation \n")
    sr_links = BacGWES::perform_MI_computation(snp.dat = snp.dat, hdw = hdw, cds_var = cds_var, ncores = ncores,
                                               lr_save_path = lr_save_path, sr_save_path = sr_save_path,
                                               plt_folder = dset, sr_dist = sr_dist, lr_retain_level = lr_retain_level,
                                               max_blk_sz = max_blk_sz, srp_cutoff = srp_cutoff, runARACNE = T,
                                               perform_SR_analysis_only = perform_SR_analysis_only)
  }

  if(file.exists(lr_save_path)){
    lr_links = read.table(lr_save_path) # This is written as a tsv file, need to load for plotting
    colnames(lr_links) = c("pos1", "pos2", "c1", "c2", "len", "MI")
  } else {
    lr_links = NULL # in case lr_links are not available
  }

  if(nrow(sr_links) == 0){
    stop("No potentially important sr_links were identified! Cannot continue analysis...")
  }


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

  cat(paste("\n\n ** All done in", round(difftime(Sys.time(), t_global, units = "mins"), 3), "m ** \n"))

}
