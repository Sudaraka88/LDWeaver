#' LDWeaver
#'
#' Function to run the LDWeaver pipeline
#'
#' @importFrom parallel detectCores
#' @importFrom utils timestamp packageVersion
#'
#' @param dset name of the dataset, all outputs will be saved to the folder <dset>
#' @param aln_path path to the multi fasta alignment - can be a SNP only or full alignment
#' @param aln_has_all_bases specify whether the alignment has all bases in the reference genome (default = T). For a SNP only alignment, set to F and provide SNP positions in pos (see example below)
#' @param pos numeric vector of positions for each base in the alignment (default = NULL). Only required for SNP alignments (e.g. Output from snp-sites or FastaR https://github.com/Sudaraka88/FastaR)
#' @param gbk_path path to genbank annotations file (default = NULL). Only provide one of genbank or gff3
#' @param gff3_path path to gff3 annotations file (default = NULL). Only provide one of genbank or gff3
#' @param ref_fasta_path path to Reference fasta file. This must be in fasta format and contain exactly one sequence! Only required for gff3 annotations
#' @param validate_ref_ann_lengths check if the reference sequence length matches with fasta alignment. (default = T for full alignments, = F for snp-only alignments.). Set to F for alignments that only contain a portion of the genome
#' @param snp_filt_method specify SNP filtering method: 'relaxed' or 'default' (default = 'default'). Unlike default, relaxed considers ambiguous/gap characters (N) as minor alleles when applying the
#' maf_freq filter. Eg: Under default filter values, a site with allele frequencies A:0.85, C:0.0095, N:0.1405 will be respectively, dropped and allowed by 'default' and 'relaxed'
#' @param gap_freq sites with a gap frequency >gap_greq will be dropped (default = 0.15)
#' @param maf_freq sites with a minor allele frequency <maf_freq will be dropped (default = 0.01)
#' @param hdw_threshold Hamming distance similarity threshold (default = 0.1, i.e. 10\%) - lower values will force stricter population structure control at the cost of masking real signal
#' @param perform_SR_analysis_only specify whether to only perform LDWeaver (short-range link) analysis. (default = FALSE) will analyse all links (i.e., LDWeaver + SpydrPick)
#' @param SnpEff_Annotate specify whether to perform annotations using SnpEff (default = TRUE)
#' @param sr_dist links shorter than <sr_dist> apart are considered 'short range' (default = 20000), range (1000 - 100000) bp
#' @param lr_retain_links specify the maximum number of long-range MI links to retain (default = 1000000) - in each block, only a top subset of links will be saved such approximately this many links will be retained
#' @param max_tophits specify the maximum number of short range links to save as <sr_tophits.tsv>. Note: all short-range links will be annotated (and saved separately),
#' but only the top <max_tophits> will be used for visualisation (default = 250), range 50 - 1000
#' @param num_clusts_CDS partition to genome into num_clusts_CDS regions using k-means (default = 3) range 1 - 10 - accounts for within genome heterogeneity in the short-range analysis
#' @param srp_cutoff specify the short-range -log10(p) cut-off value to discard short-range links before returning the data.frame (default = 3, i.e., p = 10^-3). This setting has no impact on modelling since all links are used.
#' Larger values will reduce memory usage, plotting time and ARACNE run time. If all links are required, set to 0 (i.e. p = 10^0 = 1), range 0 - 5
#' @param tanglegram_break_segments specify the number of genome tanglegram segments to prepare (default = 5 segments). Set NULL to skip tanglegram. range 1 - 10.
#' @param write_gwesExplorer specify whether to write output for GWESExplorer (default = T)
#' @param multicore specify whether to use parallel processing (default = T)
#' @param ncores specify the number of cores to use for parallel processing. Auto detect (detect = NULL)
#' @param max_blk_sz specify maximum block size for MI computation (default = 10000), larger sizes require more RAM, range 1000 - 100000
#' @param save_additional_outputs specify whether to save outputs such as extracted SNPs and Hamming distance weights. Recommended for very large datasets to save time on re-computation (default = F)
#' @param mega_dset specify whether the datasets is megascale. This mode requires spam and spam64 packages. This is  >5 times slower, set to TRUE only if the normal analysis fails (default = F)

#'
#' @return All generated outputs will be saved to folder <dset>.
#'
#' @examples
#' \dontrun{
#' # These examples use data included with the package.
#'
#' # Example 1 - using a full alignment that includes non SNP sites
#' dset <- "full_dset"
#' gbk_path <- system.file("extdata", "sample.gbk", package = "LDWeaver")
#' aln_path <- system.file("extdata", "sample.aln.gz", package = "LDWeaver")
#' LDWeaver::LDWeaver(dset = dset,  aln_path = aln_path,  gbk_path = gbk_path,
#'                    validate_ref_ann_lengths = F)
#'
#' # Example 2- using a SNP only alignment
#' dset <- "snp_dset"
#' gbk_path <- system.file("extdata", "sample.gbk", package = "LDWeaver")
#' aln_path <- system.file("extdata", "snp_sample.fa.gz", package = "LDWeaver")
#' pos <- as.numeric(readLines(system.file("extdata", "snp_sample.fa.pos",
#'                   package = "LDWeaver")))
#' LDWeaver::LDWeaver(dset = dset,  aln_path = aln_path, aln_has_all_bases = F,
#'                   pos = pos, gbk_path = gbk_path)
#'
# Example 3 - Redoing the full analysis as a mega scale dataset
#' dset <- "full_dset_spam"
#' gbk_path <- system.file("extdata", "sample.gbk", package = "LDWeaver")
#' aln_path <- system.file("extdata", "sample.aln.gz", package = "LDWeaver")
#' LDWeaver::LDWeaver(dset = dset,  aln_path = aln_path,  gbk_path = gbk_path,
#'                    validate_ref_ann_lengths = F, mega_dset = T)
#' }
#' @export
LDWeaver = function(dset, aln_path, aln_has_all_bases = T, pos = NULL, gbk_path = NULL, gff3_path = NULL,
                    ref_fasta_path = NULL, validate_ref_ann_lengths = T, snp_filt_method = "default",
                    gap_freq = 0.15, maf_freq = 0.01, hdw_threshold = 0.1, perform_SR_analysis_only = F,
                    SnpEff_Annotate = T, sr_dist = 20000, lr_retain_links = 1e6,
                    max_tophits = 250, num_clusts_CDS = 3, srp_cutoff = 3, tanglegram_break_segments = 5,
                    write_gwesExplorer = T, multicore = T, max_blk_sz = 10000, ncores = NULL,
                    save_additional_outputs = F, mega_dset = F){

  # Build blocks
  # BLK1: Extract SNPs and create sparse Mx from MSA (fasta)
  # BLK2: Parse GBK or GFF+REF
  # BLK3: Estimate diversity within each CDS, cluster and paint < # possible inputs on methods>
  # BLK4: Compute Hamming Distance weights
  # BLK5: Compute MI between all links, sr_links model fitter, ARACNE
  # BLK6: Genomewide LD Map
  # BLK7: GWES_plots
  # BLK8: Snpeff annotation pipeline, determine tophits
  # BLK9: Tanglegram (depends: chromoMap)
  # BLK10: GWESExplorer (depends: GWESExplorer)
  # BLK11: Cleanup

  #TODO: Provide the option to skip SNP extraction and use the whole provided alignment (redundant if pre-filtered)
  #TODO: Add the option to provide GFF file without reference sequence
  #TODO: Count through blocks and automate the displayed BLOCK NUMBER
  #TODO: Add Hamming Distance plot, can we have a SNP Tree + Hamming Distance weights to show population structure control?
  #TODO: Benchmark spam VS matrixExtra

  #NOTE: SnpEff does not parse the GBK and GFF3 file from the same refseq reference genome the same way. There might be differences between annotations/tophits/etc.
  # # Welcome message # #

  # Sanity checks
  # annotations
  if( (is.null(gbk_path) & is.null(gff3_path)) | (!is.null(gbk_path) & !is.null(gff3_path)) ) stop("Either gbk_path or gff3_path must be provided") # only one of gbk or gff can be NULL
  if(!is.null(gff3_path) & is.null(ref_fasta_path)) stop("Reference fasta file must be provided for gff3 annoations") # only one of gbk or gff can be NULL

  if(SnpEff_Annotate == T) {
    # Added snpEff to inst/extdata
    snpeff_jar_path = system.file("extdata", "snpEff.jar", package = "LDWeaver")
    ######################################## These checks must be unnecessary now ########################################
    # if(is.null(snpeff_jar_path)) stop("<snpeff_jar_path> must be provided for annotations. To run without annotations, set SnpEff_Annotate = F")
    # if(!file.exists(snpeff_jar_path)) stop(paste("<SnpEff.jar> not found at:", snpeff_jar_path, "please check the path provided"))
    ######################################################################################################################
    order_links = F # sr_links should be ordered at the end after adding annotations
  } else {
    snpeff_jar_path = NULL
    order_links = T # sr_links will be ordered and saved without annotations
  }

  # alignment (now with SNP-only alignment support 2023/10/12)
  if(aln_has_all_bases == F){ # snp-only alignment, POS must be provided
    # sanity checks
    if(is.null(pos)) stop("A numeric vector of 'positions' <pos> must be provided if aln_has_all_bases = F")
    if(!is.numeric(pos)) stop("Provided pos must be numeric!")
    if(any(duplicated(pos))) stop("Provided pos contains duplicates!")
  } else { # This is a full alignment, pos must be NULL
    if(!is.null(pos)) stop("pos cannot be provided for alignments with all bases! Depending on the use case, either set pos = NULL or aln_has_all_bases = T")
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
  if(sr_dist < 1000 | sr_dist > 100000) {
    tmp = get_sr_threshold(sr_dist)
    warning(paste("Unable to use the provided value for <sr_dist>:", sr_dist, ", instead using", tmp, "\n"))
    sr_dist = tmp
  }

  if(lr_retain_links <= 1e3 | lr_retain_links >= 1e10) {
    warning(paste("Unable to use the provided value for <lr_retain_links>, using", 1e6))
    lr_retain_links = 1e6
  }

  if(lr_retain_links > 1e6) warning("The given lr_retain_links value may generate a very large lr_links.tsv file!")

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

  if(!is.null(tanglegram_break_segments)){
    if(tanglegram_break_segments < 0 | tanglegram_break_segments > 10) {
      warning(paste("Unable to use the provided value for <tanglegram_break_segments>, using", 5))
      tanglegram_break_segments = 5
    }
  }
  if(max_blk_sz < 1000 | max_blk_sz > 100000) {
    warning(paste("Unable to use the provided value for <max_blk_sz>, using", 10000, "...!If this value is causing the function to crash, consider reducing!..."))
    max_blk_sz = 10000
  }

  if(aln_has_all_bases == F){ # For snp-only alignments, reference and alignment length checks will fail, stop checking
    validate_ref_ann_lengths = F
  }

  if(mega_dset){ # This is a mega dataset, we need spam and spam64 packages
    if(!requireNamespace("spam") & !requireNamespace("spam64")){
      message("mega_dset is set to TRUE but spam and spam64 dependencies are not installed.")
      return(invisible())
    }
    message("mega_dset is selected. Warning! This mode has a much slower run time. Setting spam.force64=TRUE (see https://cran.r-project.org/web/packages/spam64/spam64.pdf)")
    options(spam.force64 = TRUE)
  }


  ######################################################################################################################
  # normalise_input_paths
  aln_path = normalizePath(aln_path)
  if(!is.null(gbk_path)) gbk_path = normalizePath(gbk_path)
  if(!is.null(gff3_path)) gff3_path = normalizePath(gff3_path)
  if(!is.null(ref_fasta_path)) ref_fasta_path = normalizePath(ref_fasta_path)
  if(!is.null(snpeff_jar_path)) snpeff_jar_path = normalizePath(snpeff_jar_path)


  # setup paths
  if(!file.exists(dset)) dir.create(dset) # save everything in here

  # Save console output as a text file
  info_file = file.path(dset, paste("LDW_run_",format(Sys.time(), "%Y%m%d%H%M%S"), ".txt", sep = ""))
  suppressWarnings(sink(file= NULL))
  sink(info_file, split = T)

  add_path = file.path(dset, "Additional_Outputs") # Additional Outputs
  if(save_additional_outputs) {
    if(!file.exists(add_path)) dir.create(file.path(dset, "Additional_Outputs"))
  }

  # default paths (will be moved at the end of the run)
  ACGTN_snp_path = file.path(dset, "snp_ACGTN.rds")
  if(!is.null(gbk_path)) parsed_gbk_path = file.path(dset, "parsed_gbk.rds")
  if(!is.null(gff3_path)) parsed_gff_path = file.path(dset, "parsed_gff3.rds")
  cds_var_path = file.path(dset, "cds_var.rds")
  hdw_path = file.path(dset, "hdw.rds")
  gwLDplt_path = file.path(dset, "LD_plot.png")

  # a previous computation might exist
  if(file.exists(file.path(add_path, "snp_ACGTN.rds"))) ACGTN_snp_path = file.path(add_path, "snp_ACGTN.rds")
  if(!is.null(gbk_path)) if(file.exists(file.path(add_path, "parsed_gbk.rds"))) parsed_gbk_path = file.path(add_path, "parsed_gbk.rds")
  if(!is.null(gff3_path)) if(file.exists(file.path(add_path, "parsed_gff3.rds"))) parsed_gff_path = file.path(add_path, "parsed_gff3.rds")
  if(file.exists(file.path(add_path, "cds_var.rds"))) cds_var_path = file.path(add_path, "cds_var.rds")
  if(file.exists(file.path(add_path, "hdw.rds"))) hdw_path = file.path(add_path, "hdw.rds")

  clust_plt_path = file.path(dset, "CDS_clustering.png")

  # These files can exist from a previous run, with or without cleaning, choose the correct path
  lr_save_path = file.path(dset, "lr_links.tsv")
  sr_save_path = file.path(dset, "sr_links.tsv")
  tophits_path = file.path(dset, "sr_tophits.tsv") # This is for sr only

  if(file.exists(file.path(dset, "Temp/lr_links.tsv"))) lr_save_path = file.path(dset, "Temp/lr_links.tsv")
  if(file.exists(file.path(dset, "Temp/sr_links.tsv"))) sr_save_path = file.path(dset, "Temp/sr_links.tsv")
  if(file.exists(file.path(dset, "Tophits/sr_tophits.tsv"))) tophits_path = file.path(dset, "Tophits/sr_tophits.tsv") # This is for sr only

  ######## Welcome message ########
  {
    timestamp()
    cat("\n ***** This is LDWeaver", as.character(packageVersion(pkg = "LDWeaver")), " *****\n")
    test_openmp()
    if(ncores > 1) cat(paste("\n\n Performing GWES analysis on:", dset, " - using", ncores, "cores\n\n"))
    if(ncores == 1) cat(paste("\n\n Performing GWES analysis on:", dset, "\n\n"))
    if(perform_SR_analysis_only) cat("Only short-range analysis requested. \n")
    cat(paste("All outputs will be saved to:", normalizePath(dset), "\n"))
    cat(paste("\n *** Input paths *** \n\n"))
    if(mega_dset) cat(paste("* Mega Alignment:", aln_path, "\n")) else cat(paste("* Alignment:", aln_path, "\n"))
    if(!is.null(gbk_path)) {
      cat(paste("* GenBank Annotation:", gbk_path, "\n"))
      cat(paste("* Parser built using genbankr source (https://github.com/gmbecker/genbankr) \n"))
    }
    if(!is.null(gff3_path)) cat(paste("* GFF3 Annotation:", gff3_path, "\n"))
    if(!is.null(snpeff_jar_path)) cat(paste("* SnpEff Annotations will be performed on short-range links. SnpEff path:", snpeff_jar_path, "\n"))

    cat(paste("\n *** Parameters *** \n\n"))

    if(snp_filt_method == "default") {
      cat(paste("Default SNP filtering: sites with gap_freq <", gap_freq, "and non-gap minor allele freq >", maf_freq, "will be retained. \n"))
    } else {
      cat(paste("Relaxed SNP filtering: sites with gap_freq <", gap_freq, "and minor allele freq >", maf_freq, "will be retained. \n"))
    }
    cat(paste("Hamming distance calculation weight:", hdw_threshold, "\n"))
    cat(paste("Links <=", sr_dist, "bp-apart will be classified as short-range (sr-links) \n"))
    if(!perform_SR_analysis_only) cat(paste("Approx. top", lr_retain_links, "long range links will be saved \n"))
    cat(paste("Top sr-links with -log10(p) >", srp_cutoff, "will be saved \n"))
    if(!is.null(tanglegram_break_segments)) cat(paste("Tanglegram/GWESExplorer outputs will illustrate upto:", max_tophits, "top sr-links \n"))
    cat(paste("MI Computation will use a max block size of:", max_blk_sz, "x", max_blk_sz, "SNPs! Reduce <max_blk_sz> if RAM is scarce!\n\n"))
    cat(paste("~~~~~ https://github.com/Sudaraka88/LDWeaver/ ~~~~~"))
  }
  ######## <Welcome message> ########
  t_global = Sys.time() # global timer

  # BLK1
  cat("\n\n #################### BLOCK 1 #################### \n\n")
  if(!file.exists(ACGTN_snp_path)) {
    t0 = Sys.time()
    cat(paste("Parsing Alignment:", aln_path ,"\n"))

    # Adding support for SNP-only alignments
    if(aln_has_all_bases == T){
      snp.dat = LDWeaver::parse_fasta_alignment(aln_path = aln_path, method = snp_filt_method, gap_freq = gap_freq, maf_freq = maf_freq, mega_dset = mega_dset)
      if(save_additional_outputs){
        cat("Step 5: Savings snp.dat...")
        saveRDS(snp.dat, ACGTN_snp_path)
      }
    } else {
      snp.dat = LDWeaver::parse_fasta_SNP_alignment(aln_path = aln_path, pos = pos, method = snp_filt_method, gap_freq = gap_freq, maf_freq = maf_freq, mega_dset = mega_dset)
      # Note that snp.dat$g = NULL (we cannot measure this, need to get it from the genbank file)
      # we cannot save snp.dat here due to absent snp.dat$g, moving downstream (block 2)
    }


    cat(paste("BLOCK 1 complete in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s \n"))
  }else{
    cat("Loading previous snp matrix \n")
    snp.dat = readRDS(ACGTN_snp_path)
  }

  # BLK2
  cat("\n\n #################### BLOCK 2 #################### \n\n")
  if(!is.null(gbk_path)){ # genbank file format
    gff = NULL # alternative set to NULL
    if(!file.exists(parsed_gbk_path)) {
      cat(paste("Reading the GBK file, validate_length_check = ", validate_ref_ann_lengths, "\n"))
      gbk_ = LDWeaver::parse_genbank_file(gbk_path = gbk_path, g = snp.dat$g, length_check = validate_ref_ann_lengths) # will return 1 if fails
      gbk = gbk_$gbk # gbk is coming in gbk_$gbk now
      # For snp-only alignments, g = snp.dat$g = NULL and validate_ref_ann_lengths = F, won't be an issue...
      if(save_additional_outputs){
        saveRDS(gbk, parsed_gbk_path)
      }
    } else {
      cat("Loading parsed gbk file \n")
      gbk = readRDS(parsed_gbk_path)
    }
  }
  if(!is.null(gff3_path)){ # gff3 file format
    gbk = NULL # alternative set to NULL
    if(!file.exists(parsed_gff_path)) {
      cat("Reading the gff3 file \n")
      gff = LDWeaver::parse_gff_file(gff3_path = gff3_path, ref_fasta_path = ref_fasta_path, perform_length_check = validate_ref_ann_lengths) # will stop() if fails!
      # gff will contain the reference sequence within it, we need to insert the ref to gbk if it is provided separately
      if(save_additional_outputs){
        saveRDS(gff, parsed_gff_path)
      }
    } else {
      cat("Loading parsed gff3 file \n")
      gff = readRDS(parsed_gff_path)
    }
  }

  # This is a patch to add snp.dat$g back in case this is a SNP only alignment
  if(is.null(snp.dat$g)){
    if(!is.null(gbk_path)){ # input reference is genbank format
      snp.dat$g = gbk_$ref_g
      cat("Extracted ref genome length", snp.dat$g, "from genbank...")
    }
    if(!is.null(gff3_path)){ # input reference is gff3 format
      snp.dat$g = gff$g
    }
    if(save_additional_outputs){
      cat("saving snp.dat...\n")
      saveRDS(snp.dat, ACGTN_snp_path)
      #### This saved file will have snp.dat$g = NULL ### TODO: Ensure no issues downstream
    }
  }

  # BLK3
  cat("\n\n #################### BLOCK 3 #################### \n\n")
  if(!file.exists(cds_var_path)) {
    cat("Estimating the variation in CDS \n")
    cds_var = LDWeaver::estimate_variation_in_CDS(gbk = gbk, gff = gff, snp.dat = snp.dat, ncores = ncores, num_clusts_CDS = num_clusts_CDS, clust_plt_path = clust_plt_path, mega_dset = mega_dset)
    if(save_additional_outputs){
      saveRDS(cds_var, cds_var_path)
    }
  } else {
    cat("Loading previous CDS variation estimates \n")
    cds_var = readRDS(cds_var_path)
  }

  # BLK4
  cat("\n\n #################### BLOCK 4 #################### \n\n")
  if(!file.exists(hdw_path)) {
    cat("Estimating per sequence Hamming distance \n")
    # hdw = perform_pop_struct_correction_sparse(snp.matrix = snp.dat$snp.matrix, nsnp = snp.dat$nsnp)
    hdw = LDWeaver::estimate_Hamming_distance_weights(snp.dat = snp.dat, threshold = hdw_threshold, mega_dset = mega_dset)
    if(save_additional_outputs){
      saveRDS(hdw, hdw_path)
    }
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
    sr_links = LDWeaver::perform_MI_computation(snp.dat = snp.dat, hdw = hdw, cds_var = cds_var, ncores = ncores,
                                                lr_save_path = lr_save_path, sr_save_path = sr_save_path,
                                                plt_folder = dset, sr_dist = sr_dist, lr_retain_links = lr_retain_links,
                                                max_blk_sz = max_blk_sz, srp_cutoff = srp_cutoff, runARACNE = T,
                                                perform_SR_analysis_only = perform_SR_analysis_only,
                                                order_links = order_links,mega_dset = mega_dset)
  }



  # BLK6
  # there is no way to provide a genomewide LD map if LR analysis is not provided with the pipeline, can plot later with spydrpick links
  if(!perform_SR_analysis_only){
    cat("\n\n #################### BLOCK 6 #################### \n\n")
    LDWeaver::genomewide_LDMap(lr_links_path = lr_save_path, sr_links_path = sr_save_path,
                               plot_title = paste("GW-LD:", dset),
                               plot_save_path = gwLDplt_path, links_from_spydrpick = F)
  } else {
    cat("Genomewide LD map cannot be plotted with only the short_range analysis. If SpydrPick links are avaialble, use LDWeaver::genomewide_LDMap()")
  }


  if(nrow(sr_links) == 0){
    suppressWarnings(sink(file= NULL)) ### output info to text file
    stop("No potentially important sr_links were identified! Cannot continue analysis...")
  }


  # BLK6
  cat("\n\n #################### BLOCK 7 #################### \n\n")
  # order links will be false for analyses that require further analysis, make_gwes_plots() will order the links before plotting
  LDWeaver::make_gwes_plots(lr_links = NULL, sr_links = sr_links, plt_folder = dset, are_srlinks_ordered = order_links)

  # BLK7
  cat("\n\n #################### BLOCK 8 #################### \n\n")
  if(SnpEff_Annotate == F){
    LDWeaver::cleanup(dset = dset, delete_after_moving = F)
    cat(paste("\n\n ** All done in", round(difftime(Sys.time(), t_global, units = "mins"), 3), "m ** \n"))
    return()
  }

  if(!file.exists(tophits_path)){
    tophits = LDWeaver::perform_snpEff_annotations(dset_name = dset, annotation_folder = file.path(getwd(), dset),
                                                   snpeff_jar = snpeff_jar_path, gbk = gbk, gbk_path = gbk_path,
                                                   cds_var = cds_var, gff = gff, links_df = sr_links, snp.dat = snp.dat,
                                                   tophits_path = tophits_path, max_tophits = max_tophits)
  } else {
    cat("Loading previous top hits \n")
    tophits = LDWeaver::read_TopHits(top_hits_path = tophits_path)
  }


  # BLK8
  if(!is.null(tanglegram_break_segments)){
    # tanglegram
    tanglegram_path = file.path(dset, "SR_Tanglegram")
    if(!file.exists(tanglegram_path)) dir.create(tanglegram_path)
    cat("\n\n #################### BLOCK 9 #################### \n\n")
    LDWeaver::create_tanglegram(tophits = tophits, gbk = gbk, gff = gff, tanglegram_folder = tanglegram_path, break_segments = tanglegram_break_segments)
  }
  # BLK9
  if(write_gwesExplorer){
    # GWESExplorer
    gwesexplorer_path = file.path(dset, "SR_GWESExplorer")
    if(!file.exists(gwesexplorer_path)) dir.create(gwesexplorer_path)
    cat("\n\n #################### BLOCK 10 #################### \n\n")
    if(mega_dset) {
      message("GWES Explorer output currently not generated for mega datasets\n")
      } else LDWeaver::write_output_for_gwes_explorer(snp.dat = snp.dat, tophits = tophits, gwes_explorer_folder = gwesexplorer_path)
  }


  # BLK10
  # Additional paths if annotations are requested
  # NetworkPlot
  netplot_path = file.path(dset, "SR_network_plot.png")

  cat("\n\n #################### BLOCK 11 #################### \n\n")
  LDWeaver::create_network(tophits = tophits, netplot_path = netplot_path, plot_title = paste("Networks in short-range tophits for", dset))

  if(!perform_SR_analysis_only){
    # BLK11
    cat("\n\n #################### BLOCK 12 #################### \n\n")
    if(SnpEff_Annotate){
      if( !(  (file.exists(file.path(dset, "lr_tophits.tsv"))) | (file.exists(file.path(dset, "Tophits/lr_tophits.tsv"))) ) ) { # if the annotated_links file exists, no need to run again
        LDWeaver::analyse_long_range_links(dset = dset, lr_links_path =  lr_save_path, sr_links_path = sr_save_path, SnpEff_Annotate = T, snpeff_jar_path = snpeff_jar_path,
                                           gbk_path = gbk_path, gff3_path = gff3_path, snp.dat = snp.dat, cds_var = cds_var, ref_fasta_path = ref_fasta_path,
                                           validate_ref_ann_lengths = validate_ref_ann_lengths, mega_dset = mega_dset)
      } else {
        cat("Results from previous LR anlayis exist!")
      }
    } else {
      if( !(  (file.exists(file.path(dset, "lr_gwes.png"))) | (file.exists(file.path(dset, "GWESPlots/lr_gwes.png"))) ) ) { # if the lr_gwes plot exist, no need to run again
        LDWeaver::analyse_long_range_links(dset = dset, lr_links_path =  lr_save_path, sr_links_path = sr_save_path, SnpEff_Annotate = F, mega_dset = mega_dset)
      } else {
        cat("Results from previous LR anlayis exist!")
      }
    }
  }

  LDWeaver::cleanup(dset = dset, delete_after_moving = F)

  cat(paste("\n\n ** All done in", round(difftime(Sys.time(), t_global, units = "mins"), 3), "m ** \n"))
}
