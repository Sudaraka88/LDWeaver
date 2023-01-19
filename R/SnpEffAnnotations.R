#' perform_snpEff_annotations
#'
#' Function to generate long and short range GWES plots
#'
#' @importFrom utils write.table
#'
#' @param plt_folder specify the folder to save generated plots (default = NULL, will be saved to a folder called PLOTS in getwd())
#' @param lr_links data.frame containing lr_links or the path to saved lr_links TSV file. Usually output from perform_MI_computation()
#' @param sr_links data.frame containing sr_links or the path to saved sr_links TSV file. Usually output from perform_MI_computation()
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
perform_snpEff_annotations = function(dset_name, annotation_folder, snpeff_jar, gbk, gbk_path, cds_var, sr_links, snp.dat){
  genome_name = gbk@genes@seqinfo@genome
  snpeff_ready = prep_snpEff(RUN_SNPEFF = TRUE, dset = dset_name, genome_name = genome_name, snpeff_jar = snpeff_jar, work_dir = annotation_folder, gbk_path = gbk_path)

  if(snpeff_ready==1){ # ready for annotation?
    # Write the interested sites to a VCF file
    snps_to_ann = sort(unique(c(sr_links$pos1, sr_links$pos2))) # These are the SNPs that need annotation
    snps_to_ann_idx = sapply(snps_to_ann, function(x) which(snp.dat$POS %in% x)) # These are their original positions

    vcf_sr_path = file.path(annotation_folder, "sr_snps.vcf")

    append_vcf_header(vcf_sr_path, snp.dat$g) # This creates a fresh VCF file
    create_vcf_file(genome_name, snps_to_ann, cds_var$ref[snps_to_ann_idx], cds_var$alt[snps_to_ann_idx], vcf_sr_path)

    # write.table(x = vcf, file = vcf_sr_path, append = T, row.names = F, col.names = F, quote = F, sep = '\t')
  } else {
    # Return with snpeff error
    warning("snpEff error!")
    return(-1)
  }
  vcf_annotated_path = file.path(annotation_folder, "sr_snps_ann.vcf")
  stat_sr_path = file.path(annotation_folder, "sr_annotated_stats.html")

  snpeff_complete = run_snpeff(dset = dset_name, genome_name = genome_name, snpeff_jar = snpeff_jar, work_dir = annotation_folder,
                               stat_sr_path = stat_sr_path, vcf_sr_path = vcf_sr_path, vcf_annotated_path = vcf_annotated_path)

}


prep_snpEff = function(RUN_SNPEFF = TRUE, dset, genome_name, snpeff_jar, work_dir, gbk_path){
  # java -jar snpEff.jar build -genbank -v Enterococcus_faecalis_gcf_000007785.1
  # java -Xmx16G -jar snpEff.jar -v -stats efaecalis_link_annotated.html Enterococcus_faecalis_gcf_000007785.1 ../bac_efaecalis/efaecalis_link_snps.vcf > efaecalis_link_annotated.vcf
  # Enterococcus_faecalis_gcf_000007785.1.genome : Enterococcus_faecalis_gcf_000007785.1
  # Enterococcus_faecalis_gcf_000007785.1.NC_004668.1.codonTable : Bacterial_and_Plant_Plastid

  if(!RUN_SNPEFF) return(0) else status = 1 # run snpeff?
  print(paste("Setting up snpeff for:", dset))
  t0 = Sys.time()

  if(!file.exists(snpeff_jar)){ # does it point to <snpEff.jar>
    warning(paste("<snpEff.jar> not available at:", snpeff_jar))
    status = -1
  }

  options(warn = -1)
  jout = system(paste('java -jar', snpeff_jar), ignore.stderr = T)
  options(warn = 0)

  if(jout!=255){ # file missing warning must be given by snpEff, there must be a better way to check this?
    warning(paste(snpeff_jar, "not functional, please confirm that it is configured correctly..."))
    status = -1
  }

  if(status != -1){
    # move files around as required for snpEff
    snpeff_config = file.path(work_dir, "snpEff.config")
    #TODO: How best to place the config file within an R package?

    if(file.exists(snpeff_config)) unlink(snpeff_config)
    # file.copy("snpEff.template", snpeff_config)
    snpeff_template_path = system.file("extdata", "snpEff.template", package = "BacGWES")
    file.copy(snpeff_template_path, snpeff_config)

    write.table(x = paste(dset, '.genome : ', dset, sep = ""), file = snpeff_config, append = T, col.names = F, row.names = F, quote = F)
    write.table(x = paste(dset, '.', genome_name, '.codonTable : Bacterial_and_Plant_Plastid', sep = ""), file = snpeff_config, append = T, col.names = F, row.names = F, quote = F)

    snpeff_data = file.path(work_dir, "snpEff_data")
    if(file.exists(snpeff_data))  unlink(snpeff_data, recursive = T)
    dir.create(snpeff_data) # create the data file
    dir.create(file.path(snpeff_data, dset))
    file.copy(gbk_path, file.path(snpeff_data, dset, "genes.gbk"))


    system(paste('java -jar', snpeff_jar, 'build -genbank -config', snpeff_config, '-dataDir', snpeff_data, '-v', dset)) # Build index for snpEff
  }
  print(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
  return(status)
}

append_vcf_header = function(path, g){
  write("##fileformat=VCF4.1", file = path, append = F)
  write(paste("##contig=<ID=1,length=",g,">", sep = ""), file = path, append = T)
  write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file = path, append = T)
  write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO", file = path, append = T)
}

create_vcf_file = function(genome_name, snps_to_ann, REF, ALT, vcf_sr_path){
  vcf = data.frame(CHROM = genome_name,
                   POS = sprintf("%.0f", snps_to_ann),
                   ID = '.',
                   REF = REF,
                   ALT = ALT,
                   QUAL = '.',
                   FILTER = '.',
                   INFO = '.')
  write.table(x = vcf, file = vcf_sr_path, append = T, row.names = F, col.names = F, quote = F, sep = '\t')
}

# annotate with snp_eff
run_snpeff = function(RUN_SNPEFF = TRUE, dset = NULL, genome_name = NULL, snpeff_jar = NULL, work_dir = NULL, stat_sr_path = NULL, vcf_sr_path = NULL, vcf_annotated_path = NULL){
  #java -Xmx16G -jar snpEff.jar -v -stats Aus0004_efcm_link_annotated.html Enterococcus_faecium_gcf_000250945.1 ../bac_Efaecium/efcm_link_snps.vcf > Aus0004_efcm_link_annotated.vcf
  snpeff_config = file.path(work_dir, "snpEff.config")
  snpeff_data = file.path(work_dir, "snpEff_data")
  system(paste('java -Xmx16G -jar', snpeff_jar, '-v -dataDir', snpeff_data, '-config', snpeff_config,'-stats', stat_sr_path, dset, vcf_sr_path, '>', vcf_annotated_path))
  status = 1
  return(status)
}
