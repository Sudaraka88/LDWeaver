#' perform_snpEff_annotations
#'
#' Function to generate long and short range GWES plots
#'
#' @importFrom utils write.table
#'
#' @param dset_name name of the dataset
#' @param annotation_folder folder to save annotatations
#' @param snpeff_jar path to <snpEff.jar>
#' @param gbk output from parsing the genbank file using BacGWES::parse_genbank_file()
#' @param gbk_path path to genbank file
#' @param cds_var output from BacGWES::estimate_variation_in_CDS()
#' @param sr_links data.frame containing sr_links output from perform_MI_computation()
#' @param snp.dat output from parsing the multi fasta alignment using BacGWES::parse_fasta_alignment()
#' @param max_tophits specify the number of top hits to return (all links will be annotated and saved to the annotations_folder)
#'
#' @return R data frame with the top short range GWES links with snpEff annotated embedded
#'
#' @examples
#' \dontrun{
#' sr_links_red = perform_snpEff_annotations(dset_name, annotation_folder, snpeff_jar, gbk,
#' gbk_path, cds_var, sr_links, snp.dat, max_tophits = 250)
#' }
#'
#' @export
perform_snpEff_annotations = function(dset_name, annotation_folder, snpeff_jar, gbk, gbk_path, cds_var, sr_links, snp.dat, max_tophits = 250){
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
  annotations_path = file.path(annotation_folder, "annotations.tsv")

  snpeff_complete = run_snpeff(dset = dset_name, genome_name = genome_name, snpeff_jar = snpeff_jar, work_dir = annotation_folder,
                               stat_sr_path = stat_sr_path, vcf_sr_path = vcf_sr_path, vcf_annotated_path = vcf_annotated_path)

  if(snpeff_complete == 1){
    ann = convert_vcfann_to_table(vcf_annotated_path, annotations_path, snps_to_ann_idx, cds_var$allele_table, snp.dat$nseq)
  } else {
    # Return with snpeff error
    warning("snpEff error!")
    return(-1)
  }
  srlinks_annotated_path = file.path(annotation_folder,  "sr_links_annotated.tsv")
  srlinks_annotated = add_annotations_to_links(sr_links_red = sr_links, ann = ann, srlinks_annotated_path = srlinks_annotated_path)

  tophits_path = file.path(annotation_folder, "tophits.tsv")
  top_srlinks_annotated = detect_top_hits(l1_a1_d = srlinks_annotated, max_tophits = max_tophits, tophits_path = tophits_path)
  return(top_srlinks_annotated)
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

convert_vcfann_to_table = function(vcf_annotated_path, annotations_path, snps_to_ann_idx, allele_table, nseq){
  ann = read.table(vcf_annotated_path, quote = "") # we can read the annotated VCF as a table into R
  ann = ann[,c(2,4,5,8)]
  colnames(ann) = c("pos", "REF", "ALT", "ANN")

  ann = cbind(ann, annotation =  sapply(ann$ANN, function(x) unlist(strsplit(x, "\\|"))[2]))
  ann = cbind(ann, description =  sapply(ann$ANN, function(x) paste(unique(unlist(strsplit(x, "\\|"))[c(4,5,10,11)]), collapse = ":") ))
  ann = cbind(ann, cds =  sapply(ann$ANN,   function(x) unlist(strsplit(x, "\\|"))[5]))

  code = rep("ns", nrow(ann)) # start as everything is ns
  code[grep("synonymous_variant", ann$annotation)] = "sy"
  code[grep("stop_retained_variant", ann$annotation)] = "sy"
  code[grep("downstream_gene_variant", ann$annotation)] = "ig"
  code[grep("upstream_gene_variant", ann$annotation)] = "ig"

  ann$code = code
  table(ann$code)

  ann = ann[,-4]

  ann$allele_dist = getAlleleDistribution(ad_ = allele_table[, snps_to_ann_idx], nseq = nseq)

  rownames(ann) = NULL

  write.table(x = ann, file = annotations_path, quote = F, row.names = F, sep = '\t', col.names = T)
  return(ann)
}

getAlleleDistribution = function(ad_, nseq){
  # Get allele distribution for SNPs in ad_ in reportable format
  allele_dist = rep(NA, ncol(ad_))
  for(i in 1:ncol(ad_)){
    ad = sort(ad_[ad_[,i]>0, i], decreasing = T)
    n = names(ad)
    allele_dist[i] = paste(n, ad/nseq, sep = ":", collapse = ", ")
  }
  return(allele_dist)
}

add_annotations_to_links = function(sr_links_red, ann, srlinks_annotated_path){

  pos1_idx = sapply(sr_links_red$pos1, function(x) which(ann$pos %in% x))
  pos2_idx = sapply(sr_links_red$pos2, function(x) which(ann$pos %in% x))

  l1_a1_d = data.frame(pos1 = sr_links_red$pos1,
                       pos2 = sr_links_red$pos2,
                       len = sr_links_red$len,
                       ARACNE = sr_links_red$ARACNE,
                       MI = sr_links_red$MI,
                       srp = sr_links_red$srp_max,
                       pos1_ad = ann$allele_dist[pos1_idx],
                       pos2_ad = ann$allele_dist[pos2_idx],
                       pos1_genreg = ann$cds[pos1_idx],
                       pos2_genreg = ann$cds[pos2_idx],
                       pos1_ann = ann$description[pos1_idx],
                       pos2_ann = ann$description[pos2_idx],
                       links = paste(ann$code[pos1_idx], ann$code[pos2_idx], sep = "X"))

  # saving annotated links
  write.table(l1_a1_d, file = srlinks_annotated_path, sep = '\t', col.names = T, row.names = F, quote = F)
  return(l1_a1_d)
}

detect_top_hits = function(l1_a1_d, max_tophits, tophits_path){
  l1_a1_d = l1_a1_d[which(l1_a1_d$ARACNE == 1),] # ARACNE direct only
  l1_a1_d = l1_a1_d[which(l1_a1_d$links!='syXsy'),] # remove syXsy links
  l1_a1_d = l1_a1_d[which(l1_a1_d$pos1_genreg != l1_a1_d$pos2_genreg), ] # remove links from the same region

  if(nrow(l1_a1_d) > max_tophits) l1_a1_d = l1_a1_d[1:max_tophits, ] # all top hits or max_tophits
  write.table(l1_a1_d, file = tophits_path, sep = '\t', col.names = T, row.names = F, quote = F)

  return(l1_a1_d)

}
