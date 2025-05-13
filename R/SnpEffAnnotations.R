#' perform_snpEff_annotations
#'
#' Function to generate long and short range GWES plots
#'
#' @importFrom utils write.table
#'
#' @param dset_name name of the dataset
#' @param annotation_folder folder to save annotatations
#' @param snpeff_jar path to <snpEff.jar>
#' @param snp.dat output from parsing the multi fasta alignment using LDWeaver::parse_fasta_alignment()
#' @param cds_var output from LDWeaver::estimate_variation_in_CDS()
#' @param links_df data.frame containing links (reduced links can be used to save time) output from perform_MI_computation()
#' @param gbk output from parsing the genbank file using LDWeaver::parse_genbank_file() (default = NULL), only provide either GFF3 or GBK annotation
#' @param gbk_path path to genbank file - plan to omit (default = NULL)
#' @param gff output from parsing the gff3 file using LDWeaver::parse_gff_file() (default = NULL), only provide either GFF3 or GBK annotation
#' @param max_tophits specify the number of top hits to return (all links will be annotated and saved to the annotations_folder)
#' @param tophits_path specify the path to save tophit links (default = NULL). If NULL, will be saved to annotations_folder/tophits.tsv
#' @param links_type specify the links type long-range "LR" or short-range "SR" (default = "SR")
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
perform_snpEff_annotations = function(dset_name, annotation_folder, snpeff_jar, snp.dat, cds_var, links_df, gbk = NULL,
                                      gbk_path = NULL, gff = NULL,  tophits_path = NULL, max_tophits = 250, links_type = "SR"){

  # UPDATE: 20250505: github issue: https://github.com/Sudaraka88/LDWeaver/issues/11#issue-2996489972

  # only one of gbk or gff can be NULL
  if( (is.null(gbk) & is.null(gff)) | (!is.null(gbk) & !is.null(gff)) ) stop("Provide either one of gbk or gff")
  if(!is.null(gbk) & is.null(gbk_path)) stop("GBK path must also be provided")

  if(links_type == "LR"){
    vcf_annotated_path = file.path(annotation_folder, "lr_snps_ann.vcf")
    vcf_write_path = file.path(annotation_folder, "lr_snps.vcf")
    # stat_path = file.path(annotation_folder, "lr_annotated_stats.html")
    annotations_path = file.path(annotation_folder, "lr_annotations.tsv")
    links_annotated_path = file.path(annotation_folder,  "lr_links_annotated.tsv")
    if(is.null(tophits_path)) tophits_path = file.path(annotation_folder, "lr_tophits.tsv")
  } else if(links_type == "SR"){
    vcf_annotated_path = file.path(annotation_folder, "sr_snps_ann.vcf")
    vcf_write_path = file.path(annotation_folder, "sr_snps.vcf")
    # stat_path = file.path(annotation_folder, "sr_annotated_stats.html")
    annotations_path = file.path(annotation_folder, "sr_annotations.tsv")
    links_annotated_path = file.path(annotation_folder,  "sr_links_annotated.tsv")
    if(is.null(tophits_path)) tophits_path = file.path(annotation_folder, "sr_tophits.tsv")
  } else {
    stop("Links type must be LR or SR")
  }

  if(!is.null(gbk)) {
    genome_name = gbk@genes@seqinfo@genome
    snpeff_ready = prep_snpEff(dset = dset_name, genome_name = genome_name, snpeff_jar = snpeff_jar, work_dir = annotation_folder, gbk_path = gbk_path)
  }
  if(!is.null(gff)) {
    genome_name = as.character(gff$gff$seqid[1])
    snpeff_ready = prep_snpEff(dset = dset_name, genome_name = genome_name, snpeff_jar = snpeff_jar, work_dir = annotation_folder, gff_path = gff$gff_path, ref_path = gff$ref_path)
  }


  if(snpeff_ready==1){ # ready for annotation?
    # Write the interested sites to a VCF file
    cat("Preparing VCF files ... ")
    snps_to_ann = sort(unique(c(links_df$pos1, links_df$pos2))) # These are the SNPs that need annotation
    snps_to_ann_idx = sapply(snps_to_ann, function(x) which(snp.dat$POS %in% x)) # These are their original positions

    append_vcf_header(vcf_write_path, snp.dat$g) # This creates a fresh VCF file
    create_vcf_file(genome_name, snps_to_ann, cds_var$ref[snps_to_ann_idx], cds_var$alt[snps_to_ann_idx], vcf_write_path)
    cat(" Done!\n")
  } else {
    # Return with snpeff error
    stop("snpEff error!")
    # return(-1)
  }

  snpeff_complete = run_snpeff(dset = dset_name, genome_name = genome_name, snpeff_jar = snpeff_jar, work_dir = annotation_folder,
                               vcf_write_path = vcf_write_path, vcf_annotated_path = vcf_annotated_path)

  if(snpeff_complete == 1){
    cat("Convert VCF to table ...")
    ann = convert_vcfann_to_table(vcf_annotated_path, annotations_path, snps_to_ann_idx, cds_var$allele_table, snp.dat$nseq)
    cat("Done!\n")
  } else {
    # Return with snpeff error
    stop("snpEff error!")
    # return(-1)
  }
  cat("Adding annotations to links ...\n")

  links_annotated = add_annotations_to_links(links_red = links_df, ann = ann, links_annotated_path = links_annotated_path,
                                             links_type = links_type)



  top_links_annotated = detect_top_hits(l1_a1_d = links_annotated, max_tophits = max_tophits, tophits_path = tophits_path)
  # cat("Done! \n")
  return(top_links_annotated)
}


prep_snpEff = function(dset, genome_name, snpeff_jar, work_dir, gbk_path = NULL, gff_path = NULL, ref_path = NULL){

  if(length(grep("windows", .Platform$OS.type))) isWin = T else isWin = F # We might need more hacks for other OSes

  # only one of gbk or gff can be NULL
  if( (is.null(gbk_path) & is.null(gff_path)) | (!is.null(gbk_path) & !is.null(gff_path)) ) stop("Provide either one of gbk_path or gff_path")

  # java -jar snpEff.jar build -genbank -v Enterococcus_faecalis_gcf_000007785.1
  # java -Xmx16G -jar snpEff.jar -v -stats efaecalis_link_annotated.html Enterococcus_faecalis_gcf_000007785.1 ../bac_efaecalis/efaecalis_link_snps.vcf > efaecalis_link_annotated.vcf
  # Enterococcus_faecalis_gcf_000007785.1.genome : Enterococcus_faecalis_gcf_000007785.1
  # Enterococcus_faecalis_gcf_000007785.1.NC_004668.1.codonTable : Bacterial_and_Plant_Plastid

  cat(paste("Setting up snpeff for:", dset, "\n"))
  t0 = Sys.time()

  if(!file.exists(snpeff_jar)){ # does it point to <snpEff.jar>
    stop(paste("<snpEff.jar> not available at:", snpeff_jar))
  }

  options(warn = -1)
  if(isWin){
    cat("COMMAND :=:=:", paste('java -jar', shQuote(snpeff_jar)), '\n')
    jout = shell(paste('java -jar', shQuote(snpeff_jar)), ignore.stderr = T, translate = T)
  } else {
    cat("COMMAND :=:=:", paste('java -jar', shQuote(snpeff_jar)), '\n')
    jout = system(paste('java -jar', shQuote(snpeff_jar)), ignore.stderr = T)
  }

  options(warn = 0)

  if(jout==1){ # for Error:Unable to access jarfile, return value seems to be 1, check for better way (maybe a file.exists would work?)
    stop(paste(snpeff_jar, "not functional, please confirm that it is configured correctly..."))
  }

  # move files around as required for snpEff
  snpeff_config = file.path(work_dir, "snpEff.config")
  #TODO: How best to place the config file within an R package?

  if(file.exists(snpeff_config)) unlink(snpeff_config) # delete the configuration file
  # file.copy("snpEff.template", snpeff_config)
  snpeff_template_path = system.file("extdata", "snpEff.template", package = "LDWeaver")
  file.copy(snpeff_template_path, snpeff_config) # add a new one

  # add the data for this dataset
  write.table(x = paste(dset, '.genome : ', dset, sep = ""), file = snpeff_config, append = T, col.names = F, row.names = F, quote = F)
  write.table(x = paste(dset, '.', genome_name, '.codonTable : Bacterial_and_Plant_Plastid', sep = ""), file = snpeff_config, append = T, col.names = F, row.names = F, quote = F)


  snpeff_data = file.path(work_dir, "snpEff_data")
  if(file.exists(snpeff_data))  unlink(snpeff_data, recursive = T) # delete the folder
  dir.create(snpeff_data) # create the data folder
  dir.create(file.path(snpeff_data, dset)) # create the folder with dset name

  if(!is.null(ref_path)) file.copy(ref_path, file.path(snpeff_data, dset, "sequences.fa")) # copy the reference fasta if provided (will this clash if the gbk/gff also has the sequence?)

  if(!is.null(gbk_path)) {
    file.copy(gbk_path, file.path(snpeff_data, dset, "genes.gbk")) # copy the gff annotation file if not null

    if(isWin){
      cat("COMMAND :=:=:", paste('java -jar',
                                 shQuote(normalizePath(snpeff_jar)), 'build -genbank -config',
                                 shQuote(normalizePath(snpeff_config)), '-dataDir', "snpeff_data",
                                 '-v', shQuote(dset)), '\n')

      shell(paste('java -jar',
                  shQuote(normalizePath(snpeff_jar)), 'build -genbank -config',
                  shQuote(normalizePath(snpeff_config)), '-dataDir', "snpeff_data",
                  '-v', shQuote(dset)), translate = T) # Build index for snpEff
    } else {
      cat("COMMAND :=:=:", paste('java -jar', shQuote(snpeff_jar),
                                 'build -genbank -config', shQuote(snpeff_config),
                                 '-dataDir', shQuote(snpeff_data),
                                 '-v', shQuote(dset)), '\n')

      system(paste('java -jar', shQuote(snpeff_jar),
                   'build -genbank -config', shQuote(snpeff_config),
                   '-dataDir', shQuote(snpeff_data),
                   '-v', shQuote(dset)
      )) # Build index for snpEff
    }

  }
  if(!is.null(gff_path)) {
    file.copy(gff_path, file.path(snpeff_data, dset, "genes.gff")) # copy the gff annotation file if not null
    if(isWin){
      cat("COMMAND :=:=:", paste('java -jar', shQuote(normalizePath(snpeff_jar)), 'build -gff3 -noCheckCds -noCheckProtein -config',
                                 shQuote(normalizePath(snpeff_config)), '-dataDir', "snpeff_data",
                                 '-v', shQuote(dset)), '\n')

      shell(paste('java -jar', shQuote(normalizePath(snpeff_jar)), 'build -gff3 -noCheckCds -noCheckProtein -config',
                  shQuote(normalizePath(snpeff_config)), '-dataDir', "snpeff_data",
                  '-v', shQuote(dset))) # Build index for snpEff
    } else {

      cat("COMMAND :=:=:", paste('java -jar', shQuote(snpeff_jar), 'build -gff3 -noCheckCds -noCheckProtein -config',
                                 shQuote(snpeff_config),
                                 '-dataDir', shQuote(snpeff_data),
                                 '-v', shQuote(dset)), '\n')

      system(paste('java -jar', shQuote(snpeff_jar), 'build -gff3 -noCheckCds -noCheckProtein -config',
                   shQuote(snpeff_config),
                   '-dataDir', shQuote(snpeff_data),
                   '-v', shQuote(dset))) # Build index for snpEff
    }

  }

  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
  return(1)
}

append_vcf_header = function(path, g){
  write("##fileformat=VCF4.1", file = path, append = F)
  write(paste("##contig=<ID=1,length=",g,">", sep = ""), file = path, append = T)
  write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", file = path, append = T)
  write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO", file = path, append = T)
}

create_vcf_file = function(genome_name, snps_to_ann, REF, ALT, vcf_write_path){
  vcf = data.frame(CHROM = genome_name,
                   POS = sprintf("%.0f", snps_to_ann),
                   ID = '.',
                   REF = REF,
                   ALT = ALT,
                   QUAL = '.',
                   FILTER = '.',
                   INFO = '.')
  write.table(x = vcf, file = vcf_write_path, append = T, row.names = F, col.names = F, quote = F, sep = '\t')
}

# annotate with snp_eff
run_snpeff = function(dset = NULL, genome_name = NULL, snpeff_jar = NULL, work_dir = NULL, vcf_write_path = NULL, vcf_annotated_path = NULL){
  #java -Xmx16G -jar snpEff.jar -v -stats Aus0004_efcm_link_annotated.html Enterococcus_faecium_gcf_000250945.1 ../bac_Efaecium/efcm_link_snps.vcf > Aus0004_efcm_link_annotated.vcf
  if(length(grep("windows", .Platform$OS.type))) isWin = T else isWin = F # We might need more hacks for other OSes
  snpeff_config = file.path(work_dir, "snpEff.config")
  snpeff_data = file.path(work_dir, "snpEff_data")

  if(isWin){
    cat("COMMAND :=:=:", paste('java -Xmx16G -jar', shQuote(normalizePath(snpeff_jar)), '-v -dataDir', "snpeff_data", '-config',
                               shQuote(normalizePath(snpeff_config)),
                               shQuote(dset),
                               shQuote(vcf_write_path), '>', shQuote(vcf_annotated_path)), '\n')


    shell(paste('java -Xmx16G -jar', shQuote(normalizePath(snpeff_jar)), '-v -dataDir', "snpeff_data", '-config',
                shQuote(normalizePath(snpeff_config)),
                shQuote(dset),
                shQuote(vcf_write_path), '>', shQuote(vcf_annotated_path)), translate = T)
  } else {
    cat("COMMAND :=:=:", paste('java -Xmx16G -jar', shQuote(snpeff_jar), '-v -dataDir',
                               shQuote(snpeff_data), '-config',
                               shQuote(snpeff_config),
                               shQuote(dset),
                               shQuote(vcf_write_path), '>', shQuote(vcf_annotated_path)), '\n')

    system(paste('java -Xmx16G -jar', shQuote(snpeff_jar), '-v -dataDir',
                 shQuote(snpeff_data), '-config',
                 shQuote(snpeff_config),
                 shQuote(dset),
                 shQuote(vcf_write_path), '>', shQuote(vcf_annotated_path)))
  }

  status = 1
  return(status)
}

convert_vcfann_to_table = function(vcf_annotated_path, annotations_path, snps_to_ann_idx, allele_table, nseq){
  ann = read.table(vcf_annotated_path, quote = "") # we can read the annotated VCF as a table into R
  ann = ann[,c(2,4,5,8)]
  colnames(ann) = c("pos", "REF", "ALT", "ANN")
  # There can be weird quotes in some annotations, remove them from the ANN field
  for(x in 1:nrow(ann)){
    ann$ANN[x] = gsub('[\"]', '',  ann$ANN[x])
  }

  ann = cbind(ann, annotation =  sapply(ann$ANN, function(x) unlist(strsplit(x, "\\|"))[2]))
  ann = cbind(ann, description =  sapply(ann$ANN, function(x) paste(unique(unlist(strsplit(x, "\\|"))[c(4,5,10,11)]), collapse = ":") ))
  ann = cbind(ann, cds =  sapply(ann$ANN,   function(x) unlist(strsplit(x, "\\|"))[5]))

  # There can still be some weird quotes in some annotations, remove them from each field
  for(x in 1:nrow(ann)){
    # Possibly redundant?
    ann$annotation[x] = gsub('[\"]', '',  ann$annotation[x])
    ann$description[x] = gsub('[\"]', '',  ann$description[x])
    ann$cds[x] = gsub('[\"]', '',  ann$cds[x])
  }


  code = rep("ns", nrow(ann)) # start as everything is ns
  code[grep("synonymous_variant", ann$annotation)] = "sy"
  code[grep("stop_retained_variant", ann$annotation)] = "sy"
  code[grep("downstream_gene_variant", ann$annotation)] = "ig"
  code[grep("upstream_gene_variant", ann$annotation)] = "ig"

  ann$code = code
  # table(ann$code)

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

add_annotations_to_links = function(links_red, ann, links_annotated_path, links_type = "SR"){

  # t0 = Sys.time()
  # a bit faster now
  pp_ = 0
  vidx = 0
  nlinks = nrow(links_red)
  pos1_idx = rep(NA, nlinks)
  pos2_idx = rep(NA, nlinks)

  pb = utils::txtProgressBar(min = 1, max = nlinks, initial = 1)

  for(idx in 1:nlinks){
    utils::setTxtProgressBar(pb,idx)
    pp = links_red$pos1[idx]
    if(pp != pp_){
      vidx = pos1_idx[idx] = which(ann$pos %in% pp)
      pp_ = pp
    } else {
      pos1_idx[idx] = vidx
    }

    pos2_idx[idx] = which(ann$pos %in% links_red$pos2[idx])
  }

  if(links_type == "SR"){
    l1_a1_d = data.frame(pos1 = links_red$pos1,
                         pos2 = links_red$pos2,
                         len = links_red$len,
                         ARACNE = links_red$ARACNE,
                         MI = links_red$MI,
                         srp = links_red$srp_max,
                         pos1_ann = ann$description[pos1_idx],
                         pos2_ann = ann$description[pos2_idx],
                         pos1_genreg = ann$cds[pos1_idx],
                         pos2_genreg = ann$cds[pos2_idx],
                         links = paste(ann$code[pos1_idx], ann$code[pos2_idx], sep = "X"),
                         pos1_ad = ann$allele_dist[pos1_idx],
                         pos2_ad = ann$allele_dist[pos2_idx])

    # saving annotated links
    # sort these links here
    l1_a1_d = l1_a1_d[order(l1_a1_d$srp, decreasing = T), ]
    rownames(l1_a1_d) = NULL
  } else if(links_type == "LR"){
    l1_a1_d = data.frame(pos1 = links_red$pos1,
                         pos2 = links_red$pos2,
                         len = links_red$len,
                         ARACNE = links_red$ARACNE,
                         MI = links_red$MI,
                         pos1_ann = ann$description[pos1_idx],
                         pos2_ann = ann$description[pos2_idx],
                         pos1_genreg = ann$cds[pos1_idx],
                         pos2_genreg = ann$cds[pos2_idx],
                         links = paste(ann$code[pos1_idx], ann$code[pos2_idx], sep = "X"),
                         pos1_ad = ann$allele_dist[pos1_idx],
                         pos2_ad = ann$allele_dist[pos2_idx])

    # saving annotated links
    # sort these links here
    l1_a1_d = l1_a1_d[order(l1_a1_d$MI, decreasing = T), ]
    rownames(l1_a1_d) = NULL

  }

  write.table(l1_a1_d, file = links_annotated_path, sep = '\t', col.names = T, row.names = F, quote = F)
  return(l1_a1_d)
}

detect_top_hits = function(l1_a1_d, max_tophits = 250, tophits_path){
  l1_a1_d = l1_a1_d[which(l1_a1_d$ARACNE == 1),] # ARACNE direct only
  l1_a1_d = l1_a1_d[which(l1_a1_d$links!='syXsy'),] # remove syXsy links
  l1_a1_d = l1_a1_d[which(l1_a1_d$pos1_genreg != l1_a1_d$pos2_genreg), ] # remove links from the same region

  if(nrow(l1_a1_d) > max_tophits) l1_a1_d = l1_a1_d[1:max_tophits, ] # all top hits or max_tophits
  write.table(l1_a1_d, file = tophits_path, sep = '\t', col.names = T, row.names = F, quote = F)

  return(l1_a1_d)

}
