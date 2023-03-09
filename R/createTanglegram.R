#' create_tanglegram
#'
#' Function to generate long and short range GWES plots
#'
#' @importFrom htmlwidgets saveWidget
#' @importFrom stats cutree hclust dist
#' @importFrom genbankr cds
#' @importFrom plyr . ddply
#' @importFrom chromoMap chromoMap
#'
#' @param tophits data frame with top short range GWES links, returned from LDWeaver::perform_snpEff_annotations()
#' @param gbk output from parsing the genbank file using LDWeaver::parse_genbank_file() (default = NULL), only provide either GFF3 or GBK annotation
#' @param gff output from parsing the gff3 file using LDWeaver::parse_gff_file() (default = NULL), only provide either GFF3 or GBK annotation
#' @param tanglegram_folder folder to save tanglegram(s)
#' @param break_segments specify the number of genome segments to prepare - one tanglegram per segment (default = 5)
#' @param links_type specify the links type long-range "LR" or short-range "SR" (default = "SR"). A tanglegram may not be suitable to visualise long range links!
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' gbk = LDWeaver::parse_genbank_file("path_to_genbank_file", length_check = F)
#' create_tanglegram(tophits, gbk, "save_folder")
#' }
#'
#' @export
create_tanglegram = function(tophits, gbk = NULL, gff = NULL, tanglegram_folder, break_segments = 5, links_type = "SR"){
  if( (is.null(gbk) & is.null(gff)) | (!is.null(gbk) & !is.null(gff)) ) stop("Provide either one of gbk or gff")

  # if(!is.null(gff)) return(0) # temporary hold for gff until the function is ready

  p1a = p2a = NULL # avoid plyr NSE issue (https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/)
  if(!file.exists(tanglegram_folder)) dir.create(tanglegram_folder)
  cat("Preparing Tanglegrams ... ")
  t0 = Sys.time()
  tophits$dummy_chrom = cutree(hclust(dist(tophits$pos1)), break_segments); #print(table(tophits$dummy_chrom))

  # we should swap labels depending on min
  clst_brk_ord = order(sapply(1:break_segments, function(x) min(tophits$pos1[tophits$dummy_chrom == x])))
  dc_tmp = tophits$dummy_chrom
  change = F
  for(i in 1:break_segments){
    if(i == clst_brk_ord[i]) {
      next
    } else{ # swap needed
      dc_tmp[tophits$dummy_chrom == i] = clst_brk_ord[i]
      change = T
    }
  }
  if(change) tophits$dummy_chrom = dc_tmp

  # gbkv = as.data.frame(genbankr::cds(gbk)) # extract regions from gbk

  chr_view = list()
  chr_file = list()
  ann_file = list()
  link_data = list()


  for(clustidx in 1:break_segments){
    tophits_dumchr = tophits[which(tophits$dummy_chrom == clustidx), ]
    if(links_type == "SR"){
      df = data.frame(p1a = tophits_dumchr$pos1_genreg,
                      p2a = tophits_dumchr$pos2_genreg,
                      w = tophits_dumchr$srp)
    } else if(links_type == "LR"){
      df = data.frame(p1a = tophits_dumchr$pos1_genreg,
                      p2a = tophits_dumchr$pos2_genreg,
                      w = tophits_dumchr$MI)
    } else {
      stop("Links type must be SR or LR")
    }


    df_uq =  plyr::ddply(df, plyr::.(p1a,p2a), nrow)
    # get the max srp for each link
    for(i in 1:nrow(df_uq)){
      df_uq$w[i] = max(df$w[which(df_uq$p1a[i] == df$p1a & df_uq$p2a[i] == df$p2a)])
      df_uq$p1a[i] = unique(df$p1a[which(df_uq$p1a[i] == df$p1a)])
      df_uq$p2a[i] = unique(df$p2a[which(df_uq$p2a[i] == df$p2a)])
    }

    all_locs = unique(c(df_uq$p1a, df_uq$p2a))
    # loc_strt_end = t(unname(sapply(all_locs, function(x) {idx = grep(x, gbkv$locus_tag); return(c(unique(gbkv$start[idx])[1], unique(gbkv$end[idx])[1]))})))
    loc_strt_end = matrix(rep(NA, length(all_locs)*2), ncol = 2)

    all_locs_notfound = c()
    #### GBK VERSION
    if(!is.null(gbk)){
      for(loc_idx in 1:length(all_locs)){
        idx = c()
        # Let's search through the gbk file for the locus
        if(length(gbk@genes@elementMetadata$locus_tag)){
          idx = grep(all_locs[loc_idx], gbk@genes@elementMetadata@listData$locus_tag)
          if(length(idx) > 0) {
            loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
            loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
            next
          }
        }
        if(length(gbk@cds@elementMetadata$locus_tag)){
          idx = grep(all_locs[loc_idx], gbk@cds@elementMetadata@listData$locus_tag)
          if(length(idx) > 0) {
            loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
            loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
            next
          }
        }
        if(length(gbk@exons@elementMetadata$locus_tag)){
          idx = grep(all_locs[loc_idx], gbk@exons@elementMetadata@listData$locus_tag)
          if(length(idx) > 0) {
            loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
            loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
            next
          }
        }
        if(length(gbk@transcripts@elementMetadata$locus_tag)){
          idx = grep(all_locs[loc_idx], gbk@transcripts@elementMetadata@listData$locus_tag)
          if(length(idx) > 0) {
            loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
            loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
            next
          }
        }
        if(length(gbk@other_features@elementMetadata$locus_tag)){
          idx = grep(all_locs[loc_idx], gbk@other_features@elementMetadata@listData$locus_tag)
          if(length(idx) > 0) {
            loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
            loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
            next
          }
        }

        if(length(idx) == 0) {
          warning(paste("Could not locate", all_locs[loc_idx], 'in the genbankr parsed gbk file, these link will be dropped from the tanglegram...'))
          all_locs_notfound = c(all_locs_notfound, all_locs[loc_idx])
        }
      }

    }
    #### GFF VERSION
    if(!is.null(gff)){
      for(loc_idx in 1:length(all_locs)){
        idx = c()
        # Let's search through the gbk file for the locus
        # Note: For the GFF based annotations, missing genes are replaced with a tag starting with <GENE_>, sub it out

        idx = grep(gsub('GENE_', '', all_locs[loc_idx]), gff$gff$attributes)
        if(length(idx) > 0) {
          loc_strt_end[loc_idx, 1] = unique(gff$gff$start[idx])[1]
          loc_strt_end[loc_idx, 2] = unique(gff$gff$end[idx])[1]
          next
        }
      }

      if(length(idx) == 0) {
        warning(paste("Could not locate", all_locs[loc_idx], 'in the genbankr parsed gbk file, these link will be dropped from the tanglegram...'))
        all_locs_notfound = c(all_locs_notfound, all_locs[loc_idx])
      }
    }



    # prune df_uq if any locs are missing
    if(length(all_locs_notfound) > 0){
      prune_rows = c()
      for(rmidx in 1:length(all_locs_notfound)){
        prune_rows = c(prune_rows, c(grep(all_locs_notfound[rmidx], df_uq$p1a), grep(all_locs_notfound[rmidx], df_uq$p2a)))
      }
      prune_rows = sort(unique(prune_rows))
      df_uq = df_uq[-prune_rows, ]

      # prune all_locs
      # This is a quick hack, redo the loc_strt_end if any items end up missing in the data
      all_locs = unique(c(df_uq$p1a, df_uq$p2a))

      loc_strt_end = matrix(rep(NA, length(all_locs)*2), ncol = 2)

      #### GBK VERSION
      if(!is.null(gbk)){
        for(loc_idx in 1:length(all_locs)){
          idx = c()
          # Let's search through the gbk file for the locus
          if(length(gbk@genes@elementMetadata$locus_tag)){
            idx = grep(all_locs[loc_idx], gbk@genes@elementMetadata@listData$locus_tag)
            if(length(idx) > 0) {
              loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
              loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
              next
            }
          }
          if(length(gbk@cds@elementMetadata$locus_tag)){
            idx = grep(all_locs[loc_idx], gbk@cds@elementMetadata@listData$locus_tag)
            if(length(idx) > 0) {
              loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
              loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
              next
            }
          }
          if(length(gbk@exons@elementMetadata$locus_tag)){
            idx = grep(all_locs[loc_idx], gbk@exons@elementMetadata@listData$locus_tag)
            if(length(idx) > 0) {
              loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
              loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
              next
            }
          }
          if(length(gbk@transcripts@elementMetadata$locus_tag)){
            idx = grep(all_locs[loc_idx], gbk@transcripts@elementMetadata@listData$locus_tag)
            if(length(idx) > 0) {
              loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
              loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
              next
            }
          }
          if(length(gbk@other_features@elementMetadata$locus_tag)){
            idx = grep(all_locs[loc_idx], gbk@other_features@elementMetadata@listData$locus_tag)
            if(length(idx) > 0) {
              loc_strt_end[loc_idx, 1] = unique(gbk@genes@ranges@start[idx])[1]
              loc_strt_end[loc_idx, 2] = unique(gbk@genes@ranges@start[idx] + gbk@genes@ranges@width[idx] - 1)[1]
              next
            }
          }
        }
      }

      #### GFF VERSION
      if(!is.null(gff)){
        for(loc_idx in 1:length(all_locs)){
          idx = c()
          # Let's search through the gbk file for the locus
          # Note: For the GFF based annotations, missing genes are replaced with a tag starting with <GENE_>, sub it out

          if(nrow(gff$gff) > 0){
            idx = grep(gsub('GENE_', '', all_locs[loc_idx]), gff$gff$attributes)
            if(length(idx) > 0) {
              loc_strt_end[loc_idx, 1] = unique(gff$gff$start[idx])[1]
              loc_strt_end[loc_idx, 2] = unique(gff$gff$end[idx])[1]
              next
            }
          }
        }
      }

    }

    # let's map back to annotations
    all_ans = c()
    for(i in 1:length(all_locs)){
      pm1 = c()
      pm2 = c()
      m1 = which(df_uq$p1a %in% all_locs[i])
      m2 = which(df_uq$p2a %in% all_locs[i])
      if(length(m1)>0) pm1 = unique(df_uq$p1a[m1])
      if(length(m2)>0) pm2 = unique(df_uq$p2a[m2])
      pm = unique(c(pm1, pm2))

      if(length(pm) == 1){
        all_ans = c(all_ans, pm)
      } else {
        cat(paste("Check:", i, "\n"))
      }
    }

    chr_file[[clustidx]] = data.frame(V1 = c("p", "q"),
                                      V2 = min(loc_strt_end[,1])-1e3,
                                      V3 = max(loc_strt_end[,2])+1e3)


    ann_file[[clustidx]] = data.frame(V1 = c(paste("p_", all_ans, sep = ""), paste("q_", all_ans, sep = "")),
                                      V2 = c(rep("p", length(all_locs)),  rep("q", length(all_locs))),
                                      V3 = loc_strt_end[,1],
                                      V4 = loc_strt_end[,2])


    link_data[[clustidx]] = data.frame(V1 = paste("p_", df_uq$p1a, sep = ""), V2 = 1,
                                       V3 = paste("q_", df_uq$p2a,sep = ""), V4 = 1)

    chr_view[[clustidx]] = chromoMap::chromoMap(ch.files = list(chr_file[[clustidx]]),
                                                data.files = list(ann_file[[clustidx]]),
                                                show.links = T,
                                                loci_links = link_data[[clustidx]],
                                                links.colors = "red2",
                                                n_win.factor = 3,
                                                labels = T,
                                                ch_gap = 50,
                                                y_chr_scale = 50,
                                                top_margin = 100,
                                                chr_length = 6,
                                                chr_width = 25)

    chr_plot_path = file.path(tanglegram_folder, paste("tng_", clustidx, ".html", sep = ""))
    htmlwidgets::saveWidget( widget = chr_view[[clustidx]], file = chr_plot_path, selfcontained = T)

  }
  cat(paste("Done in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s"))
}
