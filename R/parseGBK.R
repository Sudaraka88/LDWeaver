#' parse_genbank_file
#'
#' Function to parse the genbank file for the fasta alignment. Parser built using genbankr source available here: https://github.com/gmbecker/genbankr.
#' @importFrom methods is new
#' @importFrom GenomeInfoDb seqnames seqlevels seqinfo Seqinfo
#' @importFrom Biostrings DNAString AAString extractAt AAStringSet
#' @importFrom IRanges IRanges ranges heads
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @importFrom VariantAnnotation VRanges makeVRangesFromGRanges
#' @importFrom utils stack
#' @import GenomicRanges
#' @param gbk_path path to genbank file
#' @param g sequence length, available from the LDWeaver::parse_fasta_alignment() output (default = NULL),
#' required if <length_check = T>
#' @param length_check specify whether to check if fasta and gbk sequence lengths are equal (default = T)
#'
#' @return GenBankRecord object
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "sample.aln", package = "LDWeaver")
#' snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
#' gbk_path <- system.file("extdata", "sample.gbk", package = "LDWeaver")
#' gbk <- parse_genbank_file(gbk_path, snp.dat$g)
#' }
#' @export
parse_genbank_file = function(gbk_path, g = NULL, length_check = T){

  # TODO: Fix later
  #### I can't get this function to work by only importing the required parts from GenomicRanges:
  # #'importFrom GenomicRanges GRanges mcols findOverlaps GRangesList resize
  # #'importMethodsFrom GenomicRanges  range
  # I'm importing the whole package instead
  # Fails at txslst = range(spl, ignore.strand  = F), it seems the S4 method range() does not get exported

  # TODO: This parser is almost identical to genbankr and can probably be improved for speed and reliability.
  # More imp. we only need to extract CDS start/stop positions and the fasta reference sequence for LDWeaver: should be simplified eventually...
  ##########################################

  t0 = Sys.time()
  # Check inputs
  if(length_check){ # perform the length check
    if(is.null(g)){
      # then g cannot be null
      stop("g must be provided to perform length check!\n")
      # return(-1)
    }
  }

  if(!file.exists(gbk_path)) stop(paste("Can't locate file", gbk_path))

  ### INSERT CODE HERE
  gbk = suppressWarnings(readGenBank2(gbk_path))
  refseq = gbk@sequence
  # refseq = genbankr::getSeq(gbk)


  if(length(refseq) != 1){
    stop("The GBK file should contain the reference sequence!\n")
    # return(-1)
  }
  # the length check is good to check if the alignment matches with the gbk, setting it to F will stop it
  ref_g = length(refseq[[1]]) # length of the reference sequence

  if(length_check){ # perform the length check
    if(ref_g != g){ # perform check
      stop("Genbank reference sequence length mismatches with the fasta alignment!\n")
      # return(-1)
    }

  } else {
    if(!is.null(g)){
      if(ref_g != g){
        warning("Fasta length does not match the genbank reference sequence length!\n")
      }
    } else {
      warning("Similarity between the genbank reference and fasta sequences NOT checked, ignore if <pos> was provided...\n")
    }
  }
  cat(paste("Successfully read gbk file:", gbk_path, "in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))

  # genbankr::seqinfo(gbk)
  return(list(gbk = gbk,
              ref_g = ref_g))

}

### readGenBank from genbankr
readGenBank2 = function(file, text = readLines(file), partial = NA) {

  verbose = FALSE
  # I am skipping this because we test for file validity separately
  # if(missing(text)) {
  #   if(missing(file))
  #     stop("One of text or file must be specified.")
  #   if(is(file, "GBAccession"))
  #     text = .getGBfromNuccore(file)
  #   else if (!file.exists(file))
  #     stop("file does not appear to an existing file or a GBAccession object. Valid file or text argument is required.")
  # }

  record_brk_re = "^[[:space:]]*//[[:space:]]*$" ## the spaces probably aren't needed here...


  if(methods::is(text, "character") &&
     length(recbrks <- grep(record_brk_re, text)) > 1) {
    sq = seq_along(text)
    indfac = cut(sq, breaks = recbrks[-length(recbrks)])
    text = split(text, indfac)
    stopifnot(all(sapply(text,
                         function(x) identical(substr(x[1], 1, 5), "LOCUS"))))
  }
  # if(is(text, "list"))
  #   return(lapply(text, function(txt) readGenBank(text = txt, partial = partial, ret.seq= ret.seq, verbose = verbose)))
  ## we always read in sequence because it is required for variant annotations
  ## we throw it away after if the user set ret.seq=FALSE
  prsed = parseGenBank(text = text, partial = partial, verbose = verbose,
                       ret.seq = TRUE)
  ret = make_gbrecord(rawgbk = prsed, verbose = verbose)
  # if(!ret.seq) sequence(ret) = NULL
  ret
}


parseGenBank = function(file, text = readLines(file),  partial = NA,
                        verbose = FALSE,
                        ret.anno = TRUE,
                        ret.seq = TRUE) {

  prime_field_re = "^[[:upper:]]" #+[[:space:]]+" # bakta ORIGIN has no space

  # remove fuzzy start/end points
  text <- gsub("[<>]", "", text)

  if(!ret.anno && !ret.seq)
    stop("Must return at least one of annotations or sequence.")
  bf = proc.time()["elapsed"]
  if(missing(text) && !file.exists(file))
    stop("No text provided and file does not exist or was not specified. Either an existing file or text to parse must be provided.")
  if(length(text) == 1)
    text = fastwriteread(text)

  fldlines = grepl(prime_field_re, text)
  fldfac = cumsum(fldlines)
  fldnames = gsub("^([[:upper:]]+).*", "\\1", text[fldlines])[fldfac]

  spl = split(text, fldnames)

  resthang = list(LOCUS = readLocus(spl[["LOCUS"]]))
  resthang[["FEATURES"]] = readFeatures(spl[["FEATURES"]],
                                        source.only=!ret.anno,
                                        partial = partial)
  seqtype = .seqTypeFromLocus(resthang$LOCUS)

  if(ret.seq) resthang$ORIGIN = readOrigin(spl[["ORIGIN"]], seqtype = seqtype)

  if(ret.anno) {
    resthang2 = mapply(function(field, lines, verbose) {
      switch(field,
             DEFINITION = readDefinition(lines),
             ACCESSION = readAccession(lines),
             VERSION = readVersions(lines),
             KEYWORDS = readKeywords(lines),
             SOURCE = readSource(lines),
             ## don't read FEATURES, ORIGIN, or LOCUS because they are
             ## already in resthang from above
             NULL)
    }, lines = spl, field = names(spl), SIMPLIFY=FALSE, verbose = verbose)
    resthang2$FEATURES = resthang2$FEATURES[sapply(resthang2$FEATURES,
                                                   function(x) length(x)>0)]
    resthang2 = resthang2[!names(resthang2) %in% names(resthang)]
    resthang = c(resthang, resthang2)
  }
  ##DNAString to DNAStringSet
  origin = resthang$ORIGIN
  if(ret.seq && length(origin) > 0) {
    typs = sapply(resthang$FEATURES, function(x) x$type[1])
    srcs = fill_stack_df(resthang$FEATURES[typs == "source"])
    ## dss = DNAStringSet(lapply(GRanges(ranges(srcs), function(x) origin[x])))
    dss = Biostrings::extractAt(origin, IRanges::ranges(srcs))
    names(dss) = as.character(GenomeInfoDb::seqnames(srcs))
    if(!ret.anno)
      resthang = dss
    else
      resthang$ORIGIN = dss
  } else if (!ret.anno) { ##implies ret.seq is TRUE
    stop("Asked for only sequence (ret.anno=FALSE) from a file with no sequence information")
  }
  af = proc.time()["elapsed"]
  if(verbose)
    message("Done Parsing raw GenBank file text. [ ", af-bf, " seconds ]")
  resthang

}

make_gbrecord = function(rawgbk, verbose = FALSE) {
  bf = proc.time()["elapsed"]
  feats = rawgbk$FEATURES
  sq = rawgbk$ORIGIN

  typs = sapply(feats, function(x)
    if(length(x) > 0) x$type[1] else NA_character_)
  empty = is.na(typs)
  feats = feats[!empty]
  typs = typs[!empty]
  featspl = split(feats, typs)
  srcs = fill_stack_df(featspl$source)
  circ = rep(grepl("circular", rawgbk$LOCUS), times = length(srcs))
  ##grab the versioned accession to use as the "genome" in seqinfo
  genom = gn_from_vers(rawgbk$VERSION)
  sqinfo = GenomeInfoDb::Seqinfo(GenomeInfoDb::seqlevels(srcs), width(srcs), circ, genom)
  if(verbose)
    message(Sys.time(), " Starting creation of gene GRanges")
  gns = make_genegr(featspl$gene, sqinfo)

  if(verbose)
    message(Sys.time(), " Starting creation of CDS GRanges")
  if(!is.null(featspl$CDS))
    cdss = make_cdsgr(featspl$CDS, gns, sqinfo)
  else
    cdss = GenomicRanges::GRanges(seqinfo = sqinfo)
  if(verbose)
    message(Sys.time(), " Starting creation of exon GRanges")


  exns = make_exongr(featspl$exon, cdss = cdss, sqinfo)

  if(verbose)
    message(Sys.time(), " Starting creation of variant VRanges")
  vars = make_varvr(featspl$variation, sq = sq, sqinfo)


  if(verbose)
    message(Sys.time(), " Starting creation of transcript GRanges")

  txs = make_txgr(featspl$mRNA, exons = exns, sqinfo, genes = gns)

  if(verbose)
    message(Sys.time(), " Starting creation of misc feature GRanges")

  ofeats = fill_stack_df(feats[!typs %in% c("gene", "exon", "CDS", "variation",
                                            "mRNA", "source")])

  ofeats$temp_grouping_id = NULL

  if(is.null(ofeats))
    ofeats = GenomicRanges::GRanges()
  GenomeInfoDb::seqinfo(ofeats) = sqinfo
  res = methods::new("GenBankRecord", genes = gns, cds = cdss, exons = exns,
                     transcripts = txs, variations = vars,
                     sources = fill_stack_df(feats[typs == "source"]),
                     other_features = ofeats,
                     accession = rawgbk$ACCESSION %||% NA_character_,
                     version = rawgbk$VERSION %||% NA_character_,
                     locus = rawgbk$LOCUS,
                     definition = rawgbk$DEFINITION,
                     sequence = sq)
  af = proc.time()["elapsed"]
  if(verbose)
    message(Sys.time(), " - Done creating GenBankRecord object [ ", af - bf, " seconds ]")
  res
}

## slightly specialized function to stack granges which may have different
## mcols together, filling with NA_character_ as needed.
## also collapses multiple db_xref notes into single CharacterList column and
## creates an AAStringSet for the translation field
fill_stack_df = function(dflist, cols, fill.logical = TRUE, sqinfo = NULL) {
  if(length(dflist) == 0 )
    return(NULL)

  multivalfields = c("db_xref", "EC_number", "gene_synonym", "old_locus_tag")

  dflist = dflist[sapply(dflist, function(x) !is.null(x) && nrow(x) > 0)]


  allcols = unique(unlist(lapply(dflist, function(x) names(x))))
  basenms = gsub("(.*)(\\.[[:digit:]]+)$", "\\1", allcols)
  nmtab = table(basenms)
  dupnms = names(nmtab[nmtab>1])
  if(any(!dupnms %in% multivalfields))
    warning("Got unexpected multi-value field(s) [ ",
            paste(setdiff(dupnms, multivalfields), collapse = ", "),
            " ]. The resulting column(s) will be of class CharacterList, rather than vector(s). Please contact the maintainer if multi-valuedness is expected/meaningful for the listed field(s).")
  allcols = unique(basenms)

  logcols = unique(unlist(lapply(dflist, function(x) names(x)[sapply(x, is.logical)])))
  charcols = setdiff(allcols, logcols)

  if(missing(cols))
    cols = allcols

  filled = mapply(
    function(x, i) {
      ## have to deal with arbitrary multiple columns
      ## transform them into list columns
      for(nm in dupnms) {
        locs = grep(nm, names(x))
        if(length(locs)) {
          rows = lapply(seq(along = rownames(x)),
                        function(y) unlist(x[y,locs]))

          x = x[,-locs]
        } else {
          rows = list(character())
        }

        x[[nm]] = rows
      }



      ## setdiff is not symmetric
      missnm = setdiff(charcols, names(x))
      x[,missnm] = NA_character_
      falsenm = setdiff(logcols, names(x))
      x[,falsenm] = FALSE
      x = x[,cols]
      x$temp_grouping_id = i
      x
    }, x = dflist, i = seq(along = dflist), SIMPLIFY=FALSE)
  stk = .simple_rbind_dataframe(filled, "temp")
  stk[["temp"]] = NULL

  listcols = which(sapply(names(stk),
                          function(x) methods::is(stk[[x]], "list") ||
                            x %in% multivalfields))
  stk[listcols] = lapply(listcols, function(i) as(stk[[i]], "CharacterList"))
  mc = names(stk)[!names(stk) %in% c("seqnames", "start", "end", "strand")]
  if(fill.logical) {
    logcols = which(sapply(stk, is.logical))
    stk[,logcols] = lapply(logcols, function(i) {
      dat = stk[[i]]
      dat[is.na(dat)] = FALSE
      dat
    })
  }
  grstk = GenomicRanges::GRanges(seqnames = stk$seqnames,
                                 ranges = IRanges::IRanges(start = stk$start, end = stk$end),
                                 strand = stk$strand )

  ## this may be slightly slower, but specifying mcols during
  ## creation appends mcols. to all the column names, super annoying.
  GenomicRanges::mcols(grstk) = stk[,mc]
  if("translation" %in% names(GenomicRanges::mcols(grstk))) {
    if(anyNA(grstk$translation)) {
      message("Translation product seems to be missing for ",
              sum(is.na(grstk$translation)),
              " of ", length(grstk), " ",
              grstk$type[1], " annotations. Setting to ''")
      grstk$translation[is.na(grstk$translation)] = ""
    }
    grstk$translation = Biostrings::AAStringSet(grstk$translation)
  }

  if(!is.null(sqinfo))
    GenomeInfoDb::seqinfo(grstk) = sqinfo
  grstk
}

# sec_field_re = "^( {5}|\\t)[[:alnum:]'_-]+[[:space:]]+(complement|join|order|[[:digit:]<,])"

## Functions to read in the fields of a GenBank file
## utility borrowed from Hadley et al. rlang
`%||%` <- function(a, b) if(is.null(a)) b else a

strip_fieldname = function(lines) gsub("^[[:space:]]*[[:upper:]]*[[:space:]]+", "", lines)

readLocus = function(line) {
  ## missing strip fieldname?
  spl = strsplit(line, "[\t]+")[[1]]
  spl
}

readDefinition = function(lines) {
  paste(strip_fieldname(lines), collapse = " ")
}

readAccession = function(line)  strip_fieldname(line)

readVersions = function(line) {
  txt = strip_fieldname(line)
  spl = strsplit(txt, "[[:space:]]+")[[1]]
  c(accession.version = spl[1], GenInfoID = gsub("GI:", "", spl[2]))
}

readKeywords = function(lines) {
  txt = strip_fieldname(lines)
  txt = paste(txt, collapse = " ")
  if(identical(txt, "."))
    txt = NA_character_
  txt
}

readSource = function(lines) {
  sec_field_re = "^( {5}|\\t)[[:alnum:]'_-]+[[:space:]]+(complement|join|order|[[:digit:]<,])"
  secfieldinds = grep(sec_field_re, lines)
  src = strip_fieldname(lines[1])
  lines = lines[-1] #consume SOURCE line
  org = strip_fieldname(lines[1])
  lines = lines[-1] # consume ORGANISM line
  lineage = strsplit(paste(lines, collapse = " "), split  = "; ")[[1]]
  list(source= src, organism = org, lineage = lineage)
}


chk_outer_complement = function(str) grepl("[^(]*complement\\(", str)
strip_outer_operator = function(str, op = "complement") {
  regex = paste0("[^(]*", op, "\\((.*)\\)")

  gsub(regex,"\\1", str)
}

.do_join_silliness = function(str, chr, ats, partial = NA, strand = NA) {

  #    if(grepl("^join", str))
  if(is.na(partial) || !partial) {
    part = grepl("[<>]", str)
    if(part) {
      if(is.na(partial))
        message("Partial range detected in a compound range (join or order) feature. Excluding entire annotation")
      return(
        data.frame(seqnames = character(),
                   start = numeric(),
                   end = numeric(),
                   strand = character(),
                   type = character(),
                   stringsAsFactors=FALSE)
      )
    }
  }

  sstr = substr(str, 1, 1)
  if(sstr == "j") ##join
    str = strip_outer_operator(str, "join")
  else if(sstr == "o") ## order
    str = strip_outer_operator(str, "order")
  spl = strsplit(str, ",", fixed=TRUE)[[1]]
  grs = lapply(spl, make_feat_gr, chr = chr, ats = ats,
               partial = partial, strand = strand)
  .simple_rbind_dataframe(grs)

}


make_feat_gr = function(str, chr, ats, partial = NA, strand = NA) {

  if(is.na(strand) && chk_outer_complement(str)) {
    strand = "-"
    str = strip_outer_operator(str)
  }


  sbstr = substr(str, 1, 4)
  if(sbstr == "join" || sbstr == "orde")
    return(.do_join_silliness(str = str, chr = chr,
                              ats = ats, partial = partial,
                              strand = strand))
  haslt = grepl("<", str, fixed=TRUE) || grepl(">", str, fixed=TRUE)


  ## control with option. Default to exclude those ranges entirely.
  ##    if( (haslt || hasgt ) && (is.na(partial) || !partial) ){
  if( haslt  && (is.na(partial) || !partial) ){
    if(is.na(partial))
      warning("Incomplete feature annotation detected. ",
              "Omitting feature at ", str)
    return(GenomicRanges::GRanges(character(), IRanges::IRanges(numeric(), numeric())))
  }

  ## format is 123, 123..789, or 123^124
  spl = strsplit(str, "[.^]{1,2}")[[1]]
  start = as.integer(gsub("<*([[:digit:]]+).*", "\\1", spl[1]))
  if(length(spl) == 1)
    end = start
  else
    end = as.integer(gsub(">*([[:digit:]]+).*", "\\1", spl[2]))
  if(grepl("^", str, fixed=TRUE )) {
    end = end - 1
    ats$loctype = "insert"
  } else
    ats$loctype = "normal"

  if(is.na(strand))
    strand = "+"

  gr = cheap_unsafe_data.frame(lst = c(seqnames = chr, start = start, end = end, strand = strand, ats))
  gr
}

read_feat_attr = function(line) {
  num =  grepl('=[[:digit:]]+(\\.[[:digit:]]+){0,1}$', line)
  val = gsub('[[:space:]]*/[^=]+($|="{0,1}([^"]*)"{0,1})', "\\2", line)
  mapply(function(val, num) {
    if(nchar(val)==0)
      TRUE
    else if(num)
      as.numeric(val)
    else
      val
  }, val = val, num = num, SIMPLIFY=FALSE)
}

### These two functions are just for API compatibility, replace with actual values later
getCDS = function(gbk){
  return(gbk@cds)
}

# getSeq = function(gbk){ # replaced
#   return(gbk@sequence)
# }

## XXX is the leading * safe here? whitespace is getting clobbered by line
##combination below I think it's ok
strip_feat_type = function(ln) gsub("^[[:space:]]*[[:alnum:]_'-]+[[:space:]]+((complement\\(|join\\(|order\\(|[[:digit:]<]+).*)", "\\1", ln)

readFeatures = function(lines, partial = NA, verbose = FALSE,
                        source.only = FALSE) {
  sec_field_re = "^( {5}|\\t)[[:alnum:]'_-]+[[:space:]]+(complement|join|order|[[:digit:]<,])"

  if(substr(lines[1], 1, 8) == "FEATURES")
    lines = lines[-1] ## consume FEATURES line
  fttypelins = grepl(sec_field_re, lines)
  featfactor  = cumsum(fttypelins)

  if(source.only) {
    srcfeats = which(substr(lines[fttypelins], 6, 11) == "source")
    keepinds = featfactor %in% srcfeats
    lines = lines[keepinds]
    featfactor = featfactor[keepinds]
  }

  ##scope bullshittery
  chr = "unk"

  totsources = length(grep("[[:space:]]+source[[:space:]]+[<[:digit:]]", lines[which(fttypelins)]))
  numsources = 0
  everhadchr = FALSE

  do_readfeat = function(lines, partial = NA) {

    ## before collapse so the leading space is still there
    type = gsub("[[:space:]]+([[:alnum:]_'-]+).*", "\\1", lines[1])
    ##feature/location can go across multpiple lines x.x why genbank? whyyyy
    attrstrts = cumsum(grepl("^[[:space:]]+/[^[:space:]]+($|=([[:digit:]]|\"))", lines))
    lines = tapply(lines, attrstrts, function(x) {
      paste(gsub("^[[:space:]]+", "", x), collapse="")
    }, simplify=TRUE)

    rawlocstring = lines[1]

    rngstr = strip_feat_type(rawlocstring)

    ## consume primary feature line
    lines = lines[-1]
    if(length(lines)) {
      attrs = read_feat_attr(lines)

      names(attrs) = gsub("^[[:space:]]*/([^=]+)($|=[^[:space:]].*$)", "\\1", lines)
      if(type == "source") {
        numsources <<- numsources + 1
        if("chromosome" %in% names(attrs)) {
          if(numsources > 1 && !everhadchr)
            stop("This file appears to have some source features which specify chromosome names and others that do not. This is not currently supported. Please contact the maintainer if you need this feature.")
          everhadchr <<- TRUE
          chr <<- attrs$chromosome
        } else if(everhadchr) {
          stop("This file appears to have some source features which specify chromosome names and others that do not. This is not currently supported. Please contact the maintainer if you need this feature.")
          ## this assumes that if one source has strain, they all will.
          ## Good assumption?
        } else if("strain" %in% names(attrs)) {
          chr <<- if(totsources == 1) attrs$strain else paste(attrs$strain, numsources, sep=":")
        } else {
          chr <<- if(totsources == 1) attrs$organism else paste(attrs$organism, numsources, sep=":")
        }
      }
    } else {
      attrs = list()
    }
    make_feat_gr(str = rngstr, chr = chr, ats = c(type = type, attrs),
                 partial = partial)

  }

  if(verbose)
    message(Sys.time(), " Starting feature parsing")
  resgrs = tapply(lines, featfactor, do_readfeat,
                  simplify=FALSE, partial = partial)
  if(verbose)
    message(Sys.time(), " Done feature parsing")

  resgrs

}


## seqtype: bp = base pairs (DNA/RNA), aa = amino acids (peptides/protein)

readOrigin = function(lines, seqtype = "bp") {
  ## strip spacing and line numbering
  regex = "([[:space:]]+|[[:digit:]]+|//)"

  dnachar = gsub(regex, "", lines[-1])
  chars = paste(dnachar, collapse="")
  if(any(nzchar(dnachar))) {
    switch(seqtype,
           bp = Biostrings::DNAString(chars),
           aa = Biostrings::AAString(chars),
           stop("Unknown origin sequence type: ", seqtype))
  } else
    NULL
}


fastwriteread = function(txtline) {
  f = file()
  on.exit(close(f))
  writeLines(txtline, con = f)
  readLines(f)

}



##LOCUS       ABD64816                1353 aa            linear   INV 12-MAR-2006
##LOCUS       HG941718.1           5109767 bp    DNA     linear       25-MAR-2025
.seqTypeFromLocus = function(locus) {
  gsub("^[[:space:]]*LOCUS[[:space:]]+[^[:space:]]+[[:space:]]+[[:digit:]]+[[:space:]]+([^[:space:]]+).*",
       "\\1",
       locus)

  # sub("^\\s*LOCUS\\s+\\S+\\s+\\d+\\s+bp\\s+(\\S+).*", "\\1", locus)

}



allna = function(vec) all(is.na(vec))

.getGeneIdLab = function(ranges, verbose, stoponfail = FALSE) {
  cols = names(GenomicRanges::mcols(ranges))
  if("gene_id" %in% cols && !allna(ranges$gene_id))
    "gene_id"
  else if("locus_tag" %in% cols)
    "locus_tag"
  else if ("gene" %in% cols) {
    if(verbose)
      message("Annotations don't have 'locus_tag' label, using 'gene' as gene_id column")
    "gene"
  } else if (stoponfail)
    stop("Unable to determine gene id column on feature GRanges object")
  else
    NULL
}

.getGeneIdVec = function(ranges) {
  res = .getGeneIdLab(ranges, verbose = TRUE)
  if(!is.null(res))
    res = GenomicRanges::mcols(ranges)[[res]]
  res
}

match_cds_genes = function(cds, genes) {
  ## XXX do we want "equal" here or within?
  hits = GenomicRanges::findOverlaps(cds, genes, type = "equal")
  cds$gene_id= NA_character_
  cds$gene_id[S4Vectors::queryHits(hits)] = genes$gene_id
  if(anyNA(cds$gene_id)) {
    warning("unable to determine gene for some CDS annotations.",
            " Ignoring these ", sum(is.na(cds$gene_id)), " segments")
    cds = cds[!is.na(cds$gene_id)]
  }
  cds
}

make_cdsgr = function(rawcdss, gns, sqinfo) {

  ##    rawcdss = sanitize_feats(rawcdss, cds_cols)
  rawcdss = fill_stack_df(rawcdss, sqinfo = sqinfo)
  ## double order gives us something that can be merged directly into what
  ## out of tapply

  havegenelabs = FALSE

  rawcdss$gene_id = .getGeneIdVec(rawcdss)

  if(!is.null(rawcdss$gene_id) && !anyNA(rawcdss$gene_id)) {

    o = order(order(rawcdss$gene_id))
    var = "gene_id"
  } else {
    message("genes not available for all CDS ranges, using internal grouping ids")
    var = "temp_grouping_id"
    o = order(order(rawcdss$temp_grouping_id))
  }
  idnum = unlist(tapply(rawcdss$temp_grouping_id, GenomicRanges::mcols(rawcdss)[[var]], function(x) as.numeric(factor(x)), simplify=FALSE))[o]
  newid = paste(GenomicRanges::mcols(rawcdss)[[var]], idnum, sep=".")
  if(var == "temp_grouping_id")
    newid = paste0(ifelse(is.na(GenomicRanges::mcols(rawcdss)$gene), "unknown_gene_",GenomicRanges::mcols(rawcdss)$gene),  newid)

  cdss = rawcdss
  cdss$transcript_id = newid
  cdss$temp_grouping_id = NULL
  cdss
}


.gnMappingHlpr = function(gnfeat, txfeatlst, knownlabs,
                          splcol = .getGeneIdLab(gnfeat, stoponfail=TRUE,
                                                 verbose=FALSE),
                          olaptype = "within", stopondup = TRUE) {
  unknown = which(is.na(knownlabs))
  ## is it already split?
  if(!methods::is(gnfeat, "GRangesList"))
    gnfeat = split(gnfeat, GenomicRanges::mcols(gnfeat)[[splcol]])
  hts = GenomicRanges::findOverlaps(gnfeat, txfeatlst, type=olaptype, select="all")
  ## duplicated queryHits is ok b/c gene can have more than one tx, right?
  subhits = S4Vectors::subjectHits(hts)
  if(anyDuplicated(subhits)) {
    duphits = which(duplicated(subhits))
    dupstart = start(IRanges::heads(txfeatlst[duphits], 1))

    if(stopondup)
      stop("mRNA feature(s) starting at [",
           paste(as.vector(dupstart), collapse=", "),
           "] appear to match more than one gene. If you feel the file is",
           " correct please contact the maintainer.")
    else { # throw away dup hits so they stay NA
      hts = hts[subhits %in% subhits[duphits],]
      subhits = S4Vectors::subjectHits(hts)
    }
  }
  ## don't override information that we already know.
  ## if efficiency of this routine ever becomes an issue
  ## this subsetting should be pushed up before the overlap
  ## calcs, but the indexing is simpler here and I don't
  ## think it will be too slow. We'll see ....

  keep = subhits %in% unknown

  knownlabs[ subhits[ keep ] ] = names(gnfeat)[ S4Vectors::queryHits(hts)[ keep ] ]
  knownlabs
}

assignTxsToGenes = function(rawtx, exons, genes) {

  txspl = split(rawtx, rawtx$temp_grouping_id)

  gnlabs = rep(NA_character_, times = length(txspl))

  ## exon mapping is more precise, try that first.
  ## exons should fall strictly within transcript sections
  ## not identical b/c of padding at ends. AFAIK should
  ## be identical for internal exons, but we aren't
  ## specifically checking for that here.

  gnlabs = .gnMappingHlpr(exons, txspl, gnlabs)
  if(anyNA(gnlabs)) {

    ## damn, now we have to try with genes
    gnlabs = .gnMappingHlpr(genes, txspl, gnlabs,
                            olaptype = "any",
                            stopondup=FALSE)

  }

  rawtx$gene_id= rep(gnlabs, times = lengths(txspl))
  rawtx

}

make_txgr = function(rawtxs, exons, sqinfo, genes) {
  if(length(rawtxs) > 0) {
    rawtxs = fill_stack_df(rawtxs, sqinfo = sqinfo)
    gvec = .getGeneIdVec(rawtxs)
    if(is.null(gvec)) {
      rawtxs = assignTxsToGenes(rawtxs, exons, genes)
      gvec = .getGeneIdVec(rawtxs)
    }

    rawtxs$tx_id_base = ifelse(is.na(gvec), paste("unknown_gene", cumsum(is.na(gvec)), sep="_"), gvec)
    spltxs = split(rawtxs, rawtxs$tx_id_base)
    txsraw = lapply(spltxs, function(grl) {
      grl$transcript_id = paste(grl$tx_id_base, (grl$temp_grouping_id -
                                                   min(grl$temp_grouping_id) +
                                                   1),
                                sep=".")
      grl$tx_id_base = NULL
      grl
    })
    txslst = GenomicRanges::GRangesList(txsraw)
    txs = unlist(txslst, use.names=FALSE)
    txs$gene_id = .getGeneIdVec(txs)
    txs$temp_grouping_id=NULL
  }  else if (length(exons) == 0L) {
    txs = GenomicRanges::GRanges(seqinfo=sqinfo)
  } else {
    message("No transcript features (mRNA) found, using spans of CDSs")
    spl = split(exons, exons$transcript_id %||% "unknown")
    # cat(summary(spl)) ## debug
    txslst = range(spl, ignore.strand  = F)
    # LIST = list(spl = spl, txslst = txslst)  ## debug
    # list2env(LIST, envir = .GlobalEnv) ## debug

    mcdf = GenomicRanges::mcols(unlist(IRanges::heads(spl, 1)))
    txs = unlist(txslst, use.names=FALSE)
    GenomicRanges::mcols(txs) = mcdf
    GenomeInfoDb::seqinfo(txs) = sqinfo
  }
  txs
}

make_varvr = function(rawvars, sq, sqinfo) {
  if(length(rawvars) == 0)
    return(VariantAnnotation::VRanges(seqinfo = sqinfo))
  if(is.null(sq)) {
    warning("importing variation features when origin sequence is not included in the file is not currently supported. Skipping ", length(rawvars), " variation features.")
    return(VariantAnnotation::VRanges(seqinfo = sqinfo))
  }

  vrs = fill_stack_df(rawvars, sqinfo = sqinfo)
  vrs$temp_grouping_id = NULL
  ## not all variants have /replace so we can't assume that ANY of them
  ## in a file have it (though  it is very likely at least one will)
  ## if none of them do, the column won't exist at all
  if(is.null(vrs$replace))
    vrs$replace = NA_character_

  ## makeVRangesFromGRanges seems to have a bug(?) that requires the
  ## columns used dfor the VRanges core info to be the first 6 in the
  ## granges mcols
  dels = nchar(vrs$replace) == 0L & !is.na(vrs$replace)
  if(length(dels)) {
    vrs[dels] = GenomicRanges::resize(vrs[dels], width(vrs[dels]) + 1L, fix = "end")
    vrs$replace[dels] = as.character(sq[GenomicRanges::resize(vrs[dels], 1L)])
  }
  newcols = S4Vectors::DataFrame( ref  = as.character(sq[vrs]),
                                  alt = vrs$replace,
                                  totalDepth = NA_integer_,
                                  altDepth = NA_integer_,
                                  refDepth = NA_integer_,
                                  sampleNames = NA_character_)
  GenomicRanges::mcols(vrs) = cbind(newcols, GenomicRanges::mcols(vrs))
  res  = VariantAnnotation::makeVRangesFromGRanges(vrs)
  res
}


make_exongr = function(rawex, cdss, sqinfo) {
  ##exns = sanitize_feats(rawex, exon_cols)
  exns = fill_stack_df(rawex, sqinfo = sqinfo)
  if(is.null(exns)) {

    message("No exons read from genbank file. Assuming sections of CDS are full exons")
    if(length(cdss) > 0) {
      exns = cdss
      exns$type = "exon"
    } else {
      return(GenomicRanges::GRanges(seqinfo = sqinfo))
    }

  } else {
    if(methods::is(exns, "GRangesList"))
      exns = utils::stack(exns)

    hits = GenomicRanges::findOverlaps(exns, cdss, type = "within", select = "all")
    qh = S4Vectors::queryHits(hits)
    qhtab = table(qh)
    dup = as.integer(names(qhtab[qhtab>1]))
    havdups = length(dup) > 0
    if(havdups) {
      ambig = exns[unique(dup)]
      exns = exns[-unique(dup)]
      noduphits = hits[-match(qh, dup)]
      warning("Some exons specified in genbank file have ambiguous relationship to transcript(s). ")
      ambig$transcript_id = paste(ambig$gene_id,
                                  "ambiguous", sep=".")
    } else {
      noduphits = hits
    }

    exns$transcript_id = cdss$transcript_id[S4Vectors::subjectHits(noduphits)]

  }
  exns$temp_grouping_id = NULL
  exns
}

make_genegr = function(x, sqinfo) {
  res = fill_stack_df(x, sqinfo = sqinfo)

  if(is.null(res))
    res = GenomicRanges::GRanges(seqinfo = sqinfo)
  else
    res$gene_id = .getGeneIdVec(res)

  res$temp_grouping_id = NULL
  res
}

gn_from_vers = function(x) {
  if(is.null(x))
    "unknown"
  else
    strsplit(x, " ")[[1]][1]
}

## super fast rbind of data.frame lists from Pete Haverty
.simple_rbind_dataframe <- function(dflist, element.colname) {
  numrows = vapply(dflist, nrow, integer(1))
  if (!missing(element.colname)) {
    list.name.col = factor(rep(names(dflist), numrows), levels=names(dflist))
  }
  dflist = dflist[ numrows > 0 ] # ARGH, if some data.frames have zero rows, factors become integers
  #  myunlist = function(x)  base::unlist(x, recursive=FALSE, use.names=FALSE)
  mylapply = base::lapply
  cn = names(dflist[[1]])
  inds = structure(1L:length(cn), names=cn)
  big <- mylapply(inds,
                  function(x) {
                    unlist(
                      # mylapply(dflist, function(y) { y[[x]] }),
                      mylapply(dflist, function(y) { .subset2(y, x) }),
                      recursive=FALSE, use.names=FALSE)
                  })
  if (!missing(element.colname)) {
    big[[element.colname]] = list.name.col
  }
  class(big) <- "data.frame"
  attr(big, "row.names") <- .set_row_names(length(big[[1]]))
  return(big)
}



cheap_unsafe_data.frame = function(..., lst = list(...)) {
  lens = lengths(lst)
  len = max(lens)
  if(!all(lens == len))
    lst = lapply(lst, rep, length.out = len)

  if(anyDuplicated.default(names(lst)))
    names(lst) =make.unique(names(lst))

  class(lst) = "data.frame"
  attr(lst, "row.names") = .set_row_names(length(lst[[1]]))
  lst
}


#### Classes
## File from: https://github.com/gmbecker/genbankr
#' @importClassesFrom Biostrings XStringSet
#' @importClassesFrom GenomicRanges GRangesList CompressedGRangesList
setOldClass("NULL")
setClassUnion("XStringSetOrNULL", c("XStringSet", "NULL"))

##' @title GenBank data objects.
##' @description These objects represent GenBank annotations
##'
##' @rdname GenBank-classes
##' @docType methods
setClass("GenBankRecord", slots = list(genes = "GenomicRanges", cds = "GenomicRanges",
                                       exons = "GenomicRanges",
                                       transcripts = "GenomicRanges",
                                       variations = "VRanges",
                                       sources = "GenomicRanges",
                                       other_features = "GenomicRanges",
                                       locus = "character",
                                       definition = "character",
                                       accession = "character",
                                       version = "character",
                                       source = "ANY",
                                       sequence = "XStringSetOrNULL"
)
)

##' @title GenBank File
##'
##' @description A resource class for use within the rtracklayer framework
##'
##' @rdname gbkfile
##' @docType methods
setClass("GenBankFile", contains = "GFFFile")

##' @rdname gbkfile
##' @docType methods
##' @aliases GBKFile-class
setClass("GBKFile", contains = "GenBankFile")

##' @rdname gbkfile
##' @docType methods
##' @aliases GBFile-class
setClass("GBFile", contains = "GenBankFile")

##' @title GBAccession ID class
##'
##' @description A class representing the (versioned) GenBank accession
##'
##' @rdname GBAccession
setClass("GBAccession", contains="character")

##' @rdname GBAccession
##' @param id A versioned GenBank Accession id
##' @return a \code{GBAccession} object.
GBAccession = function(id) {
  new("GBAccession", id)
}


### I'm leaving the old code below, remove eventually at a cleanup
# parse_genbank_file
#
# Function to parse the genbank file for the fasta alignment.
#
# importFrom genbankr readGenBank getSeq
#
# param gbk_path path to genbank file
# param g sequence length, available from the LDWeaver::parse_fasta_alignment() output (default = NULL),
# required if <length_check = T>
# param length_check specify whether to check if fasta and gbk sequence lengths are equal (default = T)
#
# return GenBankRecord object
#
# examples
# dontrun{
# aln_path <- system.file("extdata", "sample.aln", package = "LDWeaver")
# snp.dat <- parseFastaAlignment(aln_path, gap_freq = 0.15, maf_freq = 0.01, method = "default")
# gbk_path <- system.file("extdata", "sample.gbk", package = "LDWeaver")
# gbk <- parse_genbank_file(gbk_path, snp.dat$g)
# }
# export
# parse_genbank_file_old = function(gbk_path, g = NULL, length_check = T){
#   t0 = Sys.time()
#   # Check inputs
#   if(length_check){ # perform the length check
#     if(is.null(g)){
#       # then g cannot be null
#       stop("g must be provided to perform length check!\n")
#       # return(-1)
#     }
#   }

#   if(!file.exists(gbk_path)) stop(paste("Can't locate file", gbk_path))
#
#   # gbk = suppressWarnings(genbankr::import(gbk_path))
#   gbk = suppressWarnings(genbankr::readGenBank(gbk_path))
#   refseq = genbankr::getSeq(gbk)
#
#   if(length(refseq) != 1){
#     stop("The GBK file should contain the reference sequence!\n")
#     # return(-1)
#   }
#   # the length check is good to check if the alignment matches with the gbk, setting it to F will stop it
#   ref_g = length(refseq[[1]]) # length of the reference sequence
#
#   if(length_check){ # perform the length check
#     if(ref_g != g){ # perform check
#       stop("Genbank reference sequence length mismatches with the fasta alignment!\n")
#       # return(-1)
#     }
#
#   } else {
#     if(!is.null(g)){
#       if(ref_g != g){
#         warning("Fasta length does not match the genbank reference sequence length!\n")
#       }
#     } else {
#       warning("Similarity between the genbank reference and fasta sequences NOT checked, ignore if <pos> was provided...\n")
#     }
#   }
#   cat(paste("Successfully read gbk file:", gbk_path, "in", round(difftime(Sys.time(), t0, units = "secs"), 2), "s\n"))
#
#   # genbankr::seqinfo(gbk)
#   return(list(gbk = gbk,
#               ref_g = ref_g))
# }
