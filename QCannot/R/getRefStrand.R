# R function to take BLAT results for a set of array variants and determine the RefStrand column
# assumes BLAT has been run on TopGenomicSeq
# assumes build of Illumina SNP annotation matches build used to run BLAT
# S Nelson, UWGCC, 9/22/2015

getRefStrand <- function(snp.dat, blat.rslt, metric.sel="score") {

  # check that snp.dat argument appears to be Illumina manifest
  cols.needed <- c("Name","Chr","MapInfo", "IlmnStrand")
  cols.missing <- setdiff(cols.needed, names(snp.dat))
  if(length(cols.missing > 0)) {
    stop("SNP input is missing columns(s) ",paste(cols.missing, collapse=" "),"\n")
  }

  snp.dat <- snp.dat[,cols.needed]

  # check that blat.rslt argument appears to be output of runBlat function
  cols.needed <- c("Q_name","T_name","strand", "T_start","T_end")
  cols.missing <- setdiff(cols.needed, names(blat.rslt))
  if(length(cols.missing > 0)) {
    stop("BLAT result is missing columns(s) ",paste(cols.missing, collapse=" "),"\n")
  }

  # report info about BLAT results
  n.blat.vars <- prettyNum(length(unique(blat.rslt$Q_name)), big.mark=",") 
  n.blat <- prettyNum(nrow(blat.rslt), big.mark=",")
  message("BLAT result: ", n.blat," matches for ", n.blat.vars," total variants\n")

  # how many entries in snp.dat have BLAT match in BLAT result, initially?
  inp.vars <- unique(snp.dat$Name)
  n.inp <- prettyNum(length(inp.vars), big.mark=",")
  n.miss <- prettyNum(sum(!is.element(inp.vars, blat.rslt$Q_name)), big.mark=",")
  message("Of ", n.inp, " variants in 'snp.dat' input, ", n.miss," have no BLAT results in 'blat.rslt'\n") 

  # adjust snp.dat chromosome codes to match BLAT
  snp.dat$Chr.match <- snp.dat$Chr
  snp.dat$Chr.match[snp.dat$Chr %in% "MT"] <- "M"
  # for XY SNPs, BLAT returns a match on both X and Y - here we'll get X chrom coordinates
  snp.dat$Chr.match[snp.dat$Chr %in% "XY"] <- "X"  
  
  # narrow BLAT results down to correct location
  blat.rslt$chr <- sub("chr","", blat.rslt$T_name)

  # use Q_name to merge expected chrom and position
  blat.rslt$chr_exp <- snp.dat$Chr.match[match(blat.rslt$Q_name, snp.dat$Name)]
  blat.rslt$pos_exp <- snp.dat$MapInfo[match(blat.rslt$Q_name, snp.dat$Name)]

  # remove where match is on wrong chrom
  rm <- blat.rslt$chr_exp != blat.rslt$chr
  n.rm <- prettyNum(sum(rm), big.mark=",")
  message("Removing ", n.rm, " BLAT matches on wrong chromosome\n")
  blat.rslt <- blat.rslt[!rm,]

  # flag where SNP location is within BLAT match
  blat.rslt$snp.in.blat.match <- blat.rslt$pos_exp >= blat.rslt$T_start &
                                 blat.rslt$pos_exp <= blat.rslt$T_end

  # remove where SNP location is not within BLAT match
  n.rm <- prettyNum(sum(!blat.rslt$snp.in.blat.match), big.mark=",")
  message("Removing ", n.rm, " BLAT matches not overlapping variant position\n")
  blat.rslt <- blat.rslt[blat.rslt$snp.in.blat.match,]

  # set aside variants that now have a single BLAT match
  dup.vars <- unique(blat.rslt$Q_name[duplicated(blat.rslt$Q_name)])
  blat.r1 <- blat.rslt[!is.element(blat.rslt$Q_name, dup.vars),]

  # for remaining variants, choose based on user-specified input
  # either max score or max pct.identity
  blat.dup <- blat.rslt[is.element(blat.rslt$Q_name, dup.vars),]

  # sort by selected metric - highest to lowest
  met <- blat.dup[,metric.sel]
  blat.srt <- blat.dup[order(blat.dup$Q_name, -met), ]

  # choose the first record
  blat.r2 <- blat.srt[!duplicated(blat.srt$Q_name),]

  # combine rounds 1 and 2 of BLAT matches
  blat.keep <- c("Q_name","strand","id")
  blat.final <- rbind(blat.r1[,blat.keep], blat.r2[,blat.keep])
  names(blat.final)[2:3] <- c("blat.strand","blat.id")

  # add blat strand to snp.dat
  snp.ann <- merge(snp.dat[,c("Name","Chr","MapInfo","IlmnStrand")],
                   blat.final, by.x="Name", by.y="Q_name",
                   all=TRUE, sort=FALSE)

  # determine RefStrand from relationships of IlmnStrand and BLAT strand
  snp.ann$RefStrand <- NA
  snp.ann$RefStrand[is.element(snp.ann$IlmnStrand,c("TOP","PLUS")) &
                   snp.ann$blat.strand %in% "+"] <- "+"
  snp.ann$RefStrand[is.element(snp.ann$IlmnStrand,c("TOP","PLUS")) &
                   snp.ann$blat.strand %in% "-"] <- "-"
  snp.ann$RefStrand[is.element(snp.ann$IlmnStrand,c("BOT","MINUS")) &
                   snp.ann$blat.strand %in% "+"] <- "-"
  snp.ann$RefStrand[is.element(snp.ann$IlmnStrand,c("BOT","MINUS")) &
                   snp.ann$blat.strand %in% "-"] <- "+"

  # product - return input snp.dat, row for row, with "RefStrand" and "BLAT id"
  return(snp.ann)
}
