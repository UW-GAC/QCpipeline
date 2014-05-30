# R function to convert genome build in a SNP annotation
# uses 'liftOver' function from R package 'rtracklayer'
# SN, 5/23/2014

# Two arguments:
# 1) Open SNP annotation object (object of class 'SnpAnnotationDataFrame') with "chromosome" and "position" columns in original build
# 2) Path to UCSC chain file, from original build to target build.
# e.g., "/projects/geneva/gcc-fs2/SNP_annotation/UCSC_downloads/hg19/liftOver/hg19ToHg18.over.chain"
# 3) Verbose, logical to report out counts of successful and unsuccessful conversions

# Returns input SNP annotation as a data frame with added columns: chromosome and position in the target build

## example usage:
## library(OLGApipline)
## data.dir <- "/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1"
## olgaData <- OlgaGenotypeData(data.dir)
## snpAnnot <- getSnpAnnotation(olgaData, 21)
## chain.file <- "/projects/geneva/gcc-fs2/SNP_annotation/UCSC_downloads/hg19/liftOver/hg19ToHg18.over.chain"
## converted <- convertBuild(snpAnnot, chain.file)

###########################################

convertBuild <- function(snpAnnot,chain.file, verbose=TRUE)
  {
    # check for required fields in SNP annotation
   stopifnot(all(c("chromosome","position","rsID") %in% varLabels(snpAnnot)))

   # import chain file
   chain <- import.chain(chain.file)

   # create data frame of SNP annotation
   snp <- pData(snpAnnot)

   # allows for more than one input chromosome
   chroms <- snp$chromosome
   
   # adjust chromosome names to characters for non-autosomal SNPs
   chroms[XchromCode(snpAnnot) == snp$chromosome] <- "X"
   chroms[YchromCode(snpAnnot) == snp$chromosome] <- "Y"
   chroms[MchromCode(snpAnnot) == snp$chromosome] <- "M" 

   # create seqNames for GRanges
   chrGr <- paste("chr",chroms,sep="")

   # create IRanges object
   ir <- IRanges(start=snp$position, width=1, names=snp$rsID)
   
   # create GRanges object
   gr <-  GRanges(seqnames=chrGr,ranges=ir)

   # run liftOver function
   lift <- liftOver(gr, chain)

   if(length(lift)!=nrow(snp))
     stop("Length of GRangesList output from liftOver does not match number of input records")

   # check that elements length is never >1 (i.e. more than one mapping)  -- don't expect >1, but that would be a problem
   mult.rslt <- sum(elementLengths(lift)>1)
   if(mult.rslt>0)
     stop(mult.rslt," input positions map to more than one position in target build\n")

   # report out number of successful liftovers
   if(verbose){
   one.rslt <- sum(elementLengths(lift)==1)
   no.rslt <- sum(elementLengths(lift)==0)
   message(one.rslt," input positions successfully converted; ",no.rslt," positions failed conversion (will be set to NA)")
   pct <- round(one.rslt/nrow(snp),6)*100
   message(pct,"% successful conversions")
 }

   # extract converted chrom
   lift.chr <- as.character(seqnames(lift))

   # strip 'chr' and convert to integer code
   chr.tmp <-  gsub("chr","", lift.chr)
   chr.tmp[is.element(chr.tmp,"X")] <- XchromCode(snpAnnot)
   chr.tmp[is.element(chr.tmp,"Y")] <- YchromCode(snpAnnot)
   chr.tmp[is.element(chr.tmp,"M")] <- MchromCode(snpAnnot)
   chr.int <- as.integer(chr.tmp)

   # extract converted base pair position
   lift.start <- as.integer(start((lift)))
   lift.end <- as.integer(end((lift)))

   # start should equal end, but verify and throw warning if not
   if(!allequal(lift.start, lift.end))
     warning("Converted start and end positions are not equal - returning both")

   # merge with input
   snp$chromosome.converted <- chr.int
   snp$position.converted <- lift.start

   # if end and start were not equal, also return end
   if(!allequal(lift.start, lift.end))
     {
       snp$position.end.converted <- lift.end
     }

   # return build conversions
   return(snp)

  }
