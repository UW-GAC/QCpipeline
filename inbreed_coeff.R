##########
# Individual inbreeding coefficient
# Usage: R --args config.file < inbreed_coeff.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

nc <- NcdfGenotypeReader(config["nc_geno_file"])
genoData <- GenotypeData(nc, scanAnnot=scanAnnot)
snpID <- getSnpID(genoData)

# SNP filter - contain the condition of every two SNPs being at least 15kb apart
# pre-selected SNPs: 15kb apart, autosomal, MAF > 0%, and MCR < 5%
# also select SNPs with MAF > 5%
snp.ids <- getobj(config["out_snp_file"])
length(snp.ids)
afreq <- getobj(config["out_afreq_file"])
stopifnot(allequal(rownames(afreq), snpID))
maf.filt <- !is.na(afreq[,"all"]) & afreq[,"all"] < 0.95 & afreq[,"all"] > 0.05
snp.sel <- (snpID %in% snp.ids) & maf.filt
table(snp.sel)

# to store output
inbrd.coeff <- rep(NA, length(scanID))
names(inbrd.coeff) <- scanID

# loop through scans
N <- length(scanID)
cnt <- 0
block <- 100
for(i in 1:(N%/%block + 1)){ # read a block of scans
  message(paste("reading block", i, "of", (N%/%block + 1)))
  message(paste("start scan:", start <- (i-1)*block + 1, "\n"))
  message(paste("end scan:", end <- min(i*block, N), "\n"))
  geno <- getGenotype(genoData, snp=c(1,-1), scan=c(start,min(block, end-start+1)))

  for (j in 1:ncol(geno)) { # traverse scans read from ncdf
        
        ge <- geno[,j,drop=FALSE] # scan j in block i
        nona <- !is.na(ge) # j can have extra missing genotypes. exclude them

        filt <- snp.sel & nona 
        
        ge <- ge[filt,1, drop=FALSE]
        afreq.filt <- afreq[filt,"all", drop=FALSE]
        
        # 3 terms in numerator
        term1 <- ge * ge
        term2 <- (1 + 2*afreq.filt) * ge
        term3 <- 2*(afreq.filt * afreq.filt)
        numer <- term1 - term2 + term3
        
        # denominator
        denom <- 2*afreq.filt*(1-afreq.filt)

        # ratio of averages - Weir
        # sum over SNPs for each scan
        nume <- sum(numer)
        deno <- sum(denom)
        coeff <- nume/deno

        cnt <- cnt + 1
        inbrd.coeff[cnt] <- coeff       
      }
}

save(inbrd.coeff, file=config["out_inbrd_file"])
