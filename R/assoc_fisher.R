##########
# Fisher exact test for SNPs monomorphic in cases or controls
# Usage: R --args config.file < assoc_fisher.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "model_type", "geno_file", "outcome")
optional <- c("out_assoc_prefix", "scan_exclude")
default <- c("assoc", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

stopifnot(config["model_type"] == "logistic")

# make genotypedata
scanAnnot <- getobj(config["annot_scan_file"])

data <- GenotypeReader(config["geno_file"])
sid <- getScanID(data)
scanAnnot <- scanAnnot[match(sid, scanAnnot$scanID),] 
stopifnot(all(scanAnnot$scanID==sid))

genoData <- GenotypeData(data, scanAnnot=scanAnnot)
  
# outcome variable
(outcome <- config["outcome"])

# scan.exclude
if (!is.na(config["scan_exclude"])) {
  scan.exclude <- getobj(config["scan_exclude"])
} else {
  scan.exclude <- NULL
}

# find allele frequency in cases and controls
scanID <- getScanID(genoData)
cc <- getScanVariable(genoData, outcome) 
cc0.exclude <- union(scanID[cc == 1], scan.exclude)
cc1.exclude <- union(scanID[cc == 0], scan.exclude)

cc0.freq <- alleleFrequency(genoData, scan.exclude=cc0.exclude)[,"MAF"]
cc1.freq <- alleleFrequency(genoData, scan.exclude=cc1.exclude)[,"MAF"]

# identify SNPs monomorphic in cases or controls
cc.mono <- (!is.na(cc0.freq) & !is.na(cc1.freq)) &
    ((cc0.freq == 0 & cc1.freq > 0) | (cc0.freq > 0 & cc1.freq == 0))
if (sum(cc.mono) == 0) {
  message("No SNPs monomorphic in cases or controls.")
  q(save="no", status=0)
}

# snp.include
# can't do Fisher test for sex chromosomes on males and females combined
snp.include <- getSnpID(genoData, index=(cc.mono & getChromosome(genoData) %in% 1:22))

out <- batchFisherTest(genoData,
                       batchVar=outcome,
                       snp.include=snp.include,
                       scan.exclude=scan.exclude,
                       return.by.snp=TRUE,
                       conf.int=TRUE)

dat <- data.frame(snpID = as.integer(rownames(out$pval)),
                  Fisher.OR = out$oddsratio[,1],
                  Fisher.OR_L95 = out$confint.low[,1],
                  Fisher.OR_U95 = out$confint.high[,1],
                  Fisher.pval = out$pval[,1],
                  nA = out$allele.counts[,"nA"],
                  nA = out$allele.counts[,"nB"])
            
save(dat, file=paste0(config["out_assoc_prefix"], "_fisher.RData"))

close(genoData)
