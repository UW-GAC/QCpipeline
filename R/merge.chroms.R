##########
# Merge chromosomes for association tests
# Usage: R --args config.file start end < merge.chroms.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "model_type", "gene_action")
optional <- c("annot_snp_filtCol", "annot_snp_rsIDCol", "maf.filter.type",
              "maf.absolute.threshold", "maf.linear.threshold", "maf.logistic.threshold",
              "out_assoc_prefix")
default <- c("quality.filter", "rsID", "snp.specific", 0.02, 30, 50, "assoc")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# chromosomes to merge
if (length(args) < 3) stop("missing start and end chromosomes")
start <- as.integer(args[2])
end <- as.integer(args[3])

# merge output 
# output file name example: study.model.1.additive.chr.24_24.RData
pathprefix <- config["out_assoc_prefix"]
actions <- unlist(strsplit(config["gene_action"]," "))
model.type <- unlist(strsplit(config["model_type"]," "))


for (i in 1:length(actions)) {
  cnt <- NULL
  tmp <- list()
  for (j in start:end) {
    fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".chr.",j,"_",j,".RData", sep="")
    dat <- NULL
    if (file.exists(fname)) {
      cnt <- c(cnt, j)
      dat <- getobj(fname)
    }
    tmp[[as.character(j)]] <- dat
  }
  combined <- do.call("rbind", tmp)

  # For SNPs which were monomorphic in either cases or controls, do Fisher Exact test
  warncol <- paste("model", i, actions[i], sep=".", "warningOrError")
  mono <- combined[[warncol]] %in% c(0,1)
  print(paste(sum(mono), "SNPs monomorphic in either cases or controls"))
  if (sum(mono) > 0) {
    message("Running Fisher exact test")
    fish <- assocTestFisherExact(combined[mono,])
    snpID <- combined$snpID
    # add only new columns
    dat <- merge(combined, fish[,c(1,6:13)], all.x=TRUE, sort=FALSE)
    # put snps back in original order
    combined <- dat[match(snpID, dat$snpID),]
  } 
  
  # since this is a single-model data frame, strip out unnecessary parts of column names
  names(combined) <- sub(paste("model", i, actions[i], "", sep="."), "", names(combined), fixed=TRUE)
  names(combined) <- sub(paste("model", i, "", sep="."), "", names(combined), fixed=TRUE)
  names(combined) <- sub(".G", "", names(combined), fixed=TRUE)
  
  # add rsID, chromosome, filters
  snpAnnot <- getobj(config["annot_snp_file"])
  snpID <- getSnpID(snpAnnot)
  stopifnot(all(combined$snpID %in% snpID))
  index <- match(combined$snpID,snpID)

  combined$rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"], index=index)
  #combined$chromosome <- getChromosome(snpAnnot, char=TRUE, index=index)
  combined$chromosome <- getChromosome(snpAnnot, index=index)
  combined$chromosome[combined$chromosome == 23] <- "X"
  combined$chromosome[combined$chromosome == 24] <- "XY"
  combined$chromosome[combined$chromosome == 25] <- "Y"
  combined$chromosome[combined$chromosome == 26] <- "M"
  combined$chromosome[combined$chromosome == 27] <- "U"
  combined$chromosome[combined$chromosome == 28] <- "XYY"
  combined$composite.filter <- getVariable(snpAnnot, config["annot_snp_filtCol"], index=index)

  # MAF filter
  maf.thresh <- switch(model.type[i],
                       linear=as.numeric(config["maf.linear.threshold"]),
                       logistic=as.numeric(config["maf.logistic.threshold"]))
  
  # snp-specific filtering for logistic is calculated the number of cases or controls, whichever is smaller
  # for linear it is all samples
  N <- switch(model.type[i],
              linear=max(combined$n, na.rm=TRUE),
              logistic=max( pmin( combined$nAA.cc0 + combined$nAB.cc0 + combined$nBB.cc0,
                                  combined$nAA.cc1 + combined$nAB.cc1 + combined$nBB.cc1), na.rm=TRUE))
  
  maf.filt <- switch(config["maf.filter.type"],
                     absolute=(!is.na(combined$MAF) &
                               combined$MAF > as.numeric(config["maf.absolute.threshold"])),
                     snp.specific=(!is.na(combined$MAF) & !is.na(combined$n) &
                                   2 * combined$MAF * (1 - combined$MAF) * N > maf.thresh))
  combined$comp.maf.filter <- combined$composite.filter & maf.filt
  
  print(dim(combined))
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  cat(paste("Saving", fname, "...\n"))
  save(combined, file=fname)
}

cat("Combined data contains output for the following chromosomes: \n")
cat(cnt, "\n")

