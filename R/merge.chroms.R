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
required <- c("annot_snp_file", "out_assoc_prefix", "gene_action")
optional <- c("annot_snp_filtCol", "annot_snp_rsIDCol", "maf.filter")
default <- c("quality.filter", "rsID", 0.02)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# chromosomes to merge
if (length(args) < 3) stop("missing start and end chromosomes")
start <- as.integer(args[2])
end <- as.integer(args[3])

# merge output 
# output file name example: study.model.1.additive.chr.24_24.RData
pathprefix <- config["out_assoc_prefix"]
actions <-  config["gene_action"]
actions <- unlist(strsplit(actions," "))


for (i in 1:length(actions))
{
  cnt <- NULL
  tmp <- list()
  for (j in start:end)
  {
    fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".chr.",j,"_",j,".RData", sep="")
    dat <- NULL
    if (file.exists(fname))
    {
      cnt <- c(cnt, j)
      dat <- getobj(fname)
    }
    tmp[[as.character(j)]] <- dat
  }
  combined <- do.call("rbind", tmp)

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
  combined$chromosome <- getChromosome(snpAnnot, char=TRUE, index=index)
  combined$quality.filter <- getVariable(snpAnnot, config["annot_snp_filtCol"], index=index)
  combined$qual.maf.filter <- combined$quality.filter & (!is.na(combined$minor.allele)) & combined$MAF>as.numeric(config["maf.filter"])
  
  print(dim(combined))
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  cat(paste("Saving", fname, "...\n"))
  save(combined, file=fname)
}

cat("Combined data contains output for the following chromosomes: \n")
cat(cnt, "\n")

