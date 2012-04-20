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
required <- c("assoc_output", "gene_action")
optional <- NULL
default <- NULL
config <- setConfigDefaults(config, required, optional, default)
print(config)

# chromosomes to merge
if (length(args) < 3) stop("missing start and end chromosomes")
start <- as.integer(args[2])
end <- as.integer(args[3])

# merge output 
# output file name example: study.model.1.additive.chr.24_24.RData
pathprefix <- config["assoc_output"]
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
  print(dim(combined))
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.RData", sep="")

  cat(paste("Saving", fname, "...\n"))
  save(combined, file=fname)
}

cat("Combined data contains output for the following chromosomes: \n")
cat(cnt, "\n")

