library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# merge output for chroms 1-26
# output file name example: study.model.1.additive.chr.24_24.RData
pathprefix <- config["assoc_output"]
actions <-  config["gene_action"]
actions <- unlist(strsplit(actions," "))


for (i in 1:length(actions))
{
  cnt <- NULL
  combined <- NULL
  for (j in (args[2]:args[3]))
  {
    fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".chr.",j,"_",j,".RData", sep="")
    dat <- NULL
    if (file.exists(fname))
    {
      cnt <- c(cnt, j)
      dat <- getobj(fname)
    }
    combined <- rbind(combined, dat)
  }
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.RData", sep="")

  cat(paste("Saving", fname, "...\n"))
  save(combined, file=fname)
}

cat("Combined data contains output for the following chromosomes: \n")
cat(cnt, "\n")

