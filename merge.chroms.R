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
  combined <- NULL
  for (j in 1:26)
  {
    fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".chr.",j,"_",j,".RData", sep="")
    dat <- getobj(fname)
    combined <- rbind(combined, dat)
  }
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.RData", sep="")
  print(fname)
  save(combined, file=fname)
}


