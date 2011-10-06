##########
# Mean intensity by scan and chromosome
# Usage: R --args config.file < mean_inten.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

ncfile <- config["nc_qxy_file"]
nc <- NcdfIntensityReader(ncfile)
intenData <- IntensityData(nc)

mninten <- meanIntensityByScanChrom(intenData)

save(mninten, file=config["out_inten_file"])
close(intenData)
