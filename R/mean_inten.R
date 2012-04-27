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

# check config and set defaults
required <- c("nc_xy_file")
optional <- c("out_inten_file")
default <- c("mean_inten.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

ncfile <- config["nc_xy_file"]
nc <- NcdfIntensityReader(ncfile)
intenData <- IntensityData(nc)

mninten <- meanIntensityByScanChrom(intenData)

save(mninten, file=config["out_inten_file"])
close(intenData)
