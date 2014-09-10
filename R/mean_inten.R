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
required <- c("xy_file")
optional <- c("out_inten_file")
default <- c("mean_inten.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

xy <- IntensityReader(config["xy_file"])
intenData <- IntensityData(xy)

mninten <- meanIntensityByScanChrom(intenData)

save(mninten, file=config["out_inten_file"])
close(intenData)
