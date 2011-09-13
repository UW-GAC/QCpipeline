##########
# Mean intensity by scan and chromosome
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

ncfile <- config["nc_qxy_file"]
nc <- NcdfIntensityReader(ncfile)
intenData <- IntensityData(nc)

mninten <- meanIntensityByScanChrom(intenData)

save(mninten, file=config["out_inten_file"])
close(intenData)
