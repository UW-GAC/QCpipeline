##########
# merge HWE results to single data frame
# Usage: R --args config.file start end by < hwe_merge.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("out_hwe_prefix")
optional <- NULL
default <- NULL
config <- setConfigDefaults(config, required, optional, default)
print(config)

if (length(args) < 4) stop ("missing start end by")
start <- as.integer(args[2])
end <- as.integer(args[3])
by <- as.integer(args[4])

cstart <- seq(start, end, by)
cend <- cstart + by - 1
if (cend[length(cend)] > 23) cend[length(cend)] <- 23

tmp <- list()
for (i in 1:length(cstart)) {
  file <- paste(config["out_hwe_prefix"], ".chr.", cstart[i], "_", cend[i], ".RData", sep="")
  print(file)
  tmp[[i]] <- getobj(file)
}
hwe <- do.call("rbind", tmp)
dim(hwe)

save(hwe, file=paste(config["out_hwe_prefix"], "RData", sep="."))
