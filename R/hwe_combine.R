##########
# merge HWE results to single data frame
# Usage: R --args config.file start end < hwe_combine.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# chromosomes to merge
if (length(args) < 3) stop("missing start and end chromosomes")
start <- as.integer(args[2])
end <- as.integer(args[3])

# check config and set defaults
required <- c()
optional <- c("out_hwe_prefix")
default <- c("hwe")
config <- setConfigDefaults(config, required, optional, default)
print(config)

pathprefix <- config["out_hwe_prefix"]
tmp <- list()
for (j in start:end) {
    fname <- paste0(pathprefix, "_chr", j, ".RData")
    if (file.exists(fname)) {
        tmp[[as.character(j)]] <- getobj(fname)
    }
}
hwe <- do.call("rbind", tmp)
rm(tmp)

save(hwe, file=paste0(config["out_hwe_prefix"], ".RData"))
