##########
# Heterozygosity by scan and chromosome
# Usage: R --args config.file < het_by_scan.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("nc_geno_file")
optional <- c("out_het_file")
default <- c("het_by_scan.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc)

het <- hetByScanChrom(genoData)

save(het, file=config["out_het_file"])
close(genoData)
