##########
# Allele frequency
# Usage: R --args config.file < allele_freq.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "nc_geno_file")
optional <- c("out_afreq_file", "scan_exclude_file")
default <- c("allele_freq.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

data <- GenotypeReader(config["nc_geno_file"])
genoData <- GenotypeData(data, scanAnnot=scanAnnot)

afreq <- alleleFrequency(genoData, scan.exclude=scan.exclude)
save(afreq, file=config["out_afreq_file"])

close(genoData)
