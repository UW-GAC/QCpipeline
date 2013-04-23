##########
# Allelic frequency batch test
# Usage: R --args config.file < batch_chisq.R
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
optional <- c("annot_scan_batchCol", "annot_scan_hapmapCol", "out_fisher_file", "scan_exclude_file")
default <- c("Sample.Plate", "geno.cntl", "batch_fisher", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

data <- GenotypeReader(config["nc_geno_file"])
genoData <- GenotypeData(data, scanAnnot=scanAnnot)

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

# exclude hapmaps from allele frequency test
if (!is.na(config["annot_scan_hapmapCol"])) {
  hapmap <- getVariable(scanAnnot, config["annot_scan_hapmapCol"])
  hm.ids <- scanID[hapmap == 1]
  scan.exclude <- union(scan.exclude, hm.ids)
}
length(scan.exclude)

batchFisherTest(genoData, batchVar=config["annot_scan_batchCol"],
               scan.exclude=scan.exclude, return.by.snp=TRUE,
               outfile=config["out_fisher_file"])

close(genoData)
