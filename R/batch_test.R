##########
# Allelic frequency batch test
# Usage: R --args config.file test.type < batch_test.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "geno_file")
optional <- c("annot_scan_batchCol", "annot_scan_hapmapCol", "out_batch_prefix", "scan_exclude_file")
default <- c("Sample.Plate", "geno.cntl", "batch_test", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# check for test
if (length(args) < 2) stop("missing test type (chisq or fisher)")
type <- args[2]
  
(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

data <- GenotypeReader(config["geno_file"])
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

if (type == "chisq") {
    batchChisqTest(genoData, batchVar=config["annot_scan_batchCol"],
                   scan.exclude=scan.exclude, return.by.snp=TRUE,
                   outfile=config["out_batch_prefix"])
} else if (type == "fisher") {
    batchFisherTest(genoData, batchVar=config["annot_scan_batchCol"],
                    scan.exclude=scan.exclude, return.by.snp=TRUE,
                    outfile=config["out_batch_prefix"])
} else {
    stop("test type must be chisq or fisher")
}

close(genoData)
