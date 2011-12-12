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
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot)

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

batchFisherTest(genoData, batchVar=config["annot_scan_batchCol"],
               scan.exclude=scan.exclude, return.by.snp=TRUE,
               outfile=config["out_fisher_file"])

close(genoData)
