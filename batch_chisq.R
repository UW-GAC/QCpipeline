##########
# Allelic frequency batch test
# Usage: R --args config.file < batch_chisq.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]
print(config.table)

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

batchChisqTest(genoData, batchVar=config["annot_scan_batchCol"],
               scan.exclude=scan.exclude, return.by.snp=TRUE,
               outfile=config["out_chisq_file"])

close(genoData)
