##########
# Allele frequency
# Usage: R --args config.file < allele_freq.R
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

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

nc <- NcdfGenotypeReader(config["nc_geno_file"])
genoData <- GenotypeData(nc, scanAnnot=scanAnnot)

afreq <- alleleFrequency(genoData, scan.exclude=scan.exclude)
save(afreq, file=config["out_afreq_file"])

close(genoData)
