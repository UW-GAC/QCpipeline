##########
# Heterozygosity by SNP and sex
# Usage: R --args config.file < het_by_snp.R
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
optional <- c("scan_exclude_file", "out_het_file")
default <- c(NA, "het_by_snp.RData")
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

nc <- NcdfGenotypeReader(config["nc_geno_file"])
genoData <- GenotypeData(nc, scanAnnot=scanAnnot)

het <- hetBySnpSex(genoData, scan.exclude=scan.exclude)

save(het, file=config["out_het_file"])
close(genoData)
