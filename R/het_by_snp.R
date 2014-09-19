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
required <- c("annot_scan_file", "geno_file")
optional <- c("scan_exclude_file", "out_het_file")
default <- c(NA, "het_by_snp.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)


# netcdf should be subject-level
(data <- GenotypeReader(config["geno_file"]))
scanAnnot <- getobj(config["annot_scan_file"])
# take subset of annotation to match netCDF
scanAnnot <- scanAnnot[match(getScanID(data), getScanID(scanAnnot)),]
genoData <- GenotypeData(data, scanAnnot=scanAnnot)
scanID <- getScanID(genoData)



# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  #stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)


het <- hetBySnpSex(genoData, scan.exclude=scan.exclude)

save(het, file=config["out_het_file"])
close(genoData)
