##########
# Convert netCDF to CoreArray GDS 
# Usage: R --args config.file < gds_geno.R
##########

library(GWASTools)
library(QCpipeline)
library(gdsfmt)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "nc_geno_file", "gds_geno_file")
optional <- character()
default <- character()
config <- setConfigDefaults(config, required, optional, default)
print(config)

ncfile <- config["nc_geno_file"]
gdsfile <- config["gds_geno_file"]
snpAnnot <- getobj(config["annot_snp_file"])
scanAnnot <- getobj(config["annot_scan_file"])
# take subset of annotation to match netCDF
nc <- NcdfGenotypeReader(ncfile)
scanAnnot <- scanAnnot[match(getScanID(nc), getScanID(scanAnnot)),]
close(nc)

convertNcdfGds(ncfile, gdsfile, sample.annot=scanAnnot, snp.annot=snpAnnot)
if (!checkNcdfGds(ncfile, gdsfile)) stop("check failed")
