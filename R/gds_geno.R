##########
# Convert netCDF to CoreArray GDS 
# Usage: R --args config.file < gds_geno.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "nc_geno_file", "gds_geno_file")
optional <- character()
default <- character()
config <- setConfigDefaults(config, required, optional, default)
print(config)

ncfile <- config["nc_geno_file"]
gdsfile <- config["gds_geno_file"]
snpAnnot <- getobj(config["annot_snp_file"])

convertNcdfGds(ncfile, gdsfile, snp.annot=snpAnnot)
if (!checkNcdfGds(ncfile, gdsfile)) stop("check failed")
