##########
# Convert netCDF to CoreArray GDS 
# Usage: R --args config.file gds_geno.R
##########

library(GWASTools)
library(gdsfmt)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

ncfile <- config["nc_geno_file"]
gdsfile <- config["gds_geno_file"]

convertNcdfGds(ncfile, gdsfile)
if (!checkNcdfGds(ncfile, gdsfile)) stop("check failed")
