##########
# Subset a netCDF file
# Usage: R --args config.file < ncdf_subset.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# get scans to include
scan.include <- getobj(config["scan_include_file"])

# create a new netCDF with samples in 'scan.include'
parent.ncdf <- config["nc_file"]
sub.ncdf <- config["nc_subset_file"]
ncdfSubset(parent.ncdf, sub.ncdf, sample.include=scan.include)

# check it against the source
ncdfSubsetCheck(parent.ncdf, sub.ncdf, sample.include=scan.include)
