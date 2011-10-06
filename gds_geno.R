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
print(config)

ncfile <- config["nc_geno_file"]
gdsfile <- config["gds_geno_file"]
snpAnnot <- getobj(config["annot_snp_file"])
scanAnnot <- getobj(config["annot_scan_file"])

convertNcdfGds(ncfile, gdsfile, sample.annot=pData(scanAnnot),
               snp.annot=pData(snpAnnot),
               rsID.col=config["annot_snp_rsIDCol"],
               alleleA.col=config["annot_snp_alleleACol"],
               alleleB.col=config["annot_snp_alleleBCol"])
if (!checkNcdfGds(ncfile, gdsfile)) stop("check failed")
