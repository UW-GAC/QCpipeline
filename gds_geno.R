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
optional <- c("annot_snp_rsIDCol", "annot_snp_alleleACol", "annot_snp_alleleBCol")
default <- c(NA, NA, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

ncfile <- config["nc_geno_file"]
gdsfile <- config["gds_geno_file"]
snpAnnot <- getobj(config["annot_snp_file"])
scanAnnot <- getobj(config["annot_scan_file"])

if (!is.na(config["annot_snp_rsIDCol"]) &
    !hasVariable(snpAnnot, config["annot_snp_rsIDCol"])) {
  warning(paste(config["annot_snp_rsIDCol"], "not found in",
                config["annot_snp_file"]))
}
if (!is.na(config["annot_snp_alleleACol"]) &
    !hasVariable(snpAnnot, config["annot_snp_alleleACol"])) {
  warning(paste(config["annot_snp_alleleACol"], "not found in",
                config["annot_snp_file"]))
}
if (!is.na(config["annot_snp_alleleBCol"]) &
    !hasVariable(snpAnnot, config["annot_snp_alleleBCol"])) {
  warning(paste(config["annot_snp_alleleBCol"], "not found in",
                config["annot_snp_file"]))
}

convertNcdfGds(ncfile, gdsfile, sample.annot=pData(scanAnnot),
               snp.annot=pData(snpAnnot),
               rsID.col=config["annot_snp_rsIDCol"],
               alleleA.col=config["annot_snp_alleleACol"],
               alleleB.col=config["annot_snp_alleleBCol"])
if (!checkNcdfGds(ncfile, gdsfile)) stop("check failed")
