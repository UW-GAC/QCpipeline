##########
# Find variance in BAF for use in chromosome anomaly detection
#  and identifying low quality samples
# Usage: R --args config.file < baf_variance.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "nc_bl_file", "nc_geno_file")
optional <- c("annot_snp_missingCol", "out_baf_mean_file", "out_baf_med_file",
              "out_baf_sd_file", "snp_exclude_file")
default <- c("missing.n1", "baf_mean_by_scan_chrom.RData", "median_baf_sd_by_scan.RData",
             "baf_sd_by_scan_chrom.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

ncfile <- config["nc_geno_file"]
ncgds <- GenotypeReader(ncfile)
genoData <- GenotypeData(ncgds, snpAnnot=snpAnnot)

blfile <- config["nc_bl_file"]
blnc <- NcdfIntensityReader(blfile)
blData <- IntensityData(blnc, snpAnnot=snpAnnot)

# any snps to exclude?
if (!is.na(config["snp_exclude_file"])) {
  snp.exclude <- getobj(config["snp_exclude_file"])
  print(length(snp.exclude))
} else {
  snp.exclude <- c()
}
missing <- getVariable(snpAnnot, config["annot_snp_missingCol"])
snp.exclude <- union(snp.exclude, getSnpID(snpAnnot)[missing == 1])
print(length(snp.exclude))

baf <- sdByScanChromWindow(blData, genoData, var="BAlleleFreq",
                           snp.exclude=snp.exclude, return.mean=TRUE)
baf.sd <- baf$sd
save(baf.sd, file=config["out_baf_sd_file"])
baf.mean <- baf$mean
save(baf.mean, file=config["out_baf_mean_file"])

med.baf.sd <- medianSdOverAutosomes(baf.sd)
save(med.baf.sd, file=config["out_baf_med_file"])

close(blData)
close(genoData)
