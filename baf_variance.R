##########
# Find variance in BAF for use in chromosome anomaly detection
#  and identifying low quality samples
# Usage: R --args config.file < baf_variance.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, snpAnnot=snpAnnot)

blfile <- config["nc_bl_file"]
blnc <- NcdfIntensityReader(blfile)
blData <- IntensityData(blnc, snpAnnot=snpAnnot)

missing <- getVariable(snpAnnot, config["annot_snp_missingCol"])
snp.exclude <- getSnpID(snpAnnot)[missing == 1]

baf.sd <- sdByScanChromWindow(blData, genoData, var="BAlleleFreq", snp.exclude=snp.exclude)
save(baf.sd, file=config["out_baf_sd_file"])

med.baf.sd <- medianSdOverAutosomes(baf.sd)
save(med.baf.sd, file=config["out_baf_med_file"])

close(blData)
close(genoData)
