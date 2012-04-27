##########
# Mendelian errors
# Usage: R --args config.file < mendel_err.R
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
optional <- c("annot_scan_subjectCol", "mend_scan_exclude_file", "out_mend_file")
default <- c("subjectID", NA, "mendel_err.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
stopifnot(all(hasVariable(scanAnnot, c("family", "father", "mother", "sex"))))
scanID <- getScanID(scanAnnot)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot)

# are there any scans to exclude?
if (!is.na(config["mend_scan_exclude_file"])) {
  scan.exclude <- getobj(config["mend_scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
  scanAnnot <- scanAnnot[!(scanID %in% scan.exclude),]
}
nrow(scanAnnot)

annot <- getVariable(scanAnnot, c("family", config["annot_scan_subjectCol"],
                                  "father", "mother", "sex", "scanID"))
names(annot)[2] <- "subjectID"
men.list <- with(annot, mendelList(family, subjectID,
                                   father, mother, sex, scanID))

mendel <- mendelErr(genoData, men.list)
save(mendel, file=config["out_mend_file"])

# results
# "all.trios" is the count for all combinations of the three members of a trio (many because of multiple scans for each subject)
# "trios" takes the mean over all possible trios for a given family to give one value for each family
dim(mendel$trios)
dim(mendel$all.trios)

# get error rates
mendel$trios$Men.err.rate <- mendel$trios$Men.err.cnt/mendel$trios$Men.cnt
mendel$trios$mtDNA.err.rate <- mendel$trios$mtDNA.err/mendel$trios$mtDNA.cnt
mendel$trios[,c("fam.id", "Men.err.rate", "mtDNA.err.rate")]

mendel.err.cnt <- mendel$snp$error.cnt
mendel.sampsize <- mendel$snp$check.cnt
mendel.err.rate <- mendel.err.cnt / mendel.sampsize
mean(mendel.err.rate, na.rm=TRUE)

# what number and fraction of SNPs have a non-zero error rate
sel <- mendel.sampsize > 0
sum(sel)
sum(mendel.err.cnt[sel]>1)
sum(mendel.err.cnt[sel]>1)/sum(sel)
sum(mendel.err.cnt[sel]>2)
sum(mendel.err.cnt[sel]>2)/sum(sel)
