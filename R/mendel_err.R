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
required <- c("annot_scan_file", "annot_snp_file", "geno_subj_file")
optional <- c("annot_scan_subjectCol", "annot_snp_missingCol",
              "mend_scan_exclude_file", "out_mend_file")
default <- c("subjectID", "missing.n1", NA, "mendel_err.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

data <- GenotypeReader(config["geno_subj_file"])
scanAnnot <- getobj(config["annot_scan_file"])
stopifnot(all(hasVariable(scanAnnot, c("family", "father", "mother", "sex"))))
# take subset of annotation to match netCDF
scanAnnot <- scanAnnot[match(getScanID(data), getScanID(scanAnnot)), ]
genoData <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
scanID <- getScanID(genoData)


# are there any scans to exclude?
if (!is.na(config["mend_scan_exclude_file"])) {
  scan.exclude <- getobj(config["mend_scan_exclude_file"])
  #stopifnot(all(scan.exclude %in% scanID))
  scanAnnot <- scanAnnot[!(scanID %in% scan.exclude),]
}
nrow(scanAnnot)

annot <- getVariable(scanAnnot, c("family", config["annot_scan_subjectCol"],
                                  "father", "mother", "sex", "scanID"))
names(annot)[2] <- "subjectID"
men.list <- with(annot, mendelList(family, subjectID,
                                   father, mother, sex, scanID))

mendel <- mendelErr(genoData, men.list)

# turn mendel$snp into a data frame instead of a list
stopifnot(allequal(snpID, names(mendel$snp$error.cnt)))
stopifnot(allequal(snpID, names(mendel$snp$check.cnt)))
snp <- data.frame("snpID"=snpID, "error.cnt"=mendel$snp$error.cnt,
                  "check.cnt"=mendel$snp$check.cnt)
# set missing SNPs to NA
miss <- getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1
snp$error.cnt[miss] <- NA
snp$check.cnt[miss] <- NA
mendel$snp <- snp

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
sum(sel, na.rm=TRUE)
sum(mendel.err.cnt[sel]>0, na.rm=TRUE)
sum(mendel.err.cnt[sel]>0, na.rm=TRUE)/sum(sel, na.rm=TRUE)
sum(mendel.err.cnt[sel]>1, na.rm=TRUE)
sum(mendel.err.cnt[sel]>1, na.rm=TRUE)/sum(sel, na.rm=TRUE)
sum(mendel.err.cnt[sel]>2, na.rm=TRUE)
sum(mendel.err.cnt[sel]>2, na.rm=TRUE)/sum(sel, na.rm=TRUE)
max(mendel.err.cnt[sel], na.rm=TRUE)
