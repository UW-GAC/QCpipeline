##########
# Duplicate discordance
# Usage: R --args config.file < dup_disc.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "geno_file")
optional <- c("annot_scan_subjectCol", "annot_snp_missingCol", 
              "out_afreq_file", "out_disc_maf_file", "disc_scan_exclude_file")
default <- c("subjectID", "missing.n1", "allele_freq.RData", "dup_disc_maf.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

data <- GenotypeReader(config["geno_file"])
genoData <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# are there any scans to exclude?
if (!is.na(config["disc_scan_exclude_file"])) {
  scan.exclude <- getobj(config["disc_scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

# exclude missing SNPs
snp.exclude <- snpID[getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1]
length(snp.exclude)

afreq <- getobj(config["out_afreq_file"])
disc <- duplicateDiscordance(genoData, subjName.col=config["annot_scan_subjectCol"],
                             scan.exclude=scan.exclude, snp.exclude=snp.exclude,
                             minor.allele.only=TRUE, allele.freq=afreq[,"all"])

# give output data frame the same dimensions as snp annotation
snp <- merge(data.frame(snpID), disc$discordance.by.snp, all.x=TRUE)
disc$discordance.by.snp <- snp

save(disc, file=config["out_disc_maf_file"])
