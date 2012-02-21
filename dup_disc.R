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
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

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

disc <- duplicateDiscordance(genoData, subjName.col=config["annot_scan_subjectCol"],
                             scan.exclude=scan.exclude, snp.exclude=snp.exclude,
                             corr.by.snp=as.logical(config["corr.by.snp"]))


# by subject
# each entry is a 2x2 matrix, but only one value of each
# is important for pairs
# if there are >2 samples per subject, the first pair is chosen
npair <- length(disc$discordance.by.subject)
disc.subj <- rep(NA, npair); names(disc.subj) <- names(disc$discordance.by.subject)
corr.subj <- rep(NA, npair); names(corr.subj) <- names(disc$correlation.by.subject)
for (i in 1:npair) {
  disc.subj[i] <- disc$discordance.by.subject[[i]][1,2]
  corr.subj[i] <- disc$correlation.by.subject[[i]][1,2]
}
disc$discordance.by.pair <- disc.subj
disc$correlation.by.pair <- corr.subj

# by snp
# probability of discordance for various error rates
(N <- max(disc$discordance.by.snp$npair))
prob.disc <- duplicateDiscordanceProbability(N)

# find out how  many snps fall into each category of discordance
num <- rep(NA, 8)
discordant <- disc$discordance.by.snp$discordant
for(i in 1:8) num[i] <- length(discordant[!is.na(discordant) & discordant>(i-1)])
prob.tbl <- cbind(prob.disc, num)
disc$probability <- prob.tbl

save(disc, file=config["out_disc_file"])
