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
required <- c("annot_scan_file", "annot_snp_file", "nc_geno_file")
optional <- c("annot_scan_subjectCol", "annot_snp_missingCol", "corr.by.snp",
              "out_disc_file", "disc_scan_exclude_file")
default <- c("subjectID", "missing.n1", FALSE, "dup_disc.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
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
# give output data frame the same dimensions as snp annotation
snp <- merge(data.frame(snpID), disc$discordance.by.snp, all.x=TRUE)
disc$discordance.by.snp <- snp

# probability of discordance for various error rates
(N <- max(disc$discordance.by.snp$npair, na.rm=TRUE))
max.disc <- 7
prob.disc <- duplicateDiscordanceProbability(N, max.disc=max.disc)

# find out how  many snps fall into each category of discordance
ncat <- max.disc + 1
num <- rep(NA, ncat)
discordant <- disc$discordance.by.snp$discordant
for(i in 1:ncat) num[i] <- length(discordant[!is.na(discordant) & discordant>(i-1)])
prob.tbl <- cbind(prob.disc, num)
disc$probability <- prob.tbl

save(disc, file=config["out_disc_file"])
