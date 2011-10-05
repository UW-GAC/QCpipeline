##########
# Duplicate discordance
# Usage: R --args config.file < dup_disc.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]
print(config.table)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

# exclude missing SNPs
snp.exclude <- snpID[getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1]
length(snp.exclude)

disc <- duplicateDiscordance(genoData, subjName.col=config["annot_scan_subjectCol"],
                             scan.exclude=scan.exclude, snp.exclude=snp.exclude)


# by subject
# each entry is a 2x2 matrix, but only one value of each
# is important for pairs
# we are assuming no more than 2 samples per subject
npair <- length(disc$discordance.by.subject)
disc.subj <- rep(NA, npair); names(disc.subj) <- names(disc$discordance.by.subject)
corr.subj <- rep(NA, npair); names(corr.subj) <- names(disc$correlation.by.subject)
for (i in 1:npair) {
  disc.subj[i] <- disc$discordance.by.subject[[i]][1,2]
  corr.subj[i] <- disc$correlation.by.subject[[i]][1,2]
}
disc$discordance.by.pair <- disc.subj
disc$correlation.by.pair <- corr.subj

# plot
disc.subj <- disc.subj[order(disc.subj)]
rank <- 1:npair
pdf(config["out_disc_plot"], width=6, height=6)
plot(disc.subj, rank, xlab="discordance rate", ylab="rank",
     main=paste("Discordance in", npair, "duplicate sample pairs"))
dev.off()


# by snp
# probability of discordance for various error rates
N <- max(disc$discordance.by.snp$npair)
prob.disc <- duplicateDiscordanceProbability(N)

# find out how  many snps fall into each category of discordance
num <- rep(NA, 8)
discordant <- disc$discordance.by.snp$discordant
for(i in 1:8) num[i] <- length(discordant[!is.na(discordant) & discordant>(i-1)])
prob.tbl <- cbind(prob.disc, num)
disc$probability <- prob.tbl

save(disc, file=config["out_disc_file"])
