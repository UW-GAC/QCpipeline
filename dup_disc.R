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

# plot
disc.subj <- disc.subj[order(disc.subj)]
rank <- 1:npair
# color-code by hapmap
hapmap <- getVariable(scanAnnot, config["annot_scan_hapmapCol"])
subj <- getVariable(scanAnnot, config["annot_scan_subjectCol"])
hapmap.ids <- subj[hapmap %in% 1]
plotcol <- rep("black", length(disc.subj))
plotcol[names(disc.subj) %in% hapmap.ids] <- "red"
pdf(config["out_disc_plot"], width=6, height=6)
plot(disc.subj, rank, xlab="discordance rate", ylab="rank", col=plotcol,
     main=paste("Discordance in", npair, "duplicate sample pairs"))
legend(bestLegendPos(disc.subj, rank), c("study", "HapMap"), col=c("black", "red"), pch=c(1,1))
dev.off()


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

# snp plots
snp.conc <- 1 - disc$discordance.by.snp$discord.rate

# MAF bins
afreq <- getobj(config["out_afreq_file"])
maf <- pmin(afreq[,"all"], 1-afreq[,"all"])
maf <- maf[snpID %in% disc$discordance.by.snp$snpID]
bins <- seq(0, 0.5, 0.05)
refmaf <- bins[2:length(bins)] - 0.025
mafbin <- rep(NA, length(maf))
meanconc <- rep(NA, length(refmaf))
if (as.logical(config["corr.by.snp"])) meancorr <- rep(NA, length(refmaf))
for (i  in 1:length(bins)-1) {
  thisbin <- bins[i] < maf & maf <= bins[i+1]
  mafbin[thisbin] <- paste(bins[i], "-", bins[i+1])
  meanconc[i] <- mean(snp.conc[thisbin], na.rm=TRUE)
}
table(mafbin)

# concordance
pdf(config["out_snp_conc_plot"], width=6, height=6)
plot(refmaf, meanconc, xlab="MAF", ylab="concordance", xlim=c(0,0.5))
dev.off()

# correlation
if (as.logical(config["corr.by.snp"])) {
  snp.corr <- disc$discordance.by.snp$correlation
  meancorr <- rep(NA, length(refmaf))
  for (i  in 1:length(bins)-1) {
    thisbin <- bins[i] < maf & maf <= bins[i+1]
    meancorr[i] <- mean(snp.corr[thisbin], na.rm=TRUE)
  }
  pdf(config["out_snp_corr_plot"], width=6, height=6)
  plot(refmaf, meancorr, xlab="MAF", ylab="correlation of allelic dosage", xlim=c(0,0.5))
  dev.off()
}
