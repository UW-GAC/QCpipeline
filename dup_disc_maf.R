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

afreq <- getobj(config["out_afreq_file"])
disc <- duplicateDiscordance(genoData, subjName.col=config["annot_scan_subjectCol"],
                             scan.exclude=scan.exclude, snp.exclude=snp.exclude,
                             minor.allele.only=TRUE, allele.freq=afreq[,"all"])

save(disc, file=config["out_disc_maf_file"])

snp.conc <- 1 - disc$discordance.by.snp$discord.rate

# MAF bins
maf <- pmin(afreq[,"all"], 1-afreq[,"all"])
maf <- maf[snpID %in% disc$discordance.by.snp$snpID]
bins <- seq(0, 0.5, 0.05)
refmaf <- bins[2:length(bins)] - 0.025
mafbin <- rep(NA, length(maf))
meanconc <- rep(NA, length(refmaf))
for (i  in 1:length(bins)-1) {
  thisbin <- bins[i] < maf & maf <= bins[i+1]
  mafbin[thisbin] <- paste(bins[i], "-", bins[i+1])
  meanconc[i] <- mean(snp.conc[thisbin], na.rm=TRUE)
}
table(mafbin)

pdf(config["out_ma_conc_plot"], width=6, height=6)
plot(refmaf, meanconc, xlab="MAF", ylab="minor allele concordance", xlim=c(0,0.5))
dev.off()
