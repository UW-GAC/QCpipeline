##########
# SNP filter plots (duplicate discordance, MAF histograms)
# Usage: R --args config.file < snp_filt_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file")
optional <- c("annot_scan_hapmapCol", "annot_scan_subjectCol", "corr.by.snp", "maf.bin", "out_afreq_file", "out_disc_file", "out_disc_maf_file", "out_disc_plot", "out_ma_conc_plot", "out_maf_autosomes_plot", "out_maf_plot", "out_maf_xchrom_plot", "out_snp_conc_plot", "out_snp_corr_plot")
default <- c("geno.cntl", "subjectID", FALSE, 0.01, "allele_freq.RData", "dup_disc.RData", "dup_disc_maf.RData", "dup_disc.pdf", "snp_ma_conc.pdf", "maf_aut_hist.pdf", "maf_hist.pdf", "maf_x_hist.pdf", "snp_conc.pdf", "snp_corr.pdf")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))

(snpAnnot <- getobj(config["annot_snp_file"]))

# duplicate discordance by subject
disc <- getobj(config["out_disc_file"])
disc.subj <- disc$discordance.by.pair
disc.subj <- disc.subj[order(disc.subj)]
npair <- length(disc.subj)
rank <- 1:npair
# color-code by hapmap
hapmap <- getVariable(scanAnnot, config["annot_scan_hapmapCol"])
subj <- getVariable(scanAnnot, config["annot_scan_subjectCol"])
hapmap.ids <- unique(subj[hapmap %in% 1])
plotcol <- rep("black", length(disc.subj))
plotcol[names(disc.subj) %in% hapmap.ids] <- "red"
pdf(config["out_disc_plot"], width=6, height=6)
plot(disc.subj, rank, xlab="discordance rate", ylab="rank", col=plotcol,
     main=paste("Discordance in", npair, "duplicate sample pairs"))
if (sum(plotcol == "red") > 0) {
    legend(bestLegendPos(disc.subj, rank), c("study", "HapMap"), col=c("black", "red"), pch=c(1,1))
}
dev.off()

# summary of study data only
summary(disc.subj[!(names(disc.subj) %in% hapmap.ids)])


# snp plots
# MAF bins
afreq <- getobj(config["out_afreq_file"])
maf <- pmin(afreq[,"all"], 1-afreq[,"all"])
maf.bin <- as.numeric(config["maf.bin"])
bins <- seq(0, 0.5, maf.bin)
refmaf <- bins[2:length(bins)] - maf.bin/2
mafbin <- rep(NA, length(maf))
for (i  in 1:length(bins)-1) {
  thisbin <- bins[i] < maf & maf <= bins[i+1]
  mafbin[thisbin] <- paste(bins[i], "-", bins[i+1])
}
table(mafbin)
pdf(config["out_maf_plot"], width=6, height=6)
hist(maf, breaks=bins, xlab="MAF", main="")
dev.off()
# autosomes
chrom <- getChromosome(snpAnnot)
pdf(config["out_maf_autosomes_plot"], width=6, height=6)
hist(maf[chrom < 23], breaks=bins, xlab="MAF", main="Autosomes")
dev.off()
# x chrom
pdf(config["out_maf_xchrom_plot"], width=6, height=6)
hist(maf[chrom == 23], breaks=bins, xlab="MAF", main="X chromosome")
dev.off()


# concordance
maf <- maf[rownames(afreq) %in% disc$discordance.by.snp$snpID]
snp.conc <- 1 - disc$discordance.by.snp$discord.rate
meanconc <- rep(NA, length(refmaf))
for (i  in 1:length(bins)-1) {
  thisbin <- bins[i] < maf & maf <= bins[i+1]
  meanconc[i] <- mean(snp.conc[thisbin], na.rm=TRUE)
}
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

# minor allele concordance
if (!is.na(config["out_disc_maf_file"]) & file.exists(config["out_disc_maf_file"])) {
  disc <- getobj(config["out_disc_maf_file"])
  snp.conc <- 1 - disc$discordance.by.snp$discord.rate
  meanconc <- rep(NA, length(refmaf))
  for (i  in 1:length(bins)-1) {
    thisbin <- bins[i] < maf & maf <= bins[i+1]
    meanconc[i] <- mean(snp.conc[thisbin], na.rm=TRUE)
  }
  table(mafbin)
  pdf(config["out_ma_conc_plot"], width=6, height=6)
  plot(refmaf, meanconc, xlab="MAF", ylab="minor allele concordance", xlim=c(0,0.5))
  dev.off()
}
