##########
# PCA plots - combined with external data
# Usage: R --args config.file < pca_combined_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# scan annotation
scanAnnot <- getobj(config["annot_scan_file"])
scan1 <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"]))
names(scan1) <- c("scanID", "race")
if (sum(is.na(scan1$race)) > 0 & hasVariable(scanAnnot, config["ext_annot_scan_raceCol"])) {
  scan1$race2 <- getVariable(scanAnnot, config["ext_annot_scan_raceCol"])
  scan1$race[is.na(scan1$race)] <- scan1$race2[is.na(scan1$race)]
  scan1$race2 <- NULL
}
ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
scan2 <- getVariable(ext.scanAnnot, c("scanID", config["ext_annot_scan_raceCol"]))
names(scan2) <- c("scanID", "race")
samp <- rbind(scan1, scan2)

# get PCA results
pca <- getobj(config["out_comb_pca_file"])
samp <- samp[match(pca$sample.id, samp$scanID),]
stopifnot(allequal(pca$sample.id, samp$scanID))
table(samp$race, exclude=NULL)

# color by race
Sys.setlocale("LC_COLLATE", "C")
race <- sort(unique(samp$race))
stopifnot(all(race %in% names(config)))
for (r in race) {
  sel <- samp$race == r
  samp$plotcol[sel] <- config[r]
}
table(samp$plotcol, exclude=NULL)

# labels
(x <- pca$eigenval[1:4]/sum(pca$eigenval))
lbls <- paste("EV", 1:4, " (", format(100*x,digits=2), "%)", sep="")

# plot the first four PCs
png(config["out_comb_pairs_plot"], width=720, height=720)
par(cex=1.5, lwd=1.5, cex.lab=1.5, cex.axis=1.2)
pairs(pca$eigenvect[,1:4], labels=lbls, col=samp$plotcol)
dev.off()

# plot EV1 vs EV2
# one group at a time so we don't cover up smaller sample sets
pdf(config["out_comb_ev12_plot"], width=6, height=6)
plot(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], type="n")
tbl <- table(samp$race)
raceOrd <- names(tbl)[order(tbl, decreasing=TRUE)]
for (r in raceOrd) {
  sel <- samp$race == r
  points(pca$eigenvect[sel,1], pca$eigenvect[sel,2], col=samp$plotcol[sel])
}
legend("top", legend=race, col=config[race], pch=rep(1, length(race)))
dev.off()


#plot SNP-PC correlation
corr <- getobj(config["out_comb_corr_file"])
snpAnnot <- getobj(config["annot_snp_file"])
snp <- snpAnnot[match(corr$snp.id, snpAnnot$rsID),]
chrom.labels <- unique(getChromosome(snp, char=TRUE))

png(paste(config["out_comb_corr_plot_prefix"], "_%03d.png", sep=""), height=720, width=720)
par(cex=1.5, lwd=1.5, cex.lab=1.5, cex.axis=1.2, cex.main=1.5)
par(mfrow=c(4,1))
for(i in 1:8){
  snpCorrelationPlot(abs(corr$snpcorr[i,]), snp$chromosome,
                     chrom.labels=chrom.labels,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

