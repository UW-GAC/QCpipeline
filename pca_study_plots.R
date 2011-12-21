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
samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"]))
names(samp) <- c("scanID", "race")

# get PCA results
pca <- getobj(config["out_pca_file"])
samp <- samp[match(pca$sample.id, samp$scanID),]
stopifnot(allequal(pca$sample.id, samp$scanID))
table(samp$race, exclude=NULL)

# color by race
Sys.setlocale("LC_COLLATE", "C")
race <- as.character(sort(unique(samp$race)))
stopifnot(all(race %in% names(config)))
samp$plotcol <- "black"
for (r in race) {
  sel <- samp$race %in% r
  samp$plotcol[sel] <- config[r]
}
table(samp$plotcol, exclude=NULL)

# plot race==NA with different character
samp$plotsym <- 1
samp$plotsym[is.na(samp$race)] <- 4
table(samp$plotsym, exclude=NULL)

# labels
(x <- pca$eigenval[1:4]/sum(pca$eigenval))
lbls <- paste("EV", 1:4, " (", format(100*x,digits=2), "%)", sep="")

# plot the first four PCs
png(config["out_pairs_plot"], width=720, height=720)
par(lwd=1.5, cex.axis=1.5)
pairs(pca$eigenvect[,1:4], labels=lbls, col=samp$plotcol, pch=samp$plotsym)
dev.off()

# plot EV1 vs EV2
# one group at a time so we don't cover up smaller sample sets
pdf(config["out_ev12_plot"], width=6, height=6)
plot(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], type="n")
tbl <- table(samp$plotcol)
colOrd <- names(tbl)[order(tbl, decreasing=TRUE)]
for (r in colOrd) {
  sel <- samp$plotcol == r
  points(pca$eigenvect[sel,1], pca$eigenvect[sel,2], col=samp$plotcol[sel], pch=samp$plotsym[sel])
}
if (sum(is.na(samp$race)) > 0) {
  legend("top", legend=c(race, "NA"), col=c(config[race], "black"), pch=c(rep(1, length(race)), 4))
} else {
  legend("top", legend=race, col=config[race], pch=rep(1, length(race)))
}
dev.off()

# plot density on sides
pdf(config["out_dens_plot"], width=7, height=7)
plot2DwithHist(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], col=samp$plotcol, pch=samp$plotsym)
dev.off()

#plot SNP-PC correlation
corr <- getobj(config["out_corr_file"])
snpAnnot <- getobj(config["annot_snp_file"])
snp <- snpAnnot[match(corr$snp.id, snpAnnot$snpID),]
chrom.labels <- unique(getChromosome(snp, char=TRUE))

png(paste(config["out_corr_plot_prefix"], "_%03d.png", sep=""), height=720, width=720)
par(mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
par(mfrow=c(4,1))
for(i in 1:8){
  snpCorrelationPlot(abs(corr$snpcorr[i,]), snp$chromosome,
                     chrom.labels=chrom.labels,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

