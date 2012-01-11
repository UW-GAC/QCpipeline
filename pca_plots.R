##########
# PCA plots - combined with external data
# Usage: R --args config.file pca.type < pca_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# check for type
if (length(args) < 2) stop("missing pca type (study or combined)")
type <- args[2]

# scan annotation
if (type == "study") {
  scanAnnot <- getobj(config["annot_scan_file"])
  samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"]))
  names(samp) <- c("scanID", "race")
} else if (type == "combined") {
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
} else {
  stop("pca type must be study or combined")
}

# get PCA results
if (type == "combined") { pcafile <- config["out_comb_pca_file"]  
} else pcafile <- config["out_pca_file"]
pca <- getobj(pcafile)
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
if (type == "combined") { plotfile <- config["out_comb_pairs_plot"]
} else plotfile <- config["out_pairs_plot"]
png(plotfile, width=720, height=720)
par(lwd=1.5, cex.axis=1.5)
pairs(pca$eigenvect[,1:4], labels=lbls, col=samp$plotcol, pch=samp$plotsym)
dev.off()

# plot EV1 vs EV2
# one group at a time so we don't cover up smaller sample sets
if (type == "combined") { plotfile <- config["out_comb_ev12_plot"]
} else plotfile <- config["out_ev12_plot"]
pdf(plotfile, width=6, height=6)
plot(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], type="n")
tbl <- table(samp$plotcol)
colOrd <- names(tbl)[order(tbl, decreasing=TRUE)]
for (r in colOrd) {
  sel <- samp$plotcol == r
  points(pca$eigenvect[sel,1], pca$eigenvect[sel,2], col=samp$plotcol[sel], pch=samp$plotsym[sel])
}
if (sum(is.na(samp$race)) > 0) {
  legend(bestLegendPos(pca$eigenvect[,1], pca$eigenvect[,2]), legend=c(race, "NA"), col=c(config[race], "black"), pch=c(rep(1, length(race)), 4))
} else {
  legend(bestLegendPos(pca$eigenvect[,1], pca$eigenvect[,2]), legend=race, col=config[race], pch=rep(1, length(race)))
}
dev.off()

# plot density on sides
if (type == "study")  {
  pdf(config["out_dens_plot"], width=7, height=7)
  plot2DwithHist(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], col=samp$plotcol, pch=samp$plotsym)
  dev.off()
}

#plot SNP-PC correlation
if (type == "combined")  {
  corrfile <- config["out_comb_corr_file"]
  plotfile <- config["out_comb_corr_plot_prefix"]
  idCol <- config["annot_snp_rsIDCol"]
  
} else {
  corrfile <- config["out_corr_file"]
  plotfile <- config["out_corr_plot_prefix"]
  idCol <- "snpID"
}
corr <- getobj(corrfile)
snpAnnot <- getobj(config["annot_snp_file"])
snp <- snpAnnot[match(corr$snp.id, getVariable(snpAnnot, idCol)),]
chrom.labels <- unique(getChromosome(snp, char=TRUE))

png(paste(plotfile, "_%03d.png", sep=""), height=720, width=720)
par(mfrow=c(4,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
for(i in 1:8){
  snpCorrelationPlot(abs(corr$snpcorr[i,]), snp$chromosome,
                     chrom.labels=chrom.labels,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

