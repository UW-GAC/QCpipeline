##########
# PCA plots
# Usage: R --args config.file pca.type < pca_plots.R
##########

library(GWASTools)
library(QCpipeline)
library(MASS)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_scan_raceCol", "annot_snp_file",
              "out_corr_file", "out_pca_file")
optional <- c("annot_scan_ethnCol", "annot_snp_rsIDCol", "ext_annot_scan_file",
              "ext_annot_scan_raceCol",
              "num_evs_to_plot", "out_comb_annot_snp_file",
              "out_corr_plot_prefix", "out_corr_pruned_plot_prefix",
              "out_dens_plot", "out_ev12_plot", "out_pairs_plot", "out_scree_plot",
              "out_parcoord_plot",
              "out_ev12_plot_hapmap", "out_ev12_plot_study")
default <- c(NA, "rsID", NA, "pop.group", 12, NA, "pca_corr", NA, "pca_dens.pdf",
             "pca_ev12.pdf", "pca_pairs.png", "pca_scree.pdf",
             "pca_parcoord.pdf",
             "pca_ev12_hapmap.pdf", "pca_ev12_study.pdf")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# check for type
if (length(args) < 2) stop("missing pca type (study or combined)")
type <- args[2]

# scan annotation
if (type == "study") {
  scanAnnot <- getobj(config["annot_scan_file"])
  samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"]))
  names(samp) <- c("scanID", "race")
  if (!is.na(config["annot_scan_ethnCol"])) {
    samp$ethnicity <- getVariable(scanAnnot, config["annot_scan_ethnCol"])
  } else samp$ethnicity <- NA
} else if (type == "combined") {
  scanAnnot <- getobj(config["annot_scan_file"])
  scan1 <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"]))
  names(scan1) <- c("scanID", "race")
  if (!is.na(config["annot_scan_ethnCol"])) {
    scan1$ethnicity <- getVariable(scanAnnot, config["annot_scan_ethnCol"])
  } else scan1$ethnicity <- NA
  if (sum(is.na(scan1$race)) > 0 & hasVariable(scanAnnot, config["ext_annot_scan_raceCol"])) {
    scan1$race2 <- getVariable(scanAnnot, config["ext_annot_scan_raceCol"])
    scan1$race[is.na(scan1$race)] <- scan1$race2[is.na(scan1$race)]
    scan1$race2 <- NULL
  }
  ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
  scan2 <- getVariable(ext.scanAnnot, c("scanID", config["ext_annot_scan_raceCol"]))
  names(scan2) <- c("scanID", "race")
  scan2$ethnicity <- NA
  samp <- rbind(scan1, scan2)
} else {
  stop("pca type must be study or combined")
}

# get PCA results
pca <- getobj(config["out_pca_file"])
samp <- samp[match(pca$sample.id, samp$scanID),]
stopifnot(allequal(pca$sample.id, samp$scanID))
table(samp$race, samp$ethnicity, useNA="ifany")

# color by race
table(samp$race, useNA="ifany")
Sys.setlocale("LC_COLLATE", "C")
samp$plotcol <- "black"
race <- as.character(sort(unique(samp$race)))
if (length(race) > 0) {
  stopifnot(all(race %in% names(config)))
  for (r in race) {
    sel <- samp$race %in% r
    samp$plotcol[sel] <- config[r]
  }
}
table(samp$plotcol, useNA="ifany")

# plot symbol by ethnicity
table(samp$ethnicity, useNA="ifany")
samp$plotsym <- 1
ethn <- as.character(sort(unique(samp$ethnicity)))
if (length(ethn) > 0) {
  stopifnot(all(ethn %in% names(config)))
  for (e in ethn) {
    sel <- samp$ethn %in% e
    samp$plotsym[sel] <- as.integer(config[e])
  }
}
table(samp$plotsym, useNA="ifany")

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
legend(bestLegendPos(pca$eigenvect[,1], pca$eigenvect[,2]), legend=c(race, ethn),
       col=c(config[race], rep("black", length(ethn))),
       pch=c(rep(1, length(race)), as.integer(config[ethn])))
dev.off()

# parallel coordinates plot
# should eventually choose alpha more intelligently based on sample size.
pdf(config["out_parcoord_plot"], width=8, height=6)
# plot by order of number of samples in each race
samp$nrace <- table(samp$race, useNA="ifany")[samp$race]
samp$nrace[is.na(samp$nrace)] <- sum(is.na(samp$nrace))
i_ord <- order(-samp$nrace)
samp$alpha <- ifelse(samp$nrace < 10, 1,
                     ifelse(samp$nrace < 100, 0.5,
                            ifelse(samp$nrace < 1000, 0.3, 0.1))) * 255
# adjust margins so legend can be on the left
par(mar=c(5.1, 10.1, 3.1, 2.1), xpd=TRUE)
# need to subset both the eigenvector matrix and the colors here
parcoord(pca$eigenvect[i_ord, 1:8], bty="L",
  col=rgb(t(col2rgb(samp$plotcol)), alpha=samp$alpha, maxColorValue=255)[i_ord])
title(xlab="Eigenvector")
if(any(is.na(samp$race))) {
  legendNames <- c(race, "Unknown")
  legendCols <- c(config[race], "black")
} else {
  legendNames <- race
  legendCols <- config[race]
}
legend("left", inset=c(-0.32,0), legend=legendNames, col=legendCols, lty=1, bty="n")
dev.off()
## to do: figure out how to plot two different parcoord plots on the same range for combined.


if (type == "combined"){
	# plot hapmaps separately from study subjects
	pdf(config["out_ev12_plot_hapmap"], width=6, height=6)
	plot(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], type="n")
	ext.sel <- (samp$scanID %in% ext.scanAnnot$scanID) | (scanAnnot$geno.cntl[match(samp$scanID, scanAnnot$scanID)] %in% 1) # get external or study hapmaps
	tbl <- table(samp$plotcol[ext.sel])
	colOrd <- names(tbl)[order(tbl, decreasing=TRUE)]
	ext.race <- as.character(sort(unique(samp$race[ext.sel])))
	ext.ethn <- as.character(sort(unique(samp$ethn[ext.sel])))
	for (r in colOrd) {
	  sel <- (samp$plotcol == r) & (ext.sel)
	  points(pca$eigenvect[sel,1], pca$eigenvect[sel,2], col=samp$plotcol[sel], pch=samp$plotsym[sel])
	}
	legend(bestLegendPos(pca$eigenvect[,1], pca$eigenvect[,2]), legend=c(ext.race, ext.ethn),
		   col=c(config[ext.race], rep("black", length(ext.ethn))),
		   pch=c(rep(1, length(ext.race)), as.integer(config[ext.ethn])))
	dev.off()

	pdf(config["out_ev12_plot_study"], width=6, height=6)
	plot(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], type="n")
	ext.sel <- (samp$scanID %in% scanAnnot$scanID) & (scanAnnot$geno.cntl[match(samp$scanID, scanAnnot$scanID)] %in% 0) # get only study subjects (not study hapmaps)
	tbl <- table(samp$plotcol[ext.sel])
	colOrd <- names(tbl)[order(tbl, decreasing=TRUE)]
	ext.race <- as.character(sort(unique(samp$race[ext.sel])))
	ext.ethn <- as.character(sort(unique(samp$ethn[ext.sel])))
	for (r in colOrd) {
	  sel <- (samp$plotcol == r) & (ext.sel)
	  points(pca$eigenvect[sel,1], pca$eigenvect[sel,2], col=samp$plotcol[sel], pch=samp$plotsym[sel])
	}
	legend(bestLegendPos(pca$eigenvect[,1], pca$eigenvect[,2]), legend=c(ext.race, ext.ethn),
		   col=c(config[ext.race], rep("black", length(ext.ethn))),
		   pch=c(rep(1, length(ext.race)), as.integer(config[ext.ethn])))
	dev.off()
}

# plot density on sides
pdf(config["out_dens_plot"], width=7, height=7)
plot2DwithHist(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], col=samp$plotcol, pch=samp$plotsym)
dev.off()

#plot SNP-PC correlation
if (type == "combined")  {
  snpAnnot <- getobj(config["out_comb_annot_snp_file"])
} else {
  snpAnnot <- getobj(config["annot_snp_file"])
}
corr <- getobj(config["out_corr_file"])
snp <- snpAnnot[match(corr$snp.id, getSnpID(snpAnnot)),]
chrom <- getChromosome(snp, char=TRUE)

nev <- as.integer(config["num_evs_to_plot"])

png(paste(config["out_corr_plot_prefix"], "_%03d.png", sep=""), height=720, width=720)
par(mfrow=c(4,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
for(i in 1:nev){
  snpCorrelationPlot(abs(corr$snpcorr[i,]), chrom,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

if (!is.na(config["out_corr_pruned_plot_prefix"])) {
  snps.pruned <- getobj(config["out_pruned_file"])
  ind <- getSnpID(snp) %in% snps.pruned
  
  png(paste(config["out_corr_pruned_plot_prefix"], "_%03d.png", sep=""), height=720, width=720)
  par(mfrow=c(4,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
  for(i in 1:nev){
    snpCorrelationPlot(abs(corr$snpcorr[i,ind]), chrom[ind],
                       main=paste("Eigenvector",i), ylim=c(0,1))
  }
  dev.off()
}
  
# scree plot
x <- pca$eigenval/sum(pca$eigenval)
pdf(config["out_scree_plot"], width=6, height=6)
plot(1:nev,100*x[1:nev], xlab="Eigenvector", ylab="Percent of variance accounted for")
dev.off()
