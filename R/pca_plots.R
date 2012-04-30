##########
# PCA plots
# Usage: R --args config.file pca.type < pca_plots.R
##########

library(GWASTools)
library(QCpipeline)
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
              "num_evs_to_plot", "out_corr_plot_prefix", "out_corr_pruned_plot_prefix",
              "out_dens_plot", "out_ev12_plot", "out_pairs_plot", "out_scree_plot")
default <- c(NA, "rsID", NA, "pop.group", 12, "pca_corr", NA, "pca_dens.pdf",
             "pca_ev12.pdf", "pca_pairs.png", "pca_scree.pdf")
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

# plot density on sides
pdf(config["out_dens_plot"], width=7, height=7)
plot2DwithHist(pca$eigenvect[,1], pca$eigenvect[,2], xlab=lbls[1], ylab=lbls[2], col=samp$plotcol, pch=samp$plotsym)
dev.off()

#plot SNP-PC correlation
if (type == "combined")  {
  idCol <- config["annot_snp_rsIDCol"]
} else {
  idCol <- "snpID"
}
corr <- getobj(config["out_corr_file"])
snpAnnot <- getobj(config["annot_snp_file"])
snp <- snpAnnot[match(corr$snp.id, getVariable(snpAnnot, idCol)),]
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