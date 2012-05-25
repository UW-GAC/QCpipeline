##########
# Anomaly stats and plots
# Usage: R --args config.file < anom_stats.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "build", "nc_bl_file", "nc_geno_file",
              "out_anom_dir", "out_eligible_snps", "out_plot_dir", "project")
optional <- c("annot_snp_IntensityOnlyCol", "chromXY", "out_plot_prefix", "plot.win",
              "thresh.indiv", "thresh.sum", "scan_exclude_file")
default <- c(NA, FALSE, "long_plot", 1, 5, 10, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# create GenotypeData and IntensityData objects
(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

bl.file <- config["nc_bl_file"]
blnc <- NcdfIntensityReader(bl.file)
blData <-  IntensityData(blnc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

geno.file <- config["nc_geno_file"]
genonc <- NcdfGenotypeReader(geno.file)
genoData <-  GenotypeData(genonc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# select SNPs
snp.ok <- getobj(config["out_eligible_snps"])
length(snp.ok)

file <- file.path(config["out_anom_dir"], paste(config["project"], "BAF.filtered.all.RData", sep="."))
BAF <- getobj(file); nrow(BAF)
if (!is.null(BAF)) if (nrow(BAF) == 0) BAF <- NULL

file <- file.path(config["out_anom_dir"], paste(config["project"], "LOH.filtered.all.RData", sep="."))
LOH <- getobj(file); nrow(LOH)
if (!is.null(LOH)) if (nrow(LOH) == 0) LOH <- NULL

# combine BAF and LOH
if (!is.null(BAF)) {
  BAF$method <- "BAF"
  BAF$mad.fac <- NA
}
if (!is.null(LOH)) {
  LOH$method <- "LOH"
}
if (!is.null(BAF) & !is.null(LOH)) {
  cols <- intersect(names(BAF), names(LOH))
  anoms <- rbind(BAF[,cols], LOH[,cols])
} else if (!is.null(BAF)) {
  anoms <- BAF
} else if (!is.null(LOH)) {
  anoms <- LOH
} else {
  stop("no anomalies detected")
}
anoms <- anoms[order(anoms$scanID, anoms$chromosome, anoms$left.index),]; dim(anoms)
anoms$anom.id <- 1:nrow(anoms)

# do XY separately
if (as.logical(config["chromXY"])) {
  anoms.XY <- anoms[anoms$chromosome == XYchromCode(blData),]
  if (nrow(anoms.XY) > 0) { 
    anoms <- anoms[anoms$chromosome != XYchromCode(blData),]
  }
}

centromeres <- get(data(list=paste("centromeres", config["build"], sep=".")))
stats <- anomSegStats(blData, genoData, snp.ids=snp.ok, anom=anoms, centromere=centromeres)
saveas(stats, paste(config["project"],".anom.stats",sep=""), config["out_anom_dir"])

if (as.logical(config["chromXY"])) {
  if (nrow(anoms.XY) > 0) { 
    if (!is.na(config["annot_snp_IntensityOnlyCol"])) {
      io <- getVariable(snpAnnot, config["annot_snp_IntensityOnlyCol"])
    } else {
      io <- rep(0, length(snpID))
    }
    sXY <- is.element(chrom, "XY") & is.element(io, 0)
    sel.XY  <-  snpID[sXY]
    XY <- centromeres[is.element(centromeres$chrom,"X"),]
    XY$chrom <- "XY"
    rn <- row.names(centromeres)
    centromeres <- rbind(centromeres,XY)
    row.names(centromeres) <- c(rn,"XY")
    stats.XY <- anomSegStats(blData, genoData, snp.ids=sel.XY, anom=anoms.XY, centromere=centromeres)
    saveas(stats.XY, paste(config["project"],".anom.stats.chXY",sep=""), config["out_anom_dir"])
  }
}

# select long anomalies
mb <- 1000000
# all anomalies for sample-chromosomes with sum > thresh.sum
long.thresh <- as.numeric(config["thresh.sum"]) * mb
psan <- paste(stats$scanID,stats$chromosome)
usc <- unique(psan)
sc.chk <- NULL
for(sc in usc){
  an <- stats$nbase[psan==sc]
  lens <- sum(an)
  if(lens >= long.thresh) sc.chk <- c(sc.chk,sc)
}
# any individual anomalies > thresh.indiv
indiv.chk <- stats$nbase >= as.numeric(config["thresh.indiv"]) * mb
sum.chk <- is.element(psan,sc.chk)
table(sum.chk, indiv.chk)
any.chk <- sum.chk | indiv.chk

if (sum(any.chk) > 0) {
  long.chk <- stats[any.chk,]
  long.chk <- long.chk[order(long.chk$scanID, long.chk$chromosome, long.chk$left.base),]

  # put one anomaly per chromosome first
  psl <- paste(long.chk$scanID, long.chk$chromosome)
  w <- duplicated(psl)
  lc1 <- long.chk[!is.element(psl,psl[w]),]
  lc2 <- long.chk[is.element(psl,psl[w]),]
  long.chk <- rbind(lc1, lc2)
  n <- nrow(long.chk)
  
  datx <- data.frame("scanID"=long.chk$scanID,"chromosome"=long.chk$chromosome, "anom.id"=long.chk$anom.id,"plot.num"=1:n, filter="", "comment"="", stringsAsFactors=FALSE)
  csv.file <- file.path(config["out_plot_dir"], paste(config["out_plot_prefix"], ".csv", sep=""))
  write.csv(datx, file=csv.file, quote=FALSE, row.names=FALSE)

  # plot (non-XY)
  snp.ineligible <- setdiff(snpID, snp.ok)
  png.file <- file.path(config["out_plot_dir"], paste(config["out_plot_prefix"], "_%003d.png", sep=""))
  png(png.file, width=720, height=720)
  anomStatsPlot(blData, genoData, anom.stats=long.chk, snp.ineligible=snp.ineligible,
                win=as.integer(config["plot.win"]), centromere=centromeres, cex=0.25,
                cex.axis=1.7, cex.main=1.7, cex.lab=1.7)
  dev.off()
} else {
  message("no long anomalies found")
}
