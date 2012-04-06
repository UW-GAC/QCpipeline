##########
# Anomaly stats and plots
# Usage: R --args config.file maf < anom_stats.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# read MAF threshold
if (length(args) > 1) {
  maf <- as.numeric(args[2])
} else {
  maf <- NULL
}
print(maf)

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
chrom <- getChromosome(snpAnnot, char=TRUE)
pos <- getPosition(snpAnnot)
hla.df <- get(data(list=paste("HLA", config["build"], sep=".")))
hla <- chrom == "6" & pos >= hla.df$start.base & pos <= hla.df$end.base
xtr.df <- get(data(list=paste("pseudoautosomal", config["build"], sep=".")))
xtr <- chrom == "X" & pos >= xtr.df["X.XTR", "start.base"] & pos <= xtr.df["X.XTR", "end.base"]
centromeres <- get(data(list=paste("centromeres", config["build"], sep=".")))
gap <- rep(FALSE, length(snpID))
for (i in 1:nrow(centromeres)) {
  ingap <- chrom == centromeres$chrom[i] & pos > centromeres$left.base[i] &
    pos < centromeres$right.base[i]
  gap <- gap | ingap
}
table(chrom, gap)

#ignore includes intensity-only and failed snps
ignore <- getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1 
snp.exclude <- ignore | hla | xtr | gap
table(snp.exclude)

# maf threshold
if (!is.null(maf)) {
  afreq <- getobj(config["out_afreq_file"])
  stopifnot(allequal(rownames(afreq), snpID))
  maf.filt <- is.na(afreq[,"all"]) | afreq[,"all"] <= maf | afreq[,"all"] >= (1-maf)
  snp.exclude <- snp.exclude | maf.filt
}
table(snp.exclude)

snp.ok <- snpID[!snp.exclude]
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
if (is.na(config["thresh.sum"])) config["thresh.sum"] <- 10
if (is.na(config["thresh.indiv"])) config["thresh.indiv"] <- 5
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
  snp.ineligible <- snpID[snp.exclude]
  png.file <- file.path(config["out_plot_dir"], paste(config["out_plot_prefix"], "_%003d.png", sep=""))
  png(png.file, width=720, height=720)
  anomStatsPlotIdeogram(blData, genoData, anom.stats=long.chk, snp.ineligible=snp.ineligible,
                win=as.integer(config["out_plot_win"]), centromere=centromeres, cex=0.25)
  dev.off()
} else {
  message("no long anomalies found")
}
