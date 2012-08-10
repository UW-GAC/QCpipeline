##########
# Sample quality checks
# Usage: R --args config.file < sample_quality.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_scan_raceCol", "annot_snp_file", "baf_mean_file", "baf_sd_file", "nc_bl_file")
optional <- c("annot_scan_hetACol", "annot_snp_missingCol", "ibd_con_file", "out_baf_asym", "out_baf_asym_plot_prefix", "out_baf_sd_boxplot", "out_baf_sd_outliers", "out_baf_sd_plot_prefix", "out_het_boxplot", "out_het_outliers", "out_het_plot_prefix", "out_ibd_plot_prefix", "plot_all_unknown", "race_unknown", "range_het", "range_sd")
default <- c("het.A", "missing.n1", NA, "baf_aysm.RData", "baf_asym", "baf_sd.pdf", "baf_sd_outl.RData", "baf_sd_outl", "het_outl.pdf", "het_outl.RData", "het_outl", "high_con", TRUE, NA, 1.5, 1.5)
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

# exclude missing SNPs from plots
(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)
snp.exclude <- snpID[getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1]

# IntensityData
blnc <- NcdfIntensityReader(config["nc_bl_file"])
blData <- IntensityData(blnc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)


# SD of BAF
baf.sd <- getobj(config["baf_sd_file"])

# definition of outliers - range.sd*i(nterquartile range)
range.sd <- as.numeric(config["range_sd"])

# create a matrix
bsd <- matrix(nrow=nrow(baf.sd[[1]]), ncol=length(baf.sd),
              dimnames=list(rownames(baf.sd[[1]]), names(baf.sd)))
for (i in 1:ncol(bsd)) {
  bsd[,i] <- baf.sd[[i]][,1]
}
pdf(config["out_baf_sd_boxplot"], width=6, height=6)
bp <- boxplot(bsd[,colnames(bsd) %in% 1:22], xlab="chromosome", ylab="SD of BAF",
              las=2, range=range.sd)
dev.off()

# get boxplot hi sd outliers by chromosome
out <- list()
for (i in names(baf.sd)[names(baf.sd) %in% 1:22]) {
  x <- baf.sd[[i]]
  med <- median(x)
  tmp <- boxplot(x, plot=FALSE, range=range.sd)$out
  tmp <- tmp[tmp > med]
  mino <- min(tmp)
  out[[i]] <- as.integer(row.names(x)[x >= mino])
}
range(unlist(lapply(out,length)))

# get the number of times each id appears (i.e. how many chromosomes for which it is an outlier)
y <- table(unlist(out))
table(y)
baf.outl <- as.integer(names(y)[y > 1])
(nbaf <- length(baf.outl))

# review BAF/LRR of chr 1 for all samples for which >1 chromosome is an outlier to look for contamination or other artifact
# we don't necessarily need to look at the chromosome(s) flagged since we are looking for a sample-wide artifact
# use plot without any color-coding

# if we don't have chromosome 1, just take the first chrom
if ("1" %in% names(baf.sd)) {
  chr.sel <- 1
} else {
  chr.sel <- names(baf.sd)[1]
}

if (nbaf > 0) {
  save(baf.outl, file=config["out_baf_sd_outliers"])

  png(paste(config["out_baf_sd_plot_prefix"], "_%03d.png", sep=""), width=720, height=720)
  chromIntensityPlot(blData, scan.ids=baf.outl, chrom.ids=rep(chr.sel, nbaf),
                     info=rep("BAF SD", nbaf), snp.exclude=snp.exclude, cex=0.25,
                     ideogram=FALSE)
  dev.off()
}


# asymmetry

# mean of BAF
baf.mean <- getobj(config["baf_mean_file"])

bmn <- matrix(nrow=nrow(baf.mean[[1]]), ncol=length(baf.mean),
              dimnames=list(rownames(baf.mean[[1]]), names(baf.mean)))
for (i in 1:ncol(bmn)) {
  bmn[,i] <- baf.mean[[i]][,1]
}
range(bmn[,colnames(bmn) %in% 1:22])

# get all the samples with asymmetry > 0.05
asym <- list()
for(i in names(baf.mean)[names(baf.mean) %in% 1:22]){
  x <- baf.mean[[i]]
  tmp <- row.names(x)[x>0.55 | x <0.45]
  asym[[i]] <- tmp
}
y <- table(unlist(asym))
table(y)
asym <- as.integer(names(y)[y > 1])
(nasym <- length(asym))

if (nasym > 0) {
  save(asym, file=config["out_baf_asym"])

  png(paste(config["out_baf_asym_plot_prefix"], "_%03d.png", sep=""), width=720, height=720)
  chromIntensityPlot(blData, scan.ids=asym, chrom.ids=rep(chr.sel, nasym),
                     info=rep("Asymmetry", nasym), snp.exclude=snp.exclude, cex=0.25)
  dev.off()
}


# Find boxplot outliers for heterozygosity, both hi and lo
# heterozygosity by race
samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"], config["annot_scan_hetACol"]))
names(samp) <- c("scanID", "race", "het.A")
samp <- samp[!is.na(samp$race),]

# include unknowns?
if (!is.na(config["race_unknown"]) & as.logical(config["plot_all_unknown"])) {
  unknown.ids <- samp$scanID[samp$race %in% config["race_unknown"]]
  samp <- samp[!(samp$race %in% config["race_unknown"]),]
  (nunk <- length(unknown.ids)) 
  if (nunk > 0) {   
    png(paste(config["out_het_plot_prefix"], "_raceunknown_%03d.png", sep=""), width=720, height=720)
    chromIntensityPlot(blData, scan.ids=unknown.ids, chrom.ids=rep(chr.sel, nunk),
                       info=rep("race unknown", nunk), snp.exclude=snp.exclude, cex=0.25)
    dev.off()
  }
}

# definition of outliers 
range.het <- as.numeric(config["range_het"])

pdf(config["out_het_boxplot"], width=6, height=6)
bp <- boxplot(samp$het.A ~ as.factor(samp$race), range=range.het, varwidth=TRUE, las=2,
              main=config["annot_scan_raceCol"], ylab="Autosomal heterozygosity")
dev.off()

outliers <- samp[samp$het.A %in% bp$out,]
(nhet <- nrow(outliers))
table(outliers$race)

if (nhet > 0) {
  # high or low?
  outliers$code <- NA
  races <- unique(samp$race)
  for (race in races) {
    med <- median(samp$het.A[samp$race == race], na.rm=TRUE)
    outliers$code[outliers$race == race & outliers$het.A > med] <- "high het.A"
    outliers$code[outliers$race == race & outliers$het.A < med] <- "low het.A"
  }
  save(outliers, file=config["out_het_outliers"])
  
  png(paste(config["out_het_plot_prefix"], "_%03d.png", sep=""), width=720, height=720)
  chromIntensityPlot(blData, scan.ids=outliers$scanID, chrom.ids=rep(chr.sel, nhet),
                     info=outliers$code, snp.exclude=snp.exclude, cex=0.25)
  dev.off()
}


# IBD connectivity
if (!is.na(config["ibd_con_file"])) {
  ibd.con <- getobj(config["ibd_con_file"])
  stopifnot(all(ibd.con %in% scanID))
  ncon <- length(ibd.con)
  
  png(paste(config["out_ibd_plot_prefix"], "_%03d.png", sep=""), width=720, height=720)
  chromIntensityPlot(blData, scan.ids=ibd.con, chrom.ids=rep(chr.sel, ncon),
                     info=rep("high connectivity", ncon), snp.exclude=snp.exclude, cex=0.25)
}


close(blData)
