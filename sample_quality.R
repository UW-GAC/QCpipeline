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

# create a matrix
bsd <- matrix(nrow=nrow(baf.sd[[1]]), ncol=length(baf.sd),
              dimnames=list(rownames(baf.sd[[1]]), names(baf.sd)))
for (i in 1:ncol(bsd)) {
  bsd[,i] <- baf.sd[[i]][,1]
}
pdf(config["out_baf_sd_boxplot"], width=6, height=6)
bp <- boxplot(bsd[,colnames(bsd) %in% 1:22], xlab="chromosome", ylab="SD of BAF")
dev.off()

# get boxplot hi sd outliers by chromosome
out <- list()
for (i in names(baf.sd)[names(baf.sd) %in% 1:22]) {
  x <- baf.sd[[i]]
  med <- median(x)
  tmp <- boxplot(x, plot=FALSE)$out
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
                     code=rep("BAF SD", nbaf), snp.exclude=snp.exclude, cex=0.25)
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
                     code=rep("Asymmetry", nasym), snp.exclude=snp.exclude, cex=0.25)
  dev.off()
}


# Find boxplot outliers for heterozygosity, both hi and lo
# heterozygosity by race
samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"], config["annot_scan_hetACol"]))
names(samp) <- c("scanID", "race", "het.A")
samp <- samp[!is.na(samp$race),]

pdf(config["out_het_boxplot"], width=6, height=6)
bp <- boxplot(samp$het.A ~ as.factor(samp$race),
              xlab="race", ylab="Autosomal heterozygosity")
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
                     code=outliers$code, snp.exclude=snp.exclude, cex=0.25)
  dev.off()
}


# IBD connectivity
if (!is.na(config["ibd_con_file"])) {
  ibd.con <- getobj(config["ibd_con_file"])
  stopifnot(all(ibd.con %in% scanID))
  ncon <- length(ibd.con)
  
  png(paste(config["out_ibd_plot_prefix"], "_%03d.png", sep=""), width=720, height=720)
  chromIntensityPlot(blData, scan.ids=ibd.con, chrom.ids=rep(chr.sel, ncon),
                     code=rep("high connectivity", ncon), snp.exclude=snp.exclude, cex=0.25)
}


close(blData)
