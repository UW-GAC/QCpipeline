##########
# Sample quality checks
# Usage: R --args config.file < sample_quality.R
##########

library(GWASTools)
library(QCpipeline)
library(ggplot2)
library(reshape2)
library(plyr)
sessionInfo()

theme_set(theme_bw())

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_scan_raceCol", "annot_snp_file",
              "baf_mean_file", "baf_sd_file", "bl_file")
optional <- c("annot_scan_hetACol", "annot_snp_missingCol", "ibd_con_file",
              "out_baf_sd_boxplot", "out_het_boxplot", "out_flagged_file", "out_plot_prefix",
              "plot_all_unknown", "race_unknown", "range_het", "range_sd",
              "out_baf_density_prefix")
default <- c("het.A", "missing.n1", NA,
             "baf_sd.pdf", "het_outl.pdf", "qual_check.RData", "qual_check",
             TRUE, NA, 1.5, 1.5,
             "baf_density")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

# exclude missing SNPs from plots
(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)
snp.exclude <- snpID[getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1]

# IntensityData
bl <- IntensityReader(config["bl_file"])
blData <- IntensityData(bl, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# list for storing scanIDs of flagged samples
flagged <- list()


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
  flagged[["asymmetry"]] <- asym
}


# IBD connectivity
if (!is.na(config["ibd_con_file"])) {
  ibd.con <- getobj(config["ibd_con_file"])
  stopifnot(all(ibd.con %in% scanID))
  ncon <- length(ibd.con)

  if (ncon > 0) {
    flagged[["ibd connectivity"]] <- ibd.con
  }
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
    flagged[["race unknown"]] <- unknown.ids
  }
}

# definition of outliers 
range.het <- as.numeric(config["range_het"])

pdf(config["out_het_boxplot"], width=6, height=6)
bp <- boxplot(samp$het.A ~ as.factor(samp$race), range=range.het, varwidth=TRUE, las=2,
              main=config["annot_scan_raceCol"], ylab="Autosomal heterozygosity")
dev.off()
## ggplot cannot use range unless you precompute stats yourself.
#p <- ggplot(samp, aes(x=race, y=het.A)) + geom_boxplot() + ggtitle(config["annot_scan_raceCol"]) + ylab("Autosomal heterozygosity")
#ggsave(config["out_het_boxplot"], plot=p, width=6, height=6)

outliers <- samp[samp$het.A %in% bp$out,]
(nhet <- nrow(outliers))
table(outliers$race)

if (nhet > 0) {
  # high or low?
  outliers$code <- NA
  races <- unique(samp$race)
  for (race in races) {
    med <- median(samp$het.A[samp$race == race], na.rm=TRUE)
    outliers$code[outliers$race == race & outliers$het.A > med] <- "high"
    outliers$code[outliers$race == race & outliers$het.A < med] <- "low"
  }

  if (sum(outliers$code == "high") > 0) {
    flagged[["high heterozygosity"]] <- outliers$scanID[outliers$code == "high"]
  }
  if (sum(outliers$code == "low") > 0) {
    flagged[["low heterozygosity"]] <- outliers$scanID[outliers$code == "low"]
  }
}


# SD of BAF
baf.sd <- getobj(config["baf_sd_file"])

# definition of outliers - range.sd*i(nterquartile range)
range.sd <- as.numeric(config["range_sd"])

# create a matrix
bsd <- do.call(cbind, baf.sd)
colnames(bsd) <- names(baf.sd)

# can't use range in ggplot without computing stats yourself, so keep this as base R
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

flagged[["BAF SD"]] <- baf.outl

save(flagged, file=config["out_flagged_file"])


# review BAF/LRR of chr 1 for all samples for which >1 chromosome is an outlier to look for contamination or other artifact
# we don't necessarily need to look at the chromosome(s) flagged since we are looking for a sample-wide artifact
# use plot without any color-coding

# if we don't have chromosome 1, just take the first chrom
if ("1" %in% names(baf.sd)) {
  chr.sel <- 1
} else {
  chr.sel <- names(baf.sd)[1]
}

# scans to plot
tmp <- melt(flagged)
names(tmp) <- c("scanID", "tag")
tags <- ddply(tmp, .(scanID), summarise, tags=paste(tag, collapse="; "))
table(tags$tags)

ndigits <- floor(log10(nrow(tags)) + 1)

png(paste(config["out_plot_prefix"], "_%0", ndigits, "d.png", sep=""), width=720, height=720)
chromIntensityPlot(blData, scan.ids=tags$scanID, chrom.ids=rep(chr.sel, nrow(tags)),
                   info=tags$tags, snp.exclude=snp.exclude, cex=0.25, ideogram=FALSE)
dev.off()

# review density plot of BAF for all samples for which >1 chromosome is an outlier
# one line plotted per chromosome
# contamination will appear as intermediate bumps between 0 and 0.5 and between 0.5 and 1 in most chromsomes
# plot density of BAF by chromosome to look for contamination
dat <- getVariable(snpAnnot, c("snpID", "chromosome"))
sel <- dat$chromosome %in% 1:22
dat <- dat[sel, ]
snp.idx <- range(which(sel))
snp.start <- snp.idx[1]
snp.count <- snp.idx[2] - snp.idx[1] + 1


fname <- paste0(config["out_baf_density_prefix"], "_%0", ndigits, "d.png")
adj <- 1/20 # bandwidth esque
for (i in seq_along(tags$scanID)){
  
  if (i %% 10 == 0) message(sprintf("%s of %s flagged samples", i, nrow(tags)))
  scan.id <- tags$scanID[i]
  scan.idx <- which(scanID %in% scan.id)
  
  # read in BAF for this scan
  dat$baf <- getBAlleleFreq(blData, scan=c(scan.idx, 1), snp=c(snp.start, snp.count))  
  p <- ggplot(dat, aes(x=baf)) +
    geom_line(aes(group=chromosome), alpha=0.3, stat="density", color="black", adjust=adj, na.rm=TRUE, size=0.2) +
    #geom_line(stat="density", color="red", adjust=adj, na.rm=TRUE) +
    theme(axis.text.y=element_blank()) +
    ggtitle(sprintf("%s - %s", scan.id, tags$tags[i]))
  ggsave(sprintf(fname, i), plot=p, width=6, height=6)
}

close(blData)
