##########
# BAF anomaly detection
# Usage: R --args config.file start end maf < anom_baf.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "build", "bl_file", "geno_file",
              "out_anom_dir", "out_baf_med_file", "project")
optional <- c("annot_snp_IntensityOnlyCol", "annot_snp_missingCol", "chromXY",
              "out_afreq_file", "out_eligible_snps", "scan_exclude_file", "snp_exclude_file")
default <- c(NA, "missing.n1", FALSE, "allele_freq.RData", "snps_eligible.RData", NA, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# read start and end scan numbers
if (length(args) < 3) stop("missing start and end")
st <- as.integer(args[2])
ed <- as.integer(args[3])
if (st > ed) stop("Start is larger than End")

# read MAF threshold
if (length(args) > 3) {
  maf <- as.numeric(args[4])
} else {
  maf <- NULL
}
print(maf)

# create GenotypeData and IntensityData objects
(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

bl <- IntensityReader(config["bl_file"])
blData <-  IntensityData(bl, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

geno <- GenotypeReader(config["geno_file"])
genoData <-  GenotypeData(geno, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

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

# any user-specified snps?
# any to exclude?
if (!is.na(config["snp_exclude_file"])){
  snp_exclude <- getobj(config["snp_exclude_file"])
  usr <- getSnpID(snpAnnot) %in% snp_exclude
} else {
  usr <- rep(FALSE, length(getSnpID(snpAnnot)))
}

#ignore includes intensity-only and failed snps
ignore <- getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1 
snp.exclude <- ignore | hla | xtr | gap | usr
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
save(snp.ok, file=config["out_eligible_snps"])

if (as.logical(config["chromXY"])) {
  # remove only IO snps for XY region
  if (!is.na(config["annot_snp_IntensityOnlyCol"])) {
    io <- getVariable(snpAnnot, config["annot_snp_IntensityOnlyCol"])
  } else {
    io <- rep(0, length(snpID))
  }
  sXY <- is.element(chrom, "XY") & is.element(io, 0)
  sel.XY <- snpID[sXY]
}

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

# select scans and chromosomes
scan.ids <- scanID[st:ed]
scan.ids <- setdiff(scan.ids, scan.exclude)
chrom.ids <- intersect(unique(getChromosome(snpAnnot)), 1:23)
 
## identifying some messy
med.baf.sd <- getobj(config["out_baf_med_file"])
med.baf.sd <- med.baf.sd[!(med.baf.sd$scanID %in% scan.exclude),]
low.qual.ids <- med.baf.sd$scanID[med.baf.sd$med.sd > 0.05]
if (length(low.qual.ids) == 0) low.qual.ids <- NULL

# segmentation
baf.seg <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
                          chrom.ids=chrom.ids, snp.ids=snp.ok)

if (as.logical(config["chromXY"])) {
  baf.seg.XY <- anomSegmentBAF(blData, genoData, scan.ids=scan.ids,
                               chrom.ids=XYchromCode(blData), snp.ids=sel.XY)
  baf.seg <- rbind(baf.seg, baf.seg.XY)
}

# filter anomalies
baf.anom <- anomFilterBAF(blData, genoData, segments=baf.seg, 
  snp.ids=snp.ok, centromere=centromeres, low.qual.ids=low.qual.ids)

# save
for (n in names(baf.anom)) {
  saveas(baf.anom[[n]], paste(paste(config["project"],"BAF",n,sep="."),st,ed,sep="_"),
         config["out_anom_dir"])
}

close(blData)
close(genoData)
