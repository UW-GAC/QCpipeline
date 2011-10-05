##########
# BAF anomaly detection
# Usage: R --args config.file start end < anom_baf.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]
print(config.table)

# read start and end scan numbers
if (length(args) < 3) stop("missing start and end")
st <- as.integer(args[2])
ed <- as.integer(args[3])
if (st > ed) stop("Start is larger than End")

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
if (config["build"] == 36) {
  data(HLA.hg18)
  hla <- chrom == "6" & pos >= HLA.hg18$start.base & pos <= HLA.hg18$end.base
  data(pseudoautosomal.hg18)
  xtr <- chrom == "X" & pos >= pseudoautosomal.hg18["X.XTR", "start.base"] & 
    pos <= pseudoautosomal.hg18["X.XTR", "end.base"]
  data(centromeres.hg18)
  centromeres <- centromeres.hg18
} else if (config["build"] == 37) {
  data(HLA.hg19)
  hla <- chrom == "6" & pos >= HLA.hg19$start.base & pos <= HLA.hg19$end.base
  data(pseudoautosomal.hg19)
  xtr <- chrom == "X" & pos >= pseudoautosomal.hg19["X.XTR", "start.base"] & 
    pos <= pseudoautosomal.hg19["X.XTR", "end.base"]
  data(centromeres.hg19)
  centromeres <- centromeres.hg19
}

#ignore includes intensity-only and failed snps
ignore <- getVariable(snpAnnot, config["annot_snp_missingCol"]) == 1 
snp.exclude <- ignore | hla | xtr
snp.ok <- snpID[!snp.exclude]

if (as.logical(config["chromXY"])) {
  # remove only IO snps for XY region
  io <- getVariable(snpAnnot, config["annot_snp_IntensityOnlyCol"])
  if (is.null(io)) {
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
                               chrom.ids=24, snp.ids=sel.XY)
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
