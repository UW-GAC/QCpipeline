##########
# LOH anomaly detection
# Usage: R --args config.file start end < anom_loh.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "bl_file", "geno_file",
              "out_anom_dir", "project")
optional <- c("scan_exclude_file", "out_eligible_snps")
default <- c(NA, "snps_eligible.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

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

bl <- IntensityReader(config["bl_file"])
blData <-  IntensityData(bl, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

geno <- GenotypeReader(config["geno_file"])
genoData <-  GenotypeData(geno, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# select SNPs
snp.ok <- getobj(config["out_eligible_snps"])
length(snp.ok)

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

BAF.filename <- file.path(config["out_anom_dir"], paste(config["project"], "BAF.filtered.all.RData", sep="."))
BAF <- getobj(BAF.filename)
## these are the anomalies already found by BAF
## these are taken out before finding LOH
## most important reason: more accurately find baseline non-anom mean/sd and median/mad
## BAF is data.frame with min set of variables including sample.num,chrom, left and right
##  where left and right are snp indices

known.anoms <- BAF[is.element(BAF$scanID, scan.ids),]

# detect anomalies
loh.anom <- anomDetectLOH(blData, genoData, scan.ids=scan.ids,
  chrom.ids=chrom.ids, snp.ids=snp.ok, known.anoms=known.anoms)

# save
for (n in names(loh.anom)) {
  saveas(loh.anom[[n]], paste(paste(config["project"],"LOH",n,sep="."),st,ed,sep="_"),
         config["out_anom_dir"])
}

close(blData)
close(genoData)
