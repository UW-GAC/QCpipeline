##########
# Plot BAF/LRR for sex chrom anomalies identified by CIDR
# Usage: R --args config.file < sexchrom_anom_plot.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "build", "bl_file")
optional <- c("annot_scan_commentCol", "annot_scan_localIDCol", "annot_scan_sexCol",
              "annot_snp_IntensityOnlyCol", "out_sexchrom_prefix")
default <- c("Comment", "local.scanID", "sex", NA, "sexchrom_anom")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
comment <- getVariable(scanAnnot, config["annot_scan_commentCol"])
anomlist <- c("XXX", "XXY", "XYY", "XO")
anom.ind <- vector()
for (a in anomlist) {
  anom.ind <- c(anom.ind, grep(a, comment))
}
anom.ind <- unique(anom.ind)
# get normal samples for comparison
norm.ind <- setdiff(1:nrow(scanAnnot), anom.ind)
# take out XY/XO
xyxo <- grep("XY/XO", comment)
anom.ind <- setdiff(anom.ind, xyxo)
length(anom.ind)

(snpAnnot <- getobj(config["annot_snp_file"]))

# don't plot intensity-only
if (!is.na(config["annot_snp_IntensityOnlyCol"])) {
  snp.io <- getVariable(snpAnnot, config["annot_snp_IntensityOnlyCol"])
  snp.excl <- snp.io == 1
} else {
  snp.excl <- NULL
}

if (length(anom.ind) > 0) {
  anom.id <- getScanID(scanAnnot)[anom.ind]
  localID <- getVariable(scanAnnot, config["annot_scan_localIDCol"])[anom.ind]
  sex <- getVariable(scanAnnot, config["annot_scan_sexCol"])[anom.ind]
  anom.comment <- comment[anom.ind]
  main <- paste("Scan", anom.id, "- Local", localID, "- Sex", sex, "- Chrom X\n", anom.comment)
  
  bl.file <- config["bl_file"]
  blnc <- IntensityReader(bl.file)
  blData <-  IntensityData(blnc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  png.file <- file.path(paste(config["out_sexchrom_prefix"], "_%003d.png", sep=""))
  png(png.file, width=720, height=720)
  pseudoautoIntensityPlot(blData, scan.ids=anom.id, main=main, hg.build=config["build"],
                          snp.exclude=snp.excl)
  dev.off()

  # normal samples - 3 of each sex
  sex <- getVariable(scanAnnot, config["annot_scan_sexCol"])
  m.ind <- sample(intersect(norm.ind, which(sex == "M")), 3)
  f.ind <- sample(intersect(norm.ind, which(sex == "F")), 3)
  samp.ind <- c(f.ind, m.ind)
  sex <- sex[samp.ind]
  
  norm.id <- getScanID(scanAnnot)[samp.ind]
  localID <- getVariable(scanAnnot, config["annot_scan_localIDCol"])[samp.ind]
  norm.comment <- comment[samp.ind]
  main <- paste("Scan", norm.id, "- Local", localID, "- Sex", sex, "- Chrom X\n", norm.comment)
  
  png.file <- file.path(paste(config["out_sexchrom_prefix"], "_norm_%003d.png", sep=""))
  png(png.file, width=720, height=720)
  pseudoautoIntensityPlot(blData, scan.ids=norm.id, main=main, hg.build=config["build"],
                          snp.exclude=snp.excl)
  dev.off()
  
  close(blData)
}
