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
required <- c("annot_scan_file", "annot_snp_file", "build", "nc_bl_file")
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
  scanID <- getScanID(scanAnnot)
  anom.id <- scanID[anom.ind]

  localID <- getVariable(scanAnnot, config["annot_scan_localIDCol"])[anom.ind]
  sex <- getVariable(scanAnnot, config["annot_scan_sexCol"])[anom.ind]
  comment <- comment[anom.ind]
  main <- paste("Scan", anom.id, "- Local", localID, "- Sex", sex, "- Chrom X\n", comment)
  
  bl.file <- config["nc_bl_file"]
  blnc <- NcdfIntensityReader(bl.file)
  blData <-  IntensityData(blnc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

  png.file <- file.path(paste(config["out_sexchrom_prefix"], "_%003d.png", sep=""))
  png(png.file, width=720, height=720)
  pseudoautoIntensityPlot(blData, scan.ids=anom.id, main=main, hg.build=config["build"],
                          snp.exclude=snp.excl)
  dev.off()
  
  close(blData)
}
