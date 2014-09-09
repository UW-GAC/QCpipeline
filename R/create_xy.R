##########
# XY netCDF file
# Usage: R --args config.file <test> < create_xy.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "array_build", "array_name",
              "xy_file", "raw_path")
optional <- c("annot_scan_fileCol", "annot_scan_nameCol", "annot_snp_nameCol",
              "xy_checkFile", "xy_diagFile", "xy_file_type",
              "raw_qCol", "raw_xCol", "raw_yCol",
              "raw_colTotal", "raw_sampleCol", "raw_scanNameInFile", "raw_sepType",
              "raw_skipNum", "raw_snpCol")
default <- c("file", "Sample.Name", "rsID", "xy_check.RData", "xy_diag.RData", "gds",
             NA, 14, 15, 19, 2, 1, ",", 11, 1)
config <- setConfigDefaults(config, required, optional, default)
print(config)

getqual <- !is.na(config["raw_qCol"])

scanAnnot <- getobj(config["annot_scan_file"])

# check for test
if (length(args) > 1) {
  if ((args[2]) == "test") {
    message("testing first 5 scans")
    scanAnnot <- scanAnnot[1:5,] # testing
  }
}
nsamp <- nrow(scanAnnot)


##########
# Add data to ncdf
##########
ncfile <- config["xy_file"]
vars <- c("X","Y")
if (getqual) vars <- c("quality", vars)

(snpAnnot <- getobj(config["annot_snp_file"]))
snp.cols <- intersect(c("snpID", config["annot_snp_nameCol"], "chromosome", "position", 
                        "alleleA", "alleleB"),
                      getVariableNames(snpAnnot))
snpdf <- getVariable(snpAnnot, snp.cols)
names(snpdf)[2] <- "snpName"

scandf <- getVariable(scanAnnot, c("scanID", config["annot_scan_nameCol"], config["annot_scan_fileCol"]))
names(scandf) <- c("scanID", "scanName", "file")

col.nums <- as.integer(c(config["raw_snpCol"], config["raw_xCol"], config["raw_yCol"]))
names(col.nums) <- c("snp", "X", "Y")
if (config["raw_scanNameInFile"] == 1) {
  col.nums <- c(col.nums, "sample"=as.integer(config["raw_sampleCol"]))
}
if (getqual) {
  col.nums <- c(col.nums, "quality"=as.integer(config["raw_qCol"]))
}
skip.num <- as.integer(config["raw_skipNum"])
col.total <- as.integer(config["raw_colTotal"])
scan.name.in.file <- as.integer(config["raw_scanNameInFile"])

res <- createDataFile(path = config["raw_path"], filename = ncfile,
                      file.type = config["xy_file_type"],
                      variables = vars,
                      snp.annotation = snpdf, scan.annotation = scandf,
                      sep.type=config["raw_sepType"], skip.num=skip.num,
                      col.total=col.total,
                      col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                      array.name = config["array_name"],
                      genome.build = config["array_build"],
                      diagnostics.filename=config["xy_diagFile"])

########################################
# MANUAL REVIEW - DIAGNOSTICS
########################################
names(res)
table(res$read.file, useNA="ifany")
stopifnot(all(res$read.file == 1))

table(res$row.num, useNA="ifany")
stopifnot(all(res$row.num == nrow(snpAnnot)))

table(unlist(lapply(res$samples, length)),useNA="ifany")  
table(res$sample.match, useNA="ifany")

table(unlist(lapply(res$missg, length)), useNA="ifany") 
unique(unlist(res$missg))

table(res$snp.chk, useNA="ifany")
stopifnot(all(res$snp.chk == 1))

table(res$chk, useNA="ifany")
stopifnot(all(res$chk == 1))

########################################
# MANUAL REVIEW - DATA VALUES
########################################
(nc <- IntensityReader(ncfile))
data <- IntensityData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
scanID <- getScanID(data)
range(scanID)

n <- nsamp-4   # number of samples minus 4 to get last 5 samples

if (getqual) {
  q <- getQuality(data, snp=c(1,-1), scan=c(n,5))
  dim(q)
  nsnp <- nrow(snpAnnot)
  stopifnot(dim(q) == c(nsnp, 5))
  q[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
  q[(nsnp-4):nsnp,] # last 5 snps
}

x <- getX(data, snp=c(1,-1), scan=c(n,5))
dim(x)
nsnp <- nrow(snpAnnot)
stopifnot(dim(x) == c(nsnp, 5))
x[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
x[(nsnp-4):nsnp,] # last 5 snps

y <- getY(data, snp=c(1,-1), scan=c(n,5))
dim(y)
nsnp <- nrow(snpAnnot)
stopifnot(dim(y) == c(nsnp, 5))
y[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
y[(nsnp-4):nsnp,] # last 5 snps
close(data)


##########
# Check ncdf file
##########
res <- checkIntensityFile(path = config["raw_path"], filename = ncfile,
                          file.type = config["xy_file_type"],
                          snp.annotation = snpdf, scan.annotation = scandf,
                          sep.type=config["raw_sepType"], skip.num=skip.num,
                          col.total=col.total,
                          col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                          check.scan.index=1:nsamp, n.scans.loaded=nsamp,
                          diagnostics.filename=config["xy_checkFile"])

########################################
# MANUAL REVIEW - DIAGNOSTICS
########################################
table(res$read.file, useNA="ifany")
stopifnot(all(res$read.file == 1))

table(res$row.num, useNA="ifany")
stopifnot(all(res$row.num == nrow(snpAnnot)))

table(unlist(lapply(res$sample.names, length)),useNA="ifany")
table(res$sample.match, useNA="ifany")

table(res$snp.chk, useNA="ifany")
stopifnot(all(res$snp.chk == 1))

table(res$chk, useNA="ifany")
stopifnot(all(res$chk == 1))

table(res$snp.order, useNA="ifany")
stopifnot(all(res$snp.order == 1))

if (getqual) { 
  print(table(res$qs.chk, useNA="ifany"))
  stopifnot(all(res$chk == 1))
}

table(res$inten.chk$x, useNA="ifany")
stopifnot(all(res$inten.chk$x == 1))

table(res$inten.chk$y, useNA="ifany")
stopifnot(all(res$inten.chk$y == 1))

