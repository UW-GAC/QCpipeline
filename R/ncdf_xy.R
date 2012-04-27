##########
# XY netCDF file
# Usage: R --args config.file <test> < ncdf_xy.R
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
              "nc_xy_file", "raw_path")
optional <- c("annot_scan_fileCol", "annot_scan_nameCol", "annot_snp_nameCol",
              "nc_xy_checkFile", "nc_xy_diagFile", "raw_xCol", "raw_yCol",
              "raw_colTotal", "raw_sampleCol", "raw_scanNameInFile", "raw_sepType",
              "raw_skipNum", "raw_snpCol")
default <- c("file", "Sample.Name", "rsID", "nc_xy_check.RData", "nc_xy_diag.RData",
             14, 15, 19, 2, 1, ",", 11, 1)
config <- setConfigDefaults(config, required, optional, default)
print(config)

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
# Create ncdf and populate with snp annotation
##########
(snpAnnot <- getobj(config["annot_snp_file"]))
snpdf <- getVariable(snpAnnot, c("snpID", "chromosome", "position"))

ncfile <- config["nc_xy_file"]
system.time({
  ncdfCreate(snp.annotation = snpdf,
             ncdf.filename = ncfile,
             variables = c("X","Y"),
             n.samples = nsamp,
             precision = "single",
             array.name = config["array_name"],
             genome.build = config["array_build"])
})

# check
# validity methods of NcdfIntensityReader and IntensityData take care of
# most checks we need
(nc <- NcdfIntensityReader(ncfile))
stopifnot(nsnp(nc) == nrow(snpAnnot))
stopifnot(nscan(nc) == nsamp)
stopifnot(getAttribute(nc, "array_name") == config["array_name"])
stopifnot(getAttribute(nc, "genome_build") == config["array_build"])
data <- IntensityData(nc, snpAnnot=snpAnnot)
close(data)


##########
# Add data to ncdf
##########
snpdf <- getVariable(snpAnnot, c("snpID", config["annot_snp_nameCol"]))
names(snpdf) <- c("snpID", "snpName")

scandf <- getVariable(scanAnnot, c("scanID", config["annot_scan_nameCol"], config["annot_scan_fileCol"]))
names(scandf) <- c("scanID", "scanName", "file")

col.nums <- as.integer(c(config["raw_snpCol"], config["raw_xCol"], config["raw_yCol"]))
names(col.nums) <- c("snp", "x", "y")
if (config["raw_scanNameInFile"] == 1) {
  col.nums <- append(col.nums, as.integer(config["raw_sampleCol"]), after=1)
  names(col.nums)[2] <- "sample"
}
skip.num <- as.integer(config["raw_skipNum"])
col.total <- as.integer(config["raw_colTotal"])
scan.name.in.file <- as.integer(config["raw_scanNameInFile"])

system.time({
  res <- ncdfAddData(path = config["raw_path"], ncdf.filename = ncfile,
                     snp.annotation = snpdf, scan.annotation = scandf,
                     sep.type=config["raw_sepType"], skip.num=skip.num,
                     col.total=col.total,
                     col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                     diagnostics.filename=config["nc_xy_diagFile"])
})

########################################
# MANUAL REVIEW - DIAGNOSTICS
########################################
names(res)
table(res$read.file, exclude=NULL)
stopifnot(all(res$read.file == 1))

table(res$row.num, exclude=NULL)
stopifnot(all(res$row.num == nrow(snpAnnot)))

table(unlist(lapply(res$samples, length)),exclude=NULL)  
table(res$sample.match, exclude=NULL)

table(unlist(lapply(res$missg, length)), exclude=NULL) 
unique(unlist(res$missg))

table(res$snp.chk, exclude=NULL)
stopifnot(all(res$snp.chk == 1))

table(res$chk, exclude=NULL)
stopifnot(all(res$chk == 1))

########################################
# MANUAL REVIEW - DATA VALUES
########################################
(nc <- NcdfIntensityReader(ncfile))
data <- IntensityData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
scanID <- getScanID(data)
range(scanID)

n <- nsamp-4   # number of samples minus 4 to get last 5 samples

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
system.time({
  res <- ncdfCheckIntensity(path = config["raw_path"], ncdf.filename = ncfile,
                            snp.annotation = snpdf, scan.annotation = scandf,
                            sep.type=config["raw_sepType"], skip.num=skip.num,
                            col.total=col.total,
                            col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                            check.scan.index=1:nsamp, n.scans.loaded=nsamp,
                            diagnostics.filename=config["nc_xy_checkFile"])
})

########################################
# MANUAL REVIEW - DIAGNOSTICS
########################################
table(res$read.file, exclude=NULL)
stopifnot(all(res$read.file == 1))

table(res$row.num, exclude=NULL)
stopifnot(all(res$row.num == nrow(snpAnnot)))

table(unlist(lapply(res$sample.names, length)),exclude=NULL)
table(res$sample.match, exclude=NULL)

table(res$snp.chk, exclude=NULL)
stopifnot(all(res$snp.chk == 1))

table(res$chk, exclude=NULL)
stopifnot(all(res$chk == 1))

table(res$snp.order, exclude=NULL)
stopifnot(all(res$snp.order == 1))

table(res$inten.chk$x, exclude=NULL)
stopifnot(all(res$inten.chk$x == 1))

table(res$inten.chk$y, exclude=NULL)
stopifnot(all(res$inten.chk$y == 1))
