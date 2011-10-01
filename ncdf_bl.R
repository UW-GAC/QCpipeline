##########
# BAF/LRR netCDF file
# Usage: R --args config.file <test> < ncdf_bl.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

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

ncfile <- config["nc_bl_file"]
system.time({
  ncdfCreate(snp.annotation = snpdf,
             ncdf.filename = ncfile,
             variables = c("BAlleleFreq","LogRRatio"),
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

col.nums <- as.integer(c(config["raw_snpCol"], config["raw_bafCol"], config["raw_lrrCol"]))
names(col.nums) <- c("snp", "ballelefreq", "logrratio")
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
                     diagnostics.filename=config["nc_bl_diagFile"])
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
baf <- getBAlleleFreq(data, snp=c(1,-1), scan=c(n,5))
dim(baf)
nsnp <- nrow(snpAnnot)
stopifnot(dim(baf) == c(nsnp, 5))
baf[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
baf[(nsnp-4):nsnp,]  # last 5 snps

lrr <- getLogRRatio(data, snp=c(1,-1), scan=c(n,5))
dim(lrr)
nsnp <- nrow(snpAnnot)
stopifnot(dim(lrr) == c(nsnp, 5))
lrr[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
lrr[(nsnp-4):nsnp,]  # last 5 snps
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
                            diagnostics.filename=config["nc_bl_checkFile"])
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

table(res$inten.chk$ballelefreq, exclude=NULL)
stopifnot(all(res$inten.chk$ballelefreq == 1))

table(res$inten.chk$logrratio, exclude=NULL)
stopifnot(all(res$inten.chk$logrratio == 1))

