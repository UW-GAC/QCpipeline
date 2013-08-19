##########
# Genotype netCDF file
# Usage: R --args config.file <test> < ncdf_geno.R
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
              "nc_geno_file", "raw_path")
optional <- c("annot_scan_fileCol", "annot_scan_nameCol", "annot_snp_nameCol",
              "nc_geno_checkFile", "nc_geno_diagFile", "raw_a1Col", "raw_a2Col","raw_genoCol",
              "raw_colTotal", "raw_sampleCol", "raw_scanNameInFile", "raw_sepType",
              "raw_skipNum", "raw_snpCol")
default <- c("file", "Sample.Name", "rsID", "nc_geno_check.RData",
             "nc_geno_diag.RData", 10, 11, NA, 19, 2, 1, ",", 11, 1)
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

ncfile <- config["nc_geno_file"]
system.time({
  ncdfCreate(snp.annotation = snpdf,
             ncdf.filename = ncfile,
             variables = "genotype",
             n.samples = nsamp,
             array.name = config["array_name"],
             genome.build = config["array_build"])
})

# check
# validity methods of NcdfGenotypeReader and GenotypeData take care of
# most checks we need
(nc <- NcdfGenotypeReader(ncfile))
stopifnot(nsnp(nc) == nrow(snpAnnot))
stopifnot(nscan(nc) == nsamp)
stopifnot(getAttribute(nc, "array_name") == config["array_name"])
stopifnot(getAttribute(nc, "genome_build") == config["array_build"])
data <- GenotypeData(nc, snpAnnot=snpAnnot)
close(data)


##########
# Add data to ncdf
##########
snpdf <- getVariable(snpAnnot, c("snpID", config["annot_snp_nameCol"]))
names(snpdf) <- c("snpID", "snpName")

scandf <- getVariable(scanAnnot, c("scanID", config["annot_scan_nameCol"], config["annot_scan_fileCol"]))
names(scandf) <- c("scanID","scanName" ,"file")

skip.num <- as.integer(config["raw_skipNum"])
col.total <- as.integer(config["raw_colTotal"])
scan.name.in.file <- as.integer(config["raw_scanNameInFile"])

## if alleles are in separate columns
if (is.na(config["raw_genoCol"])) {
	col.nums <- as.integer(c(config["raw_snpCol"], config["raw_a1Col"], config["raw_a2Col"]))
	names(col.nums) <- c("snp", "a1", "a2")
	
## if alleles are in the same column
} else {
	col.nums <- as.integer(c(config["raw_snpCol"], config["raw_genoCol"]))
	names(col.nums) <- c("snp", "geno")
	
	# in test mode - print out col.nums argument
	print(col.nums)
}
	
if (config["raw_scanNameInFile"] == 1) {
  col.nums <- append(col.nums, as.integer(config["raw_sampleCol"]), after=1)
  names(col.nums)[2] <- "sample"
}

system.time({
  res <- ncdfAddData(path = config["raw_path"], ncdf.filename = ncfile,
					 snp.annotation = snpdf, scan.annotation = scandf,
					 sep.type=config["raw_sepType"], skip.num=skip.num,
					 col.total=col.total,
					 col.nums=col.nums, scan.name.in.file=scan.name.in.file,
					 diagnostics.filename=config["nc_geno_diagFile"])
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

########################################
# MANUAL REVIEW - DATA VALUES
########################################
(nc <- NcdfGenotypeReader(ncfile))
data <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
scanID <- getScanID(data)
range(scanID)

n <- nsamp-4   # number of samples minus 4 to get last 5 samples
geno <- getGenotype(data, snp=c(1,-1), scan=c(n,5))
dim(geno) #snp by sample
nsnp <- nrow(snpAnnot)
stopifnot(dim(geno) == c(nsnp, 5))
table(geno[,1], exclude=NULL)
geno[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
geno[(nsnp-4):nsnp,] # last 5 snps
close(data)


##########
# Check ncdf file
##########
system.time({
  res <-ncdfCheckGenotype(path = config["raw_path"], ncdf.filename = ncfile,
                          snp.annotation = snpdf, scan.annotation = scandf,
                          sep.type=config["raw_sepType"], skip.num=skip.num,
                          col.total=col.total,
                          col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                          check.scan.index=1:nsamp, n.scans.loaded=nsamp,
                          diagnostics.filename=config["nc_geno_checkFile"])
})

########################################
# MANUAL REVIEW - DIAGNOSTICS
########################################
table(res$geno.chk[1:nsamp])
table(res$read.file, exclude=NULL)
stopifnot(all(res$read.file == 1))

table(res$row.num, exclude=NULL)
stopifnot(all(res$row.num == nrow(snpAnnot)))

table(unlist(lapply(res$sample.names, length)),exclude=NULL)
table(res$sample.match, exclude=NULL)

table(unlist(lapply(res$missg, length)), exclude=NULL) 
unique(unlist(res$missg))

table(res$snp.chk, exclude=NULL)
stopifnot(all(res$snp.chk == 1))

table(res$chk, exclude=NULL)
stopifnot(all(res$chk == 1))

table(res$snp.order, exclude=NULL)
stopifnot(all(res$snp.order == 1))

table(res$geno.chk, exclude=NULL)
stopifnot(all(res$geno.chk == 1))


# last check that annotation matches netCDF
(nc <- NcdfGenotypeReader(ncfile))
data <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
close(data)
