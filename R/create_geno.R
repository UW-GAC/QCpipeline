##########
# Genotype netCDF file
# Usage: R --args config.file <test> <nbatches> <batch> < create_geno.R
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
              "geno_file", "raw_path")
optional <- c("annot_scan_fileCol", "annot_scan_nameCol", "annot_snp_nameCol",
              "geno_checkFile", "geno_diagFile", "geno_file_type",
              "raw_alleleCoding",
              "raw_a1Col", "raw_a2Col", "raw_genoCol",
              "raw_colTotal", "raw_sampleCol", "raw_scanNameInFile", "raw_sepType",
              "raw_skipNum", "raw_snpCol")
default <- c("file", "Sample.Name", "rsID", "geno_check.RData", "geno_diag.RData",
             "gds", "nucleotide",
             13, 14, NA, 14, NA, 0, ",", 12, 1)
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"])
ncfile <- config["geno_file"]
diagfile <- config["geno_diagFile"]
checkfile <- config["geno_checkFile"]

# check for test or batches
if (length(args) > 1) {
  if ((args[2]) == "test") {
    message("testing first 5 scans")
    scanAnnot <- scanAnnot[1:5,] # testing
  } else {
    nbatches <- as.integer(args[2])
    batch <- as.integer(args[3])
    stopifnot(batch <= nbatches)
    nsamp <- nrow(scanAnnot)
    n.per.batch <- ceiling(nsamp/nbatches)
    start <- n.per.batch*(batch-1) + 1
    end <- min(n.per.batch*batch, nsamp)
    scanAnnot <- scanAnnot[start:end,]
    ncfile <- paste0("batch", batch, ".", ncfile)
    diagfile <- paste0("batch", batch, ".", diagfile)
    checkfile <- paste0("batch", batch, ".", checkfile)
  }
}
nsamp <- nrow(scanAnnot)


##########
# Add data to ncdf
##########
(snpAnnot <- getobj(config["annot_snp_file"]))
snp.cols <- intersect(c("snpID", config["annot_snp_nameCol"], "chromosome", "position", 
                        "alleleA", "alleleB"),
                      getVariableNames(snpAnnot))
snpdf <- getVariable(snpAnnot, snp.cols)
names(snpdf)[2] <- "snpName"

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
}
	
if (config["raw_scanNameInFile"] == 1) {
  col.nums <- append(col.nums, as.integer(config["raw_sampleCol"]), after=1)
  names(col.nums)[2] <- "sample"
}

res <- createDataFile(path = config["raw_path"], filename = ncfile,
                      file.type = config["geno_file_type"],
                      variables = "genotype",
                      snp.annotation = snpdf, scan.annotation = scandf,
                      sep.type=config["raw_sepType"], skip.num=skip.num,
                      col.total=col.total,
                      col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                      allele.coding=config["raw_alleleCoding"],
                      array.name = config["array_name"],
                      genome.build = config["array_build"],
                      diagnostics.filename=diagfile)

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

########################################
# MANUAL REVIEW - DATA VALUES
########################################
(nc <- GenotypeReader(ncfile))
data <- GenotypeData(nc, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
scanID <- getScanID(data)
range(scanID)

n <- nsamp-1   # number of samples minus 1 to get last 2 samples
geno <- getGenotype(data, snp=c(1,-1), scan=c(n,2))
dim(geno) #snp by sample
nsnp <- nrow(snpAnnot)
stopifnot(dim(geno) == c(nsnp, 2))
table(geno[,1], useNA="ifany")
geno[round(nsnp/2):(round(nsnp/2)+10),] # middle 10 snps
geno[(nsnp-4):nsnp,] # last 5 snps
close(data)


##########
# Check ncdf file
##########
res <- checkGenotypeFile(path = config["raw_path"], filename = ncfile,
                         file.type = config["geno_file_type"],
                         snp.annotation = snpdf, scan.annotation = scandf,
                         sep.type=config["raw_sepType"], skip.num=skip.num,
                         col.total=col.total,
                         col.nums=col.nums, scan.name.in.file=scan.name.in.file,
                         check.scan.index=1:nsamp, n.scans.loaded=nsamp,
                         allele.coding=config["raw_alleleCoding"],
                         diagnostics.filename=checkfile)

########################################
# MANUAL REVIEW - DIAGNOSTICS
########################################
table(res$geno.chk[1:nsamp])
table(res$read.file, useNA="ifany")
stopifnot(all(res$read.file == 1))

table(res$row.num, useNA="ifany")
stopifnot(all(res$row.num == nrow(snpAnnot)))

table(unlist(lapply(res$sample.names, length)),useNA="ifany")
table(res$sample.match, useNA="ifany")

table(unlist(lapply(res$missg, length)), useNA="ifany") 
unique(unlist(res$missg))

table(res$snp.chk, useNA="ifany")
stopifnot(all(res$snp.chk == 1))

table(res$chk, useNA="ifany")
stopifnot(all(res$chk == 1))

table(res$geno.chk, useNA="ifany")
stopifnot(all(res$geno.chk == 1))

