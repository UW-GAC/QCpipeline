##########
# HWE
# Usage: R --args config.file chromosome < hwe.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# read chromosome
if (length(args) > 1) 
    (chromosome <- as.integer(args[2]))

# check config and set defaults
required <- c("annot_scan_file", "geno_file", "scan_include_file")
optional <- c("out_hwe_prefix", "block_size")
default <- c("hwe", 5000)
config <- setConfigDefaults(config, required, optional, default)
print(config)

data <- GenotypeReader(config["geno_file"])
scanAnnot <- getobj(config["annot_scan_file"])
# take subset of annotation to match netCDF
scanAnnot <- scanAnnot[match(getScanID(data), getScanID(scanAnnot)),]
genoData <- GenotypeData(data, scanAnnot=scanAnnot)

# are there any scans to exclude?
scan.include <- getobj(config["scan_include_file"])
scanID <- getScanID(genoData)
stopifnot(all(scan.include %in% scanID))
scan.exclude <- setdiff(scanID, scan.include)
length(scan.exclude)

if (length(args) > 1) {
    chrom.range <- range(which(getChromosome(genoData) %in% chromosome))
} else {
    chrom.range <- range(which(getChromosome(genoData) %in% 1:22))
}

block.size <- as.integer(config["block_size"])

hwe <- exactHWE(genoData, 
                scan.exclude = scan.exclude, 
                snpStart = chrom.range[1],
                snpEnd = chrom.range[2],
                block.size = block.size)

save(hwe, file=paste0(config["out_hwe_prefix"], "_chr", chromosome, ".RData"))

