##########
# HWE
# Usage: R --args config.file chr.start chr.end < hwe.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "nc_geno_file", "scan_include_file")
optional <- c("out_hwe_prefix", "scan_chrom_filter")
default <- c("hwe", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

data <- GenotypeReader(config["nc_geno_file"])
scanAnnot <- getobj(config["annot_scan_file"])
# take subset of annotation to match netCDF
scanAnnot <- scanAnnot[match(getScanID(data), getScanID(scanAnnot)),]
genoData <- GenotypeData(data, scanAnnot=scanAnnot)
scanID <- getScanID(genoData)
chrom <- getChromosome(genoData)

if (length(args) == 3) {
  chr.start <- as.integer(args[2])
  chr.end <- as.integer(args[3])
  chromosome.set <- intersect(seq(chr.start, chr.end), unique(chrom))
} else {
  chromosome.set <- NULL
}
chromosome.set

# are there any scans to exclude?
scan.include <- getobj(config["scan_include_file"])
stopifnot(all(scan.include %in% scanID))
scan.exclude <- setdiff(scanID, scan.include)
length(scan.exclude)

# scan-chromosome filter
if (!is.na(config["scan_chrom_filter"])) {
  scan.filt <- getobj(config["scan_chrom_filter"])
} else {
  scan.filt <- NULL
}

gwasExactHW(genoData, scan.chromosome.filter=scan.filt,
            scan.exclude = scan.exclude, outfile=config["out_hwe_prefix"],
            chromosome.set=chromosome.set)
