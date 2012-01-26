##########
# HWE
# Usage: R --args config.file < hwe.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

nc <- NcdfGenotypeReader(config["nc_subj_geno_file"])
scanAnnot <- getobj(config["annot_scan_file"])
# take subset of annotation to match netCDF
scanAnnot <- scanAnnot[match(getScanID(nc), getScanID(scanAnnot)),]
genoData <- GenotypeData(nc, scanAnnot=scanAnnot)
scanID <- getScanID(genoData)

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
            scan.exclude = scan.exclude, outfile=config["out_hwe_prefix"])
