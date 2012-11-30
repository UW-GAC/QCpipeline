##########
# Compress a VCF file with bgzip and create an index file
# Usage: R --args file < zip_vcf.R
##########

library(Rsamtools)

file <- commandArgs(trailingOnly=TRUE)
zipped <- bgzip(file, paste(file, ".gz", sep=""))
idx <- indexTabix(zipped, "vcf")
