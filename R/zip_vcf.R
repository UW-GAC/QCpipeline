library(Rsamtools)

file <- commandArgs(trailingOnly=TRUE)
zipped <- bgzip(file, paste(file, ".gz", sep=""))
idx <- indexTabix(zipped, "vcf")
