##########
# Missing call rates by SNP and scan
# Usage: R --args config.file < missing.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)

# missing.n1
miss.by.snp <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)
save(miss.by.snp, file=config["out_n1_file"])
missing.n1 <- miss.by.snp$missing.fraction
stopifnot(allequal(names(missing.n1), snpID))

# summary table
chr.type <- getChromosome(snpAnnot, char=TRUE)
chr.type[is.element(chr.type, 1:22)] <- "A"
np <- table(chr.type)

# fraction of intensity only per chrom
if (!is.na(config["annot_snp_IntensityOnlyCol"])) {
  snp.io <- getVariable(snpAnnot, config["annot_snp_IntensityOnlyCol"])
} else {
  snp.io <- NULL
}
if (!is.null(snp.io) & sum(snp.io) > 0) {
  xm <- as.matrix(table(chr.type, snp.io))
  io <- xm[,"1"] / rowSums(xm)
  nonio <- snp.io == 0
} else {
  snp.io <- NULL
  nonio <- rep(TRUE, nrow(snpAnnot))
}

# fraction of non-intensity only snps that are tech failures
xm <- as.matrix(table(chr.type[nonio], missing.n1[nonio] == 1))
tf <- xm[,"TRUE"] / rowSums(xm)

# fraction of snps with missing.n1>0.05 among passed SNPs
goodsnp <- nonio & missing.n1 < 1
xm <- as.matrix(table(chr.type[goodsnp], missing.n1[goodsnp]>0.05))
mf <- xm[,"TRUE"] / rowSums(xm)

if (!is.null(snp.io)) {
  snptbl <- rbind(np, io, tf, mf)
  row.names(snptbl) <- c("number of probes", "intensity-only", "SNP tech failures", "missing>0.05")
} else {
  snptbl <- rbind(np, tf, mf)
  row.names(snptbl) <- c("number of probes", "SNP tech failures", "missing>0.05")
}
save(snptbl, file=config["out_snp_summary"])


# missing.e1 - snps to exclude
snp.exclude <- snpID[missing.n1 == 1]
length(snp.exclude)

# missing.e1
miss.by.scan <- missingGenotypeByScanChrom(genoData, snp.exclude=snp.exclude)
missing.e1 <- miss.by.scan$missing.fraction
stopifnot(allequal(names(missing.e1), scanID))

# autosomal
auto <- colnames(miss.by.scan$missing.counts) %in% 1:22
missa <- rowSums(miss.by.scan$missing.counts[,auto]) /
  sum(miss.by.scan$snps.per.chr[auto])
stopifnot(allequal(names(missa), scanID))
miss.by.scan$missing.fraction.auto <- missa

# X chromosome
missx <- miss.by.scan$missing.counts[,"X"] / miss.by.scan$snps.per.chr["X"]
stopifnot(allequal(names(missx), scanID))
miss.by.scan$missing.fraction.xchr <- missx
save(miss.by.scan, file=config["out_e1_file"])

# plot
pdf(config["out_e1_hist"], width=6, height=6)
hist(missing.e1, xlab="Missing call rate by sample", main="")
dev.off()

# do round 2?
if (as.logical(config["round2"])) {
  # are there any scans with missing.e1 >= 0.05?
  badscans <- missing.e1 >= 0.05
  
  # missing.n2
  if (any(badscans)) {
    scan.exclude <- unique(c(scan.exclude, scanID[badscans]))
    print(length(scan.exclude))
    
    miss.by.snp <- missingGenotypeBySnpSex(genoData, scan.exclude=scan.exclude)
    save(miss.by.snp, file=config["out_n2_file"])
    missing.n2 <- miss.by.snp$missing.fraction
    stopifnot(allequal(names(missing.n2), snpID))
  } else {
    missing.n2 <- missing.n1
  }

  snp.exclude <- snpID[missing.n2 >= 0.05]
  print(length(snp.exclude))
  
  # missing.e2
  miss.by.scan <- missingGenotypeByScanChrom(genoData, snp.exclude=snp.exclude)
  missing.e2 <- miss.by.scan$missing.fraction
  stopifnot(allequal(names(missing.e2), scanID))
  
  # autosomal
  auto <- colnames(miss.by.scan$missing.counts) %in% 1:22
  missa <- rowSums(miss.by.scan$missing.counts[,auto]) /
  sum(miss.by.scan$snps.per.chr[auto])
  stopifnot(allequal(names(missa), scanID))
  miss.by.scan$missing.fraction.auto <- missa
  
  # X chromosome
  missx <- miss.by.scan$missing.counts[,"X"] / miss.by.scan$snps.per.chr["X"]
  stopifnot(allequal(names(missx), scanID))
  miss.by.scan$missing.fraction.xchr <- missx
  save(miss.by.scan, file=config["out_e2_file"])

  # plot
  pdf(config["out_e2_hist"], width=6, height=6)
  hist(missing.e2, xlab="Missing call rate by sample", main="")
  dev.off()

}

close(genoData)
