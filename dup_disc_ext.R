##########
# Duplicate discordance with external data
# Usage: R --args config.file < dup_disc_ext.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "ext_annot_scan_file",
              "ext_annot_snp_file", "ext_nc_geno_file", "nc_geno_file")
optional <- c("annot_scan_subjectCol", "annot_snp_missingCol", "annot_snp_rsIDCol",
              "disc_scan_exclude_file", "ext_annot_scan_subjectCol",
              "ext_annot_snp_missingCol", "ext_annot_snp_rsIDCol",
              "ext_scan_exclude_file", "out_disc_file", "out_disc_plot")
default <- c("subjectID", "missing.n1", "rsID", NA, "subjectID", NA, "rsID",
             NA, "dup_disc_ext.RData", "dup_disc_ext.pdf")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"]); dim(scanAnnot)
snpAnnot <- getobj(config["annot_snp_file"])
nc <- NcdfGenotypeReader(config["nc_geno_file"])
# this might be a subject-level netCDF file, so subset annotation
scanAnnot <- scanAnnot[scanAnnot$scanID %in% getScanID(nc),]; dim(scanAnnot)
(genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot))

# external netcdf
ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
ext.snpAnnot <- getobj(config["ext_annot_snp_file"])
ext.nc <- NcdfGenotypeReader(config["ext_nc_geno_file"])
(extData <- GenotypeData(ext.nc, scanAnnot=ext.scanAnnot, snpAnnot=ext.snpAnnot))

# select common SNPs
# remove missing from study
snp.study <- pData(snpAnnot)[getVariable(snpAnnot, config["annot_snp_missingCol"]) < 1,
                             c(config["annot_snp_rsIDCol"], "chromosome", "position")]
names(snp.study) <- c("rsID", "chromosome", "position")
nrow(snp.study)

# remove missing from external
if (!is.na(config["ext_annot_snp_missingCol"])) {
  ext.sel <- getVariable(ext.snpAnnot, config["ext_annot_snp_missingCol"]) < 1
} else {
  ext.sel <- rep(TRUE, nrow(ext.snpAnnot))
}
table(ext.sel)
snp.ext <- pData(ext.snpAnnot)[ext.sel,
                               c(config["ext_annot_snp_rsIDCol"], "chromosome", "position")]
names(snp.ext) <- c("rsID", "chromosome", "position")
nrow(snp.ext)

# find snps with same rsID, chrom, and position
names(snp.ext) <- c("rsID", "chromosome", "position")
snp.common <- merge(snp.study, snp.ext)
nrow(snp.common)

# are there any scans to exclude?
if (!is.na(config["disc_scan_exclude_file"])) {
  scan.exclude1 <- getobj(config["disc_scan_exclude_file"])
  stopifnot(all(scan.exclude1 %in% getScanID(genoData)))
} else {
  scan.exclude1 <- NULL
}
length(scan.exclude1)

if (!is.na(config["ext_scan_exclude_file"])) {
  scan.exclude2 <- getobj(config["ext_scan_exclude_file"])
  stopifnot(all(scan.exclude2 %in% getScanID(extData)))
} else {
  scan.exclude2 <- NULL
}
length(scan.exclude2)

discord <- duplicateDiscordanceAcrossDatasets(genoData, extData,
  subjName.cols=c(config["annot_scan_subjectCol"], config["ext_annot_scan_subjectCol"]),
  snpName.cols=c(config["annot_snp_rsIDCol"], config["ext_annot_snp_rsIDCol"]),
  scan.exclude1=scan.exclude1, scan.exclude2=scan.exclude2,
  snp.include=snp.common$rsID)

save(discord, file=config["out_disc_file"])

# print results
nsubj <- length(discord$discordance.by.subject)
allpairs <- unlist(discord$discordance.by.subject)
npairs <- length(allpairs)
subjsum <- paste(npairs, "duplicate pairs among", nsubj, "subjects")
print(subjsum)

ndisc <- sum(discord$discordance.by.snp$discordant > 0)
print(paste(ndisc, "SNPs with > 0 discordant calls"))

# plot
disc <- allpairs[order(allpairs)]
rank <- 1:length(disc)
pdf(config["out_disc_plot"], width=6, height=6)
plot(disc, rank, xlab="discordance rate", ylab="rank", main=subjsum)
dev.off()

close(nc)
close(ext.nc)
