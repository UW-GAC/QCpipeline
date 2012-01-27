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
# find snps with same rsID, chrom, and position
snp.ext <- pData(ext.snpAnnot)[,c(config["ext_annot_snp_rsIDCol"], "chromosome", "position")]
names(snp.ext) <- c("rsID", "chromosome", "position")
snp.common <- merge(snp.study, snp.ext)
nrow(snp.common)

discord <- duplicateDiscordanceAcrossDatasets(genoData, extData,
  subjName.cols=c(config["annot_scan_subjectCol"], config["ext_annot_scan_subjectCol"]),
  snpName.cols=c(config["annot_snp_rsIDCol"], config["ext_annot_snp_rsIDCol"]),
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
