##########
# Combine two netCDF/GDS files into a single GDS file
# Usage: R --args config.file < combine_gds.R
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
              "ext_annot_snp_file", "ext_geno_file", "gds_geno_file")
optional <- c("annot_scan_subjectCol", "annot_snp_rsIDCol",
              "annot_snp_alleleACol", "annot_snp_alleleBCol", 
              "comb_scan_exclude_file", "comb_snp_exclude_file",
              "ext_annot_scan_subjectCol", "ext_annot_snp_rsIDCol",
              "ext_annot_snp_alleleACol", "ext_annot_snp_alleleBCol", 
              "out_comb_prefix", "remove_discordant")
default <- c("subjectID", "rsID", "alleleA", "alleleB", NA, NA, "subjectID", "rsID",
             "alleleA", "alleleB", "comb_geno", TRUE)
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"]); dim(scanAnnot)
# this might be a subject-level file, so subset annotation
geno <- GenotypeReader(config["gds_geno_file"])
scanAnnot <- scanAnnot[match(getScanID(geno), getScanID(scanAnnot)),]; dim(scanAnnot)
snpAnnot <- getobj(config["annot_snp_file"])
# remake snpAnnot object with alleles A and B
snpAnnot <- SnpAnnotationDataFrame(pData(snpAnnot),
                                   alleleACol=config["annot_snp_alleleACol"],
                                   alleleBCol=config["annot_snp_alleleBCol"])
genoData <- GenotypeData(geno, snpAnnot=snpAnnot, scanAnnot=scanAnnot)

ext.scanAnnot <- getobj(config["ext_annot_scan_file"]); dim(ext.scanAnnot)
# this might be a subject-level file, so subset annotation
ext.geno <- GenotypeReader(config["ext_geno_file"])
ext.scanAnnot <- ext.scanAnnot[match(getScanID(ext.geno), getScanID(ext.scanAnnot)),]; dim(ext.scanAnnot)
ext.snpAnnot <- getobj(config["ext_annot_snp_file"])
ext.snpAnnot <- SnpAnnotationDataFrame(pData(ext.snpAnnot),
                                       alleleACol=config["ext_annot_snp_alleleACol"],
                                       alleleBCol=config["ext_annot_snp_alleleBCol"])
ext.genoData <- GenotypeData(ext.geno, snpAnnot=ext.snpAnnot, scanAnnot=ext.scanAnnot)

genoDataList <- list("study"=genoData, "ext"=ext.genoData)

dupids <- intersect(scanAnnot$scanID, ext.scanAnnot$scanID)
if (length(dupids) > 0) {
    stop("scan IDs are not unique - files cannot be combined")
}

# any scans to exclude?
if (!is.na(config["comb_scan_exclude_file"])) {
    scan.exclude <- getobj(config["comb_scan_exclude_file"])
    sampleList <- list("study"=setdiff(scanAnnot$scanID, scan.exclude),
                       "ext"=setdiff(ext.scanAnnot$scanID, scan.exclude))
} else {
    sampleList <- list("study"=scanAnnot$scanID,
                       "ext"=ext.scanAnnot$scanID)
}
lapply(sampleList, length)

# remove discordant snps?
if (as.logical(config["remove_discordant"])) {
    disc <- duplicateDiscordanceAcrossDatasets(genoData, ext.genoData,
                                               match.snps.on=c("name", "position", "alleles"),
                                               subjName.cols=c(config["annot_scan_subjectCol"],
                                                   config["ext_annot_scan_subjectCol"]),
                                               snpName.cols=c(config["annot_snp_rsIDCol"],
                                                   config["ext_annot_snp_rsIDCol"]),
                                               one.pair.per.subj=FALSE,
                                               scan.exclude1=setdiff(scanAnnot$scanID, sampleList[["study"]]),
                                               scan.exclude2=setdiff(ext.scanAnnot$scanID, sampleList[["ext"]]))
    nsubj <- length(disc$discordance.by.subject)
    allpairs <- unlist(disc$discordance.by.subject)
    npairs <- length(allpairs)
    print(paste(npairs, "duplicate pairs among", nsubj, "subjects"))
    print(summary(allpairs))

    ndisc <- sum(disc$discordance.by.snp$discordant > 0)
    print(paste(ndisc, "SNPs with > 0 discordant calls"))

    disc0 <- disc$discordance.by.snp$discordant == 0
    snpList <- list("study"=disc$discordance.by.snp$snpID1[disc0],
                    "ext"=disc$discordance.by.snp$snpID2[disc0])
} else {
    snpList <- list("study"=snpAnnot$snpID,
                    "ext"=ext.snpAnnot$snpID)
}

# any snps to exclude?
## TO-DO: allow exclusion of external SNPs?
if (!is.na(config["comb_snp_exclude_file"])) {
  snp.exclude <- getobj(config["comb_snp_exclude_file"])
  snpList[["study"]] <- setdiff(snpList[["study"]], snp.exclude)
}
lapply(snpList, length)


# snp names to match
snpNameList <- list("study"=config["annot_snp_rsIDCol"],
                    "ext"=config["ext_annot_snp_rsIDCol"])


# make gds files
## TO-DO - allow merging on just position or just name???
gdsMerge(genoDataList, sampleList=sampleList, snpList=snpList,
         match.snps.on=c("position", "name"), snpNameList=snpNameList,
         sortByScanID=FALSE, newSnpID=FALSE, outPrefix=config["out_comb_prefix"])

# sample information
samp.info1 <- getVariable(scanAnnot, c("scanID", config["annot_scan_subjectCol"], "sex"))
names(samp.info1)[2] <- "subjectID"
samp.info1$source <- "study"
samp.info2 <- getVariable(ext.scanAnnot, c("scanID", config["ext_annot_scan_subjectCol"], "sex"))
names(samp.info2)[2] <- "subjectID"
samp.info2$source <- "ext"
samp.info <- rbind(samp.info1, samp.info2)

scanfile <- paste0(config["out_comb_prefix"], "_scanAnnot.RData")
comb.scanAnnot <- getobj(scanfile)
samp.info <- samp.info[match(comb.scanAnnot$scanID, samp.info$scanID),]
pData(comb.scanAnnot) <- samp.info
save(comb.scanAnnot, file=scanfile)

# check
gdsMergeCheck(genoDataList, outPrefix=config["out_comb_prefix"])
