##########
# Duplicate discordance across two datasets
# Usage: R --args config.file type < dup_disc_2sets.R
# types are "all", "minor", "miss1fail", "miss2fail", "senspec"
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file_1", "annot_scan_file_2",
              "annot_snp_file_1", "annot_snp_file_2",
              "geno_file_1", "geno_file_2",
              "out_prefix")
optional <- c("annot_scan_subjCol_1", "annot_scan_subjCol_2",
              "annot_snp_snpCol_1", "annot_snp_snpCol_2",
              "scan_exclude_file_1", "scan_exclude_file_2",
              "snp_include_file")
default <- c("subjectID", "subjectID",
             "rsID", "rsID",
             NA, NA, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

if (length(args) < 2) stop("missing discordance type")
type <- args[2]

if (type == "all") {
  func.name <- "duplicateDiscordanceAcrossDatasets"
  MAonly <- FALSE
  miss.fail <- c(FALSE, FALSE)
  outfile <- paste(config["out_prefix"], "_all.RData", sep="")
  message("calculating overall discordance")
} else if (type == "minor") {
  func.name <- "duplicateDiscordanceAcrossDatasets"
  MAonly <- TRUE
  miss.fail <- c(FALSE, FALSE)
  outfile <- paste(config["out_prefix"], "_minor.RData", sep="")
  message("calculating minor allele discordance")
} else if (type == "miss1fail") {
  func.name <- "duplicateDiscordanceAcrossDatasets"
  MAonly <- TRUE
  miss.fail <- c(TRUE, FALSE)
  outfile <- paste(config["out_prefix"], "_miss1fail.RData", sep="")
  message("calculating minor allele discordance, including missing in dataset 1")
} else if (type == "miss2fail") {
  func.name <- "duplicateDiscordanceAcrossDatasets"
  MAonly <- TRUE
  miss.fail <- c(FALSE, TRUE)
  outfile <- paste(config["out_prefix"], "_miss2fail.RData", sep="")
  message("calculating minor allele discordance, including missing in dataset 2")
} else if (type == "senspec") {
  func.name <- "minorAlleleSensitivitySpecificity"
  outfile <- paste(config["out_prefix"], "_senspec.RData", sep="")
  message("calculating minor allele sensitivity and specificity")
} else {
  stop("unknown discordance type; options are all, minor, miss1fail, miss2fail, senspec")
}

# dataset 1
scanAnnot <- getobj(config["annot_scan_file_1"])
snpAnnot <- getobj(config["annot_snp_file_1"])
data <- GenotypeReader(config["geno_file_1"])
genoData1 <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
genoData1

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file_1"])) {
  scan.exclude1 <- getobj(config["scan_exclude_file_1"])
  stopifnot(all(scan.exclude1 %in% getScanID(genoData1)))
} else {
  scan.exclude1 <- NULL
}
length(scan.exclude1)

  
# dataset 2
scanAnnot <- getobj(config["annot_scan_file_2"])
scanID <- getScanID(scanAnnot)

snpAnnot <- getobj(config["annot_snp_file_2"])
snpID <- getSnpID(snpAnnot)

data <- GenotypeReader(config["geno_file_2"])
genoData2 <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
genoData2

# are there any scans to exclude?
if (!is.na(config["scan_exclude_file_2"])) {
  scan.exclude2 <- getobj(config["scan_exclude_file_2"])
  stopifnot(all(scan.exclude2 %in% getScanID(genoData2)))
} else {
  scan.exclude2 <- NULL
}
length(scan.exclude2)


# snps
if (!is.na(config["snp_include_file"])) {
  snp.include <- getobj(config["snp_include_file"])
} else {
  snp.include <- NULL
}
length(snp.include)

if (func.name == "duplicateDiscordanceAcrossDatasets") {
  disc <- duplicateDiscordanceAcrossDatasets(genoData1, genoData2,
    subjName.cols=config[c("annot_scan_subjCol_1","annot_scan_subjCol_2")],
    snpName.cols=config[c("annot_snp_snpCol_1", "annot_snp_snpCol_2")],
    one.pair.per.subj=TRUE, minor.allele.only=MAonly,
    missing.fail=miss.fail,
    scan.exclude1=scan.exclude1, scan.exclude2=scan.exclude2,
    snp.include=snp.include)
} else if (func.name == "minorAlleleSensitivitySpecificity") {
  disc <- minorAlleleDetectionAccuracy(genoData1, genoData2,
    subjName.cols=config[c("annot_scan_subjCol_1","annot_scan_subjCol_2")],
    snpName.cols=config[c("annot_snp_snpCol_1", "annot_snp_snpCol_2")],
    scan.exclude1=scan.exclude1, scan.exclude2=scan.exclude2,
    snp.include=snp.include)
} else {
  stop("unknown function name")
}

save(disc, file=outfile)

close(genoData1)
close(genoData2)
