##########
# Create filtered PLINK file
# Usage: R --args config.file < plink_unfiltered.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file",
              "annot_snp_file", "samp_geno_file", "out_plink_prefix", "out_log_prefix")
optional <- c("annot_scan_subjectCol", "annot_snp_alleleACol", "annot_snp_alleleBCol",
              "annot_snp_rsIDCol")
default <- c("subjectID", "alleleA", "alleleB", "rsID")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"])
snpAnnot <- getobj(config["annot_snp_file"])
# remake snpAnnot object with alleles A and B
snpAnnot <- SnpAnnotationDataFrame(pData(snpAnnot),
                                   alleleACol=config["annot_snp_alleleACol"],
                                   alleleBCol=config["annot_snp_alleleBCol"])
data <- GenotypeReader(config["samp_geno_file"])
genoData <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

scan.exclude <- scanAnnot$scanID[!scanAnnot$subj.plink]
length(scan.exclude)

ped <- config["out_plink_prefix"]
plinkWrite(genoData, pedFile=ped,
           individual.col=config["annot_scan_subjectCol"], rs.col=config["annot_snp_rsIDCol"],
           scan.exclude=scan.exclude)

log <- paste(config["out_log_prefix"], ".log", sep="")
res <- plinkCheck(genoData, pedFile=ped, logFile=log,
           individual.col=config["annot_scan_subjectCol"], rs.col=config["annot_snp_rsIDCol"],
           scan.exclude=scan.exclude)
# indicate an error if check was not successful
stopifnot(res)
