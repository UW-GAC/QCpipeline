##########
# Create PLINK file
# Usage: R --args config.file < plink_from_geno.R
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
              "annot_snp_file", "geno_file", "out_plink_prefix", "out_log_prefix")
optional <- c("annot_scan_subjectCol", "annot_snp_alleleACol", "annot_snp_alleleBCol",
              "annot_snp_rsIDCol")
default <- c("subjectID", "alleleA", "alleleB", "rsID")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"])
data <- GenotypeReader(config["geno_file"])
# take subset to match gds file
scanAnnot <- scanAnnot[match(getScanID(data), getScanID(scanAnnot)), ]
snpAnnot <- getobj(config["annot_snp_file"])
# remake snpAnnot object with alleles A and B
snpAnnot <- SnpAnnotationDataFrame(pData(snpAnnot),
                                   alleleACol=config["annot_snp_alleleACol"],
                                   alleleBCol=config["annot_snp_alleleBCol"])

genoData <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# exclude scans that aren't subj.plink
scan_exclude <- scanAnnot$scanID[!scanAnnot$subj.plink]
length(scan_exclude)

ped <- config["out_plink_prefix"]
plinkWrite(genoData, pedFile=ped,
           individual.col=config["annot_scan_subjectCol"],
           rs.col=config["annot_snp_rsIDCol"],
           scan.exclude=scan_exclude)

log <- paste(config["out_log_prefix"], ".log", sep="")
res <- plinkCheck(genoData, pedFile=ped, logFile=log,
           individual.col=config["annot_scan_subjectCol"],
           rs.col=config["annot_snp_rsIDCol"],
           scan.exclude=scan_exclude)
# indicate an error if check was not successful
stopifnot(res)
