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
required <- c("annot_scan_file", "annot_snp_alleleACol", "annot_snp_alleleBCol",
              "annot_snp_file", "nc_samp_geno_file", "out_plink_prefix", "out_log_prefix")
optional <- c("annot_scan_subjectCol", "annot_snp_rsIDCol")
default <- c("subjectID", "rsID")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"])
snpAnnot <- getobj(config["annot_snp_file"])
nc <- NcdfGenotypeReader(config["nc_samp_geno_file"])
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

scan.exclude <- scanAnnot$scanID[!scanAnnot$subj.plink]
length(scan.exclude)

ped <- config["out_plink_prefix"]
plinkWrite(genoData, pedFile=ped,
           individual.col=config["annot_scan_subjectCol"],
           alleleA.col=config["annot_snp_alleleACol"],
           alleleB.col=config["annot_snp_alleleBCol"],
           rs.col=config["annot_snp_rsIDCol"],
           scan.exclude=scan.exclude)

log <- paste(config["out_log_prefix"], ".log", sep="")
res <- plinkCheck(genoData, pedFile=ped, logFile=log,
           individual.col=config["annot_scan_subjectCol"],
           alleleA.col=config["annot_snp_alleleACol"],
           alleleB.col=config["annot_snp_alleleBCol"],
           rs.col=config["annot_snp_rsIDCol"],
           scan.exclude=scan.exclude)
# indicate an error if check was not successful
stopifnot(res)
