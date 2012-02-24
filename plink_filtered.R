##########
# Create filtered PLINK file
# Usage: R --args config.file < plink_filtered.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

scanAnnot <- getobj(config["annot_scan_file"])
scanAnnot <- scanAnnot[scanAnnot$subj.plink,]
snpAnnot <- getobj(config["annot_snp_file"])
nc <- NcdfGenotypeReader(config["nc_subj_geno_file"])
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

ped <- paste(config["out_plink_prefix"], "_filtered", sep="")
plinkWrite(genoData, pedFile=ped,
           individual.col=config["annot_scan_subjectCol"],
           alleleA.col=config["annot_snp_alleleACol"],
           alleleB.col=config["annot_snp_alleleBCol"],
           rs.col=config["annot_snp_rsIDCol"])

log <- paste(config["out_log_prefix"], "_filtered.log", sep="")
res <- plinkCheck(genoData, pedFile=ped, logFile=log,
           individual.col=config["annot_scan_subjectCol"],
           alleleA.col=config["annot_snp_alleleACol"],
           alleleB.col=config["annot_snp_alleleBCol"],
           rs.col=config["annot_snp_rsIDCol"])
# indicate an error if check was not successful
stopifnot(res)
