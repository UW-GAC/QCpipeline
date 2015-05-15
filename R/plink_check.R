##########
# Check PLINK file (assumes alleles coded A/B)
# Usage: R --args config.file <ABcoding> < plink_check.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "geno_file", "plink_prefix")
optional <- c("annot_scan_nameCol", "annot_snp_nameCol", "out_plink_logfile")
default <- c("Sample.Name", "rsID", "plink_check.log")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"])

# Need blank parent data for code to work
if (!hasVariable(scanAnnot, "father")) scanAnnot$father <- 0
if (!hasVariable(scanAnnot, "mother")) scanAnnot$mother <- 0

snpAnnot <- getobj(config["annot_snp_file"])

# check whether the PLINK file is coded A/B
if (length(args) > 1 & args[2] == "ABcoding") {
  message("Assuming PLINK file uses A/B coding")
  snpAnnot$alleleA <- "A"
  snpAnnot$alleleB <- "B"
}

data <- GenotypeReader(config["geno_file"])
genoData <- GenotypeData(data, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

res <- plinkCheck(genoData,
                  pedFile=config["plink_prefix"],
                  logFile=config["out_plink_logfile"],
                  family.col=config["annot_scan_nameCol"],
                  individual.col=config["annot_scan_nameCol"],
                  rs.col=config["annot_snp_nameCol"],
                  check.parents=FALSE)
close(genoData)

# indicate an error if check was not successful
stopifnot(res)
