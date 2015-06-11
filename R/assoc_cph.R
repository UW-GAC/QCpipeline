##########
# Run coxph association test
# Usage: R --args config.file chromosome < assoc_cph.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# read chromosome
if (length(args) < 2) stop("missing chromosome")
(chromosome <- as.integer(args[2]))

## if (length(args) < 3) stop("missing segment")
## (segment <- as.integer(args[3]))

# check config and set defaults
required <- c("annot_scan_file", "covars", "covars_as_factor",
              "geno_file", "event", "time_to_event")
optional <- c("effect_allele", "gene_action", "ivar", "out_assoc_prefix",
              "scan_exclude", "block_size", "snp_segment_file")
default <- c("minor", "additive", NA, "assoc", NA, 5000, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# make genotypedata

scanAnnot <- getobj(config["annot_scan_file"])

data <- GenotypeReader(config["geno_file"])
sid <- getScanID(data)
scanAnnot <- scanAnnot[match(sid, scanAnnot$scanID),] 
stopifnot(all(scanAnnot$scanID==sid))

# set categorical variables in association models as factor
if (!is.na(config["covars_as_factor"])) {
    factors <- unlist(strsplit(config["covars_as_factor"], " ", fixed=TRUE))
    for (i in factors) {
        scanAnnot[[i]] <- as.factor(scanAnnot[[i]])
    }
}

genoData <- GenotypeData(data, scanAnnot=scanAnnot)
  
# event and time_to_event variables ("outcome")
(event <- config["event"])
(time_to_event <- config["time_to_event"])

# covariates
if (!is.na(config["covars"])) {
    covars <- unlist(strsplit(config["covars"], " ", fixed=TRUE))
    ## for Y when sex needs to be taken out of model
    if (chromosome == YchromCode(genoData)) {
        covars <- setdiff(covars, "sex")
    }
} else {
    warning("Are you sure you want to run without any covariates?")
    covars <- NULL
}
covars

# interaction variables for each model
if (!is.na(config["ivar"])) {
  ivar <- config["ivar"]
  # for Y when sex needs to be taken out of model
  if (ivar == "sex" & chromosome == YchromCode(genoData)) {
    ivar <- NULL
  }
} else {
  ivar <- NULL
}
ivar

# gene actions
(gene.action <- config["gene_action"])

# effect allele
(effect.allele <- config["effect_allele"])

# scan.exclude
if (!is.na(config["scan_exclude"])) {
  scan.exclude <- getobj(config["scan_exclude"])
} else {
  scan.exclude <- NULL
}


# snp selection
snpID <- getSnpID(genoData)
chrom <- getChromosome(genoData)
if (is.na(config["snp_segment_file"])) {
    ids <- snpID[chrom == chromosome]
} else {
    stop("snp_segment_file not supported yet")
}

## only run the association test if there are SNPs in this segment in the data.
## if there are none, quit successfully.
if (length(ids) == 0) {
  message("No SNPs in this segment. Exiting gracefully.")
  q(save="no", status=0)
}

(idr <- range(ids))
index <- which(getSnpID(genoData) %in% idr)

block.size <- as.integer(config["block_size"])

out <- assocCoxPH(genoData = genoData,
                  gene.action   = gene.action,
                  event         = event,
                  time.to.event = time_to_event,
                  covar         = covars,
                  ivar          = ivar,
                  strata        = NULL,
                  scan.exclude  = scan.exclude,
                  effectAllele  = effect.allele,
                  snpStart      = index[1],
                  snpEnd        = index[2],
                  block.size    = block.size,
                  verbose       = TRUE)

close(genoData)

#outfile <- paste0(config["out_assoc_prefix"], "_chr", chromosome, "_seg", segment, ".RData")
outfile <- paste0(config["out_assoc_prefix"], "_chr", chromosome, ".RData")
save(out, file=outfile)
