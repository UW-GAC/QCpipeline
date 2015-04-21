##########
# Combine chromosome segments for association tests
# Usage: R --args config.file start end < assoc_combine.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "model_type")
optional <- c("annot_snp_filtCol", "annot_snp_rsIDCol", "maf.filter.type",
              "maf.absolute.threshold", "maf.linear.threshold", "maf.logistic.threshold",
              "out_assoc_prefix")
default <- c("quality.filter", "rsID", "snp.specific", 0.02, 30, 50, "assoc")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# chromosomes to merge
if (length(args) < 3) stop("missing start and end chromosomes")
start <- as.integer(args[2])
end <- as.integer(args[3])

pathprefix <- config["out_assoc_prefix"]
tmp <- list()
for (j in start:end) {
    fname <- paste0(pathprefix, "_chr", j, ".RData")
    if (file.exists(fname)) {
        tmp[[as.character(j)]] <- getobj(fname)
    }
}
combined <- do.call("rbind", tmp)
rm(tmp)

## add rsID, chromosome, filters
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
stopifnot(all(combined$snpID %in% snpID))
index <- match(combined$snpID,snpID)

combined$rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"], index=index)
##combined$chromosome <- getChromosome(snpAnnot, char=TRUE, index=index)
combined$chromosome <- combined$chr
combined$chromosome[combined$chromosome == 23] <- "X"
combined$chromosome[combined$chromosome == 24] <- "XY"
combined$chromosome[combined$chromosome == 25] <- "Y"
combined$chromosome[combined$chromosome == 26] <- "M"
combined$chromosome[combined$chromosome == 27] <- "U"
combined$chromosome[combined$chromosome == 28] <- "XYY"
combined$composite.filter <- getVariable(snpAnnot, config["annot_snp_filtCol"], index=index)

## MAF filter
maf.thresh <- switch(config["model_type"],
                     linear=as.numeric(config["maf.linear.threshold"]),
                     logistic=as.numeric(config["maf.logistic.threshold"]))

## snp-specific filtering for logistic is calculated the number of cases or controls, whichever is smaller
## for linear it is all samples
N <- switch(config["model_type"],
            linear=combined$n,
            logistic=pmin(combined$n0, combined$n1))

maf.filt <- switch(config["maf.filter.type"],
                   absolute=(!is.na(combined$MAF) &
                             combined$MAF > as.numeric(config["maf.absolute.threshold"])),
                   snp.specific=(!is.na(combined$MAF) &
                                 2 * combined$MAF * (1 - combined$MAF) * N > maf.thresh))
combined$comp.maf.filter <- combined$composite.filter & maf.filt

print(dim(combined))
fname <- paste0(pathprefix, "_combined_qual_filt.RData")
save(combined, file=fname)
