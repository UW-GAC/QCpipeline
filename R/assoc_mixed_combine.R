##########
# Combine assoc files for a chromosome
# Usage: R --args analysis_id dbname chromosome < combine_assoc.R
##########

library(GWASTools)
library(QCpipeline)
library(GWASbyChr)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
print(args)

date()

if (length(args) < 1) stop("missing config file")
configFile <- args[1]

if (length(args) < 2) stop("missing chromosome")
chromosome <- as.numeric(args[2])

## connect to DB
config <- readConfig(configFile)
required <- c("geno_file")
opt <- c("snp_filter_column"=NA,
         "out_assoc_prefix"="assoc_mixed/study",
         "maf.filter.type"="snp.specific",
         "maf.absolute.threshold"=0.02,
         "maf.linear.threshold"=30,
         "maf.logistic.threshold"=50,
         "annot_snp_file"=NA)
optional <- names(opt)
default <- unname(opt)
config <- setConfigDefaults(config, required, optional, default)
print(config)


if (file_test("-d", config["geno_file"])){
  genoByChr <- GenotypeDataByChr(config["geno_file"])
  snpAnnot <- getSnpAnnotation(genoByChr, chromosome=chromosome)
} else {
  snpAnnot <- getobj(config["annot_snp_file"])
}

# required columns: rsID, alleleA, alleleB
requiredColumns <- c("position", "alleleA", "alleleB")
possibleColumns <- unique(c("info", "type", "oevar", "composite.filter"))
if (!is.na(config["snp_filter_column"])) possibleColumns <- unique(c(possibleColumns, config["snp_filter_column"]))


# make sure all required columns are in the snp annotation
stopifnot(all(requiredColumns %in% varLabels(snpAnnot)))

# select only the possible columns that are in the snp annotation
possibleColumns <- possibleColumns[possibleColumns %in% varLabels(snpAnnot)]

# get list of segments for that chromosome
# may have multiple chromosomes in one snp annotation
segments <- sort(unique(snpAnnot$segment[snpAnnot$chromosome == chromosome]))

seg.files <- paste(config["out_assoc_prefix"], "_chr-", chromosome, "_seg", segments, ".RData", sep="")
stopifnot(all(file.exists(seg.files)))

assoc.list <- lapply(seg.files, getobj)
assoc <- do.call(rbind, assoc.list)
rm(assoc.list)
assoc <- assoc[order(assoc[,"snpID"]),]

# update names
# names(assoc)[names(assoc) == "Wald.Stat"] <- "Stat"
# names(assoc)[names(assoc) == "Wald.pval"] <- "pval"
# names(assoc)[names(assoc) == "Score.pval"] <- "pval"
# names(assoc)[names(assoc) == "Score.Stat"] <- "Stat"

#names(assoc)[names(assoc) == "Est"] <- "Beta"
names(assoc) <- gsub("Est", "Beta", names(assoc))
names(assoc)[names(assoc) == "chr"] <- "chromosome"

## merge assoc results with snp annotation
j <- match(assoc[,"snpID"], snpAnnot$snpID)
stopifnot(all(!is.na(j)))
snp <- pData(snpAnnot)[j, unique(c(requiredColumns, possibleColumns))]
assoc <- cbind(assoc, snp)

## character chromosome names
assoc$chromosome <- chromToChar(assoc$chromosome)

# add qual filter info
if (!is.na(config["snp_filter_column"])){
  names(assoc)[names(assoc) %in% config["snp_filter_column"]] <- "composite.filter"
} else {
  # everything passes
  assoc$composite.filter <- TRUE
}

# add comp.maf.filt

## MAF filter
maf.thresh <- switch(config["model_type"],
                     linear=as.numeric(config["maf.linear.threshold"]),
                     logistic=as.numeric(config["maf.logistic.threshold"]))


maf.filt <- switch(config["maf.filter.type"],
                              snp.specific=assoc$effN > maf.thresh,
                              absolute=assoc$MAF > as.numeric(config["maf.absolute.threshold"]))

assoc$comp.maf.filter <- assoc$composite.filter & maf.filt

# save the output
outfile <- paste(config["out_assoc_prefix"], "_chr-", chromosome, ".RData", sep="")
save(assoc, file=outfile)

# remove the intermediate segments
file.remove(seg.files)

