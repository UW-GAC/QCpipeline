library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
print(args)

date()

if (length(args) < 1) stop("missing config")
config <- readConfig(args[1])

required <- c("imputed_gds_dir", "imputed_gds_prefix",
              "out_gds_dir", "out_gds_prefix")
optional <- c("type_column",
              "type_include")
default <- c("type",
             "2 3")
config <- setConfigDefaults(config, required, optional, default)
print(config)

if (length(args) < 2) stop("missing chromosome")
(chromosome <- args[2])

types <- as.integer(strsplit(config["type_include"], " ")[[1]])


olgaData <- OlgaGenotypeData(directory=config["imputed_gds_dir"])
stopifnot(chromosome %in% getValidChromosomes(olgaData))

snpAnnot <- getSnpAnnotation(olgaData, chromosome)
snp <- getVariable(snpAnnot, c("snpID", config["type_column"]))
names(snp) <- c("snpID", "type")
snp.include <- snpAnnot$snpID[which(snp$type %in% types)]
length(snp.include)

gds.imputed.file <- paste(config["imputed_gds_dir"], "/", config["imputed_gds_prefix"], "_chr-", chromosome, ".gds", sep="")
stopifnot(file.exists(gds.imputed.file))

gds.observed.file <- paste(config["out_gds_dir"], "/", config["out_gds_prefix"], "_chr-", chromosome, ".gds", sep="")

gds.observed.file.tmp <- tempfile()
gdsSubset(parent.gds=gds.imputed.file, sub.gds=gds.observed.file.tmp, snp.include=snp.include, sub.storage="bit2", zipflag="ZIP.max", verbose=TRUE)

# add compression -- no need here, plus it gives a stream write error...
#gds <- openfn.gds(gds.observed.file.tmp, readonly=F)
#genoNode <- index.gdsn(gds, "genotype")
#compression.gdsn(genoNode, compress="ZIP.max")
#closefn.gds(gds)

cleanup.gds(gds.observed.file.tmp)

stopifnot(file.copy(gds.observed.file.tmp, gds.observed.file, overwrite=T))
file.remove(gds.observed.file.tmp)

# check it
gdsSubsetCheck(parent.gds=gds.imputed.file, sub.gds=gds.observed.file, snp.include=snp.include, sub.storage="bit2", verbose=TRUE)

# because we know these are genotypes, double-check with GdsGenotypeReader, sample by sample
#gds.imputed <- GdsGenotypeReader(gds.imputed.file)
#for (i in 1:nscan(gds.imputed)){
#  
#}

# make a new snp annotation for it
snpAnnot.observed <- snpAnnot[snpAnnot$snpID %in% snp.include, ]

# check the new snp annotation against the file
gds.observed <- GdsGenotypeReader(gds.observed.file)
genoData <- GenotypeData(gds.observed, snpAnnot.observed)
close(genoData)

# save the snp annotation
snp.observed.file <- paste(config["out_gds_dir"], "/", config["out_gds_prefix"], "_chr-", chromosome, "_snpAnnot.RData", sep="")
save(snpAnnot.observed, file=snp.observed.file)

