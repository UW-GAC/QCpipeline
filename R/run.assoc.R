##########
# Run association tests
# Usage: R --args config.file chromosome < run.assoc.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "covar.list", "covars_as_factor",
              "gene_action", "model_type", "geno_file", "outcome")
optional <- c("ivar.list", "out_assoc_prefix",
              "scan_chrom_filter", "scan_exclude", "robust")
default <- c(NA, "assoc", NA, NA, FALSE)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# read chromosome
if (length(args) < 2) stop("missing chromosome")
# set chromosome numbers given in args on the command line
chromosome.set <- as.numeric(args[2])  
chromosome.set

# make genotypedata
scanAnnot <- getobj(config["annot_scan_file"])

data <- GenotypeReader(config["geno_file"])
sid <- getScanID(data)
scanAnnot <- scanAnnot[scanAnnot$scanID %in% sid,] 
stopifnot(all(scanAnnot$scanID==sid))


# set categorical variables in association models as factor 
factors <- unlist(strsplit(config["covars_as_factor"], " ", fixed=TRUE))
idx <- which(varLabels(scanAnnot) %in% factors)
for (i in idx) {
  pData(scanAnnot)[,i] <- as.factor(pData(scanAnnot)[,i])
}
  
# sample/subject chromosome.filter sorted by scanID
if (!is.na(config["scan_chrom_filter"])) {
  filt <- getobj(config["scan_chrom_filter"])
  stopifnot(all(rownames(filt) == scanAnnot$scanID))
  chk <- apply(filt, 1, function(x) !all(x==F))
  print(table(chk)) # TRUE: to be kept for association
} else {
  filt <- NULL
}

# chromosomes to run analysis on
# outcome variables
outcome <- unlist(strsplit(config["outcome"]," "))
outcome

# covariates for each model
covar.list <- getobj(config["covar.list"])
# for Y when sex needs to be taken out of model
if (chromosome.set == 25) {
  idx <- sapply(covar.list, function(x) which(x %in% "sex"))
  for (i in 1:length(idx)) {
    if (length(idx[[i]]) > 0) {
      covar.list[[i]] <- covar.list[[i]][-idx[[i]]]
      if (length(covar.list[[i]]) == 0) covar.list[[i]] <- ""
    }
  }
}
covar.list

# interaction variables for each model
if (!is.na(config["ivar.list"])) {
  ivar.list <- getobj(config["ivar.list"])
  # for Y when sex needs to be taken out of model
  if (chromosome.set == 25) {
    idx <- sapply(ivar.list, function(x) which(x %in% "sex"))
    for (i in 1:length(idx)) {
      if (length(idx[[i]]) > 0) {
        ivar.list[[i]] <- ivar.list[[i]][-idx[[i]]]
        if (length(ivar.list[[i]]) == 0) ivar.list[[i]] <- ""
      }
    }
  }
} else {
  ivar.list <- NULL
}
ivar.list

# model types
model.type <- unlist(strsplit(config["model_type"]," "))
stopifnot(all(model.type %in% c("logistic", "linear", "Logistic", "Linear")))
model.type

# gene actions
gene.action <- unlist(strsplit(config["gene_action"]," "))
gene.action.list <- list()
for (i in 1:length(gene.action)) {
  gene.action.list [[i]] <- gene.action[i]
}
gene.action.list

# check lengths
len <- length(outcome)
len
stopifnot(length(outcome) == len & length(model.type) == len &
          length(gene.action.list) == len & length(covar.list) == len )

# output
outfile <- config["out_assoc_prefix"]

# scan.exclude
scan.exclude <- NULL
if (!is.na(config["scan_exclude"])) {
  scan.exclude <- getobj(config["scan_exclude"])
}

genoData <- GenotypeData(data, scanAnnot=scanAnnot)

assocTestRegression(genoData,
                    outcome = outcome, 
                    covar.list = covar.list, 
                    ivar.list = ivar.list,
                    model.type = model.type,
                    scan.exclude = scan.exclude,
                    gene.action.list = gene.action.list,
                    scan.chromosome.filter = filt,
                    chromosome.set = chromosome.set, 
                    outfile = outfile,
                    robust = as.logical(config["robust"]))
