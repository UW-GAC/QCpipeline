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
required <- c("annot_scan_file", "covar.list", "covars_as_factor", "gene_action", "model_type", "nc_geno_file", "outcome")
optional <- c("assoc_output", "chrom_filter", "scan_exclude")
default <- c("assoc", NA, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# read chromosome
if (length(args) < 2) stop("missing chromosome")
chromosome.set <- as.numeric(args[2])  # this sets chromosome numbers given in args on the command line
chromosome.set

# make genotypedata
scanAnnot <- getobj(config["annot_scan_file"])

nc <- NcdfGenotypeReader(config["nc_geno_file"])
sid <- getScanID(nc)
scanAnnot <- scanAnnot[scanAnnot$scanID %in% sid,] 
stopifnot(all(scanAnnot$scanID==sid))


# set categorical variables in association models as factor 
factors <- config["covars_as_factor"]
factors <- unlist(strsplit(factors," ", fixed=TRUE))
idx <- which(varLabels(scanAnnot) %in% factors)
for (i in idx)
{
  pData(scanAnnot)[,i] <- as.factor(pData(scanAnnot)[,i])
}
  
# sample/subject chromosome.filter sorted by scanID
if (!is.na(config["chrom_filter"])) {
  filt <- getobj(config["chrom_filter"])
  stopifnot(all(rownames(filt) == scanAnnot$scanID))
  chk <- apply(filt, 1, function(x) !all(x==F))
  print(table(chk)) # TRUE: to be kept for association
} else {
  filt <- NULL
}

# chromosomes to run analysis on
# outcome variables
outcome <- config["outcome"]
outcome <- unlist(strsplit(outcome," "))
outcome

# covariates for each model
covar.list <- getobj(config["covar.list"])
# for Y when sex needs to be taken out of model
if (chromosome.set == 25)
{
  idx <- sapply(covar.list, function(x) which(x %in% "sex"))
  for (i in 1:length(idx))
  {
    if (length(idx[[i]]) > 0) {
      covar.list[[i]] <- covar.list[[i]][-idx[[i]]]
      if (length(covar.list[[i]]) == 0) covar.list[[i]] <- ""
    }
  }
}
covar.list

# model types
model.type <- config["model_type"]
model.type <- unlist(strsplit(model.type," "))
stopifnot(all(model.type %in% c("logistic", "linear", "Logistic", "Linear")))
model.type

# gene actions
gene.action <- config["gene_action"]
gene.action <- unlist(strsplit(gene.action," "))
gene.action.list <- list()
for (i in 1:length(gene.action))
{
  gene.action.list [[i]] <- gene.action[i]
}
gene.action.list

# check lengths
len <- length(outcome)
len
stopifnot(length(outcome) == len & length(model.type) == len &
          length(gene.action.list) == len & length(covar.list) == len )

# output
outfile <- config["assoc_output"]

# scan.exclude
scan.exclude <- NULL
if (!is.na(config["scan_exclude"]))
{
  scan.exclude <- getobj(config["scan_exclude"])
}

genoData <- GenotypeData(nc,scanAnnot=scanAnnot)

assocTestRegression(			 genoData,
					 outcome = outcome, 
					 covar.list = covar.list, 
					 model.type = model.type,
                                         scan.exclude = scan.exclude,
                                         gene.action.list = gene.action.list,
					 scan.chromosome.filter = filt,
					 chromosome.set = chromosome.set, 
                                         outfile = outfile
)
