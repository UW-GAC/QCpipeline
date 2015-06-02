##########
# Estimate variance component for association tests
# Usage: R --args analysis_id dbname < run_assoc_estVC.R
##########

library(GWASTools)
library(QCpipeline)
library(GWASbyChr)
library(OLGAanalysis)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
print(args)

date()

if (length(args) < 1) stop("missing config file")
(configFile <- args[1])


config <- readConfig(configFile)
## check config and set defaults
required <- c("annot_scan_file", "covars", "covars_as_factor",
              "model_type", "geno_file", "outcome",
              "covar_matrix_file",
              "scan_include_file")
opt <- c("out_assoc_prefix"="assoc_mixed/study",
         "out_plot_prefix"="assoc_mixed/study",
         "block_size"=5000)
optional <- names(opt)
default <- unname(opt)
config <- setConfigDefaults(config, required, optional, default)
print(config)


## scan annotation
scanAnnot <- getobj(config["annot_scan_file"])

# make geno data
genoByChr <- GenotypeDataByChr(config["geno_file"])
chr <- getValidChromosomes(genoByChr)[1]
genoData <- getGenoData(genoByChr, chromosome=chr, snpAnnot=FALSE, scanAnnot=FALSE)
scanAnnot <- scanAnnot[match(getScanID(genoData), scanAnnot$scanID)]
close(genoData)
# this is just for checking - will be needed later
genoData <- getGenoData(genoByChr, chromosome=chr, snpAnnot=TRUE, scanAnnot=scanAnnot)
close(genoData)

# set categorical variables in association models as factor
if (!is.na(config["covars_as_factor"])) {
  factors <- unlist(strsplit(config["covars_as_factor"], " ", fixed=TRUE))
  for (i in factors) {
    scanAnnot[[i]] <- as.factor(scanAnnot[[i]])
  }
}


# outcome variable
(outcome <- config["outcome"])

# covariates
if (!is.na(config["covars"])) {
  covars <- unlist(strsplit(config["covars"], " ", fixed=TRUE))
} else {
  warning("Are you sure you want to run without any covariates?")
  covars <- NULL
}
covars

# model type
(model.type <- config["model_type"])
if (model.type == "logistic"){
  model = "glmm"
} else if (model.type == "linear") {
  model = "mixed"
} else {
  stop("'model_type' must be 'logsitic' or 'linear'")
}
model


# scan.include
if (!is.na(config["scan_include_file"])) {
  scan.include <- getobj(config["scan_include_file"])
  scan.exclude <- setdiff(scanAnnot$scanID, scan.include)
} else {
  scan.include <- NULL
  scan.exclude <- NULL
}


# covariance matrices
cov.mat <- getobj(config["covar_matrix_file"])

# compute variance component estimates
res <- fitNullMixedModel(scanAnnot = scanAnnot,
                 covMatList = cov.mat,
                 outcome = outcome,
                 covar.vec = covars,
                 scan.exclude = scan.exclude,
                 model=model)

res$VC$converged

# add workingY and save
scanAnnot$workingY <- res$VC$workingY[match(scanAnnot$scanID, rownames(res$VC$cholSigmaInv))]
outfile <- paste0(config["out_assoc_prefix"], "_scanAnnot.RData")
save(scanAnnot, file=outfile)

# save variance components
VC <- res$VC
outfile <- paste0(config["out_assoc_prefix"], "_VC.RData")
save(VC, file=outfile)

# save null model
mod <- res$mod
outfile <- paste0(config["out_assoc_prefix"], "_nullmod.RData")
save(mod, file=outfile)

if (!res$VC$converged) stop("variance components did not converge!")

