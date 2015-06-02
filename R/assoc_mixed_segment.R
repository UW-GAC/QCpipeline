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

if (length(args) < 2) stop("missing chromosome")
(chromosome <- as.numeric(args[2]))

if (length(args) < 2) stop("missing segment")
(segment <- as.numeric(args[3]))

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
scanAnnot <- getobj(paste0(config["out_assoc_prefix"], "_scanAnnot.RData"))

# set categorical variables in association models as factor
if (!is.na(config["covars_as_factor"])) {
  factors <- unlist(strsplit(config["covars_as_factor"], " ", fixed=TRUE))
  for (i in factors) {
    scanAnnot[[i]] <- as.factor(scanAnnot[[i]])
  }
}

# make geno data
genoByChr <- GenotypeDataByChr(config["geno_file"])
genoData <- getGenoData(genoByChr, chromosome=chromosome, snpAnnot=TRUE, scanAnnot=scanAnnot)


# variance components
VC <- getobj(paste0(config["out_assoc_prefix"], "_VC.RData"))


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
  test_type <- "Score"
} else if (model.type == "linear") {
  test_type <- "Wald"
} else {
  stop("'model_type' must be 'logsitic' or 'linear'")
}
test_type


# scan.include
if (!is.na(config["scan_include_file"])) {
  scan.include <- getobj(config["scan_include_file"])
  scan.exclude <- setdiff(scanAnnot$scanID, scan.include)
} else {
  scan.include <- NULL
  scan.exclude <- NULL
}

# snp selection
#snpseg <- read.csv(config["snp_segment_file"], as.is=T)
snpAnnot <- getSnpAnnotation(genoData)
ids <- snpAnnot$snpID[snpAnnot$chromosome %in% chromosome & snpAnnot$segment %in% segment]


(idr <- range(ids))
index <- which(getSnpID(genoData) %in% idr)

block.size <- as.integer(config["block_size"])
impute.geno <- as.logical(config["impute_geno"])


# compute variance component estimates
out <- assocTestMixedModel(genoData = genoData,
                           snpStart = index[1],
                           snpEnd = index[2],
                           cholSigmaInv = VC$cholSigmaInv,
                           outcome = "workingY",
                           covar.vec = covars,
                           test=test_type,
                           scan.exclude = scan.exclude,
                           impute.geno = impute.geno,
                           block.size = block.size)


# calculate effective N
if (hasVariable(snpAnnot, "oevar") & hasVariable(snpAnnot, "type")){
  oevar <- getVariable(snpAnnot, "oevar", index[1]:index[2])
  type <- getVariable(snpAnnot, "type", index[1]:index[2])
} else {
  oevar <- NULL
  type <- NULL
}

# if this is not an autosome, we need to recalculate allele frequencies
Xchromosome <- chromosome %in% XchromCode(genoData)
if (Xchromosome){
  afreq <- getAlleleFrequency(genoData, snpStart=index[1], snpEnd=index[2], scan.exclude=scan.exclude, blockSize=block.size, verbose=TRUE)
} else {
  afreq <- out[, c("MAF", "n")]
}


out$effN <- getEffectiveN(afreq, oevar=oevar, type=type, Xchromosome=chromosome %in% XchromCode(genoData))

# add effectiveN for cases and controls if logistic
if (config["model_type"] == "logistic"){
  outcomeVal <- getVariable(scanAnnot, config["outcome"])
  
  # cases
  case.exclude <- union(scan.exclude, scanAnnot$scanID[outcomeVal %in% 0])
  afreq.case <- getAlleleFrequency(genoData, snpStart=index[1], snpEnd=index[2], scan.exclude=case.exclude, blockSize = block.size)
  out$effN.case <- getEffectiveN(afreq.case, oevar=oevar, type=type, Xchromosome=Xchromosome)
  out$n1 <- afreq.case[, "n"]
  
  # controls
  ctrl.exclude <- union(scan.exclude, scanAnnot$scanID[outcomeVal %in% 1])
  afreq.ctrl <- getAlleleFrequency(genoData, snpStart=index[1], snpEnd=index[2], scan.exclude=ctrl.exclude, blockSize = block.size)
  out$effN.ctrl <- getEffectiveN(afreq.ctrl, oevar=oevar, type=type, Xchromosome=Xchromosome)
  out$n0 <- afreq.ctrl[, "n"]
  
  out$effN.all <- out$effN
  out$effN <- pmin(out$effN.case, out$effN.ctrl)
}

close(genoData)

outfile <- paste0(config["out_assoc_prefix"], "_chr-", chromosome, "_seg", segment, ".RData")
save(out, file=outfile)
