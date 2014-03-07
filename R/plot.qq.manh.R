##########
# QQ and Manhattan plots for association tests
# Usage: R --args config.file < plot.qq.manh.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "out_assoc_prefix", "covar.list", "gene_action",
              "model_type", "outcome")
optional <- c("plot_chroms", "out_plot_prefix", "signif_line", "maf.filter.type",
              "maf.absolute.threshold", "maf.linear.threshold", "maf.logistic.threshold")
default <- c(NA, "assoc", 5e-8, "snp.specific", 0.02, 30, 50)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# variables
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
(pathprefix <- config["out_assoc_prefix"])
(actions <- unlist(strsplit(config["gene_action"]," ")))
(qqfname <- config["out_plot_prefix"])
(outcome <- unlist(strsplit(config["outcome"]," ")))
(covar.list <- getobj(config["covar.list"]))
(model.type <- unlist(strsplit(config["model_type"]," ")))
stopifnot(all(model.type %in% c("logistic", "linear")))
if (!is.na(config["plot_chroms"])) {
  (plotchroms <- getobj(config["plot_chroms"]))
}

for (i in 1:length(actions)) {
  test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "-", model.type[i])
  print(test)
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  print(fname)
  combined <- getobj(fname)
  # only keep chroms in plotchroms
  if (!is.na(config["plot_chroms"])) {
    sub <- combined$snpID %in% snpID[getChromosome(snpAnnot) %in% plotchroms]
    combined <- combined[sub,]
  }

  maf.thresh <- switch(model.type[i],
                       linear=config["maf.linear.threshold"],
                       logistic=config["maf.logistic.threshold"])
  #maf.text <- switch(config["maf.filter.type"],
  #                   absolute=paste("MAF >", config["maf.absolute.threshold"]),
  #                   snp.specific=paste("2*MAF*(1-MAF)*N >", maf.thresh))
  
  # calculate MAF
  if (config["maf.filter.type"] == "absolute") {
    mafhi.text <- paste("MAF >", config["maf.absolute.threshold"])
    mafhi.string <- ""
    maflo.text <- paste("MAF <=", config["maf.absolute.threshold"])
    maflo.string <- ""
  } else if (config["maf.filter.type"] == "snp.specific") {
    N <- switch(model.type[i],
                linear=max(combined$n, na.rm=TRUE),
                logistic=max( pmin( combined$nAA.cc0 + combined$nAB.cc0 + combined$nBB.cc0,
                  combined$nAA.cc1 + combined$nAB.cc1 + combined$nBB.cc1), na.rm=TRUE))
    
    maf <- quadSolveMAF(as.numeric(maf.thresh), N)
    mafhi.text <- paste("2*MAF*(1-MAF)*N >", maf.thresh)
    mafhi.string <- paste("- MAF >", format(maf, digits=3))
    maflo.text <- paste("2*MAF*(1-MAF)*N <=", maf.thresh)
    maflo.string <- paste("- MAF <=", format(maf, digits=3))
  }

  for (type in c("LR", "Wald")) { 
    varp <- paste(type,"pval", sep=".")    
    ## no LR test for models with interactions
    if (!(varp %in% names(combined))) next
  
    ## filtering
    assoc <- combined[!is.na(combined[,varp]),
                      c(varp, "chromosome", "composite.filter", "comp.maf.filter")]
    names(assoc)[1] <- "pval"
    title.no.filt <- paste("no filter\n", nrow(assoc), "SNPs")
    qual.filt <- which(assoc$composite.filter)
    title.qual.filt <- paste("composite filter\n", length(qual.filt), "SNPs")
    mafhi.filt <- which(assoc$comp.maf.filter)
    title.mafhi.filt <- paste("composite.filter +",mafhi.text,"\n", length(mafhi.filt), "SNPs", mafhi.string)
    maflo.filt <- which(assoc$composite.filter & !assoc$comp.maf.filter)
    title.maflo.filt <- paste("composite.filter +",maflo.text,"\n", length(maflo.filt), "SNPs", maflo.string)
    filters <- list(1:nrow(assoc), qual.filt, mafhi.filt, maflo.filt)
    names(filters) <- c(title.no.filt, title.qual.filt, title.mafhi.filt, title.maflo.filt)
    
  
    # QQ plots
    outfile <- paste(qqfname,"_model_", i, "_",actions[i],"_qq_",type,".png",sep="")
    qqPlotPng(assoc$pval, filters, outfile)
  
    # Manhattan plots
    outfile <- paste(qqfname,"_model_", i, "_",actions[i],"_manh_",type,".png",sep="")
    manhattanPlotPng(assoc$pval, assoc$chromosome, filters[1:3], outfile)
  }
}


