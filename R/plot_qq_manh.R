##########
# QQ and Manhattan plots for association tests
# Usage: R --args config.file < plot.qq.manh.R
##########


## Some notes:
# to calculate df for interaction terms (if you are going to pvals for them, future reader of this code)
# Main effect: df = 1 (always)
# interaction terms (say you want interactions between genotype and race and between genotype and gengrp):
#   df for race is k-1, so if k=3 then df=2
#   df for gengrp is also k-1, so if k=6 then df=5
# joint: everything having to do with genotype
#   for a model with main effect and race, gengrp interactions:
#   df = df_mainEffect + df_raceInteraction + df_gengrpInteraction
#   e.g., df = 1 + (3-1) + (6-1) = 8
# NOTE if you are doing an interaction with a quantitative term, then df=1 for that interaction

# assocRegression returns the Wald stat as Wald.Stat = (Beta / SE)^2. It also returns a likelihood ratio statistic, LR.Stat, which is 2*[ log(likelihood of model with genotype) - log(likelihood of null) ]
# 
# That means we can use the same formula, e.g., median(Wald.Stat) / qchisq(0.5, df=1) and median(LR.Stat) / qchisq(0.5, df=1), for both types of test.

# a possibly helpful reference: PMC3093379 (Voorman et al 2011)

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("model_type", "annot_scan_file", "geno_file", "covars_as_factor")
optional <- c("ivar", "plot_chroms", "out_plot_prefix", "signif_line", "maf.filter.type",
              "maf.absolute.threshold", "maf.linear.threshold", "maf.logistic.threshold",
              "out_assoc_prefix", "scan_exclude")
default <- c(NA, NA, "assoc", 5e-8, "snp.specific", 0.02, 30, 50, "assoc", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# variables
(pathprefix <- config["out_assoc_prefix"])
(qqfname <- config["out_plot_prefix"])
(model.type <- config["model_type"])
stopifnot(all(model.type %in% c("logistic", "linear")))
if (!is.na(config["plot_chroms"])) {
    ## process plot_chroms - works with e.g., "1 2 3 5:10"
    plotchroms <- sort(unique(unlist(lapply(unlist(strsplit(config["plot_chroms"], " ")),
                                            function(x) eval(parse(text=x))))))
}

fname <- paste0(pathprefix, "_combined_qual_filt.RData")
combined <- getobj(fname)
## only keep chroms in plotchroms
if (!is.na(config["plot_chroms"])) {
    combined <- combined[combined$chr %in% plotchroms,]
}

maf.thresh <- switch(model.type,
                     linear=config["maf.linear.threshold"],
                     logistic=config["maf.logistic.threshold"])
##maf.text <- switch(config["maf.filter.type"],
##                   absolute=paste("MAF >", config["maf.absolute.threshold"]),
##                   snp.specific=paste("2*MAF*(1-MAF)*N >", maf.thresh))

## calculate MAF
if (config["maf.filter.type"] == "absolute") {
    mafhi.text <- paste("MAF >", config["maf.absolute.threshold"])
    mafhi.string <- ""
    maflo.text <- paste("MAF <=", config["maf.absolute.threshold"])
    maflo.string <- ""
} else if (config["maf.filter.type"] == "snp.specific") {
    N <- switch(model.type,
                linear=max(combined$n, na.rm=TRUE),
                logistic=max( pmin(combined$n0, combined$n1), na.rm=TRUE))
    
    maf <- quadSolveMAF(as.numeric(maf.thresh), N)
    mafhi.text <- paste("2*MAF*(1-MAF)*N >", maf.thresh)
    mafhi.string <- paste("- MAF >", format(maf, digits=3))
    maflo.text <- paste("2*MAF*(1-MAF)*N <=", maf.thresh)
    maflo.string <- paste("- MAF <=", format(maf, digits=3))
}

## for interaction, find degrees of freedom
if (!is.na(config["ivar"])) {
    ## get scan annotation with only samples used in testing
    scanAnnot <- getobj(config["annot_scan_file"])
    data <- GenotypeReader(config["geno_file"])
    sid <- getScanID(data)
    scanAnnot <- scanAnnot[match(sid, scanAnnot$scanID),]
    close(data)
    if (!is.na(config["scan_exclude"])) {
        scan.exclude <- getobj(config["scan_exclude"])
        scanAnnot <- scanAnnot[!(scanAnnot$scanID %in% scan.exclude),]
    }
    
    ## if this is quantitative, df=1, otherwise, it's the number of levels in the factor - 1
    ivar <- config["ivar"]
    factors <- unlist(strsplit(config["covars_as_factor"], " ", fixed=TRUE))
    ge.df <- ifelse(ivar %in% factors, nlevels(as.factor(scanAnnot[[ivar]])) - 1, 1)
}

for (type in c("LR", "Wald", "GxE", "Joint")) { 
    varp <- paste(type, "pval", sep=".")    
    varstat <- paste(type, "Stat", sep=".")
    
    if (!(varp %in% names(combined))) next
    
    ## filtering
    assoc <- combined[!is.na(combined[,varp]),
                      c(varp, varstat, "chromosome", "composite.filter", "comp.maf.filter")]
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
    
    
    ## The wald.stat and lr.stat have 1 degree of freeom
    ## if we were plotting joint p-values from a model with interactions, it would be higher (2 or more)
    df <- switch(type, Wald=1, LR=1, GxE=ge.df, Joint=ge.df+1)
    
    ## QQ plots
    outfile <- paste0(qqfname, "_qq_", type, ".png")
    qqPlotPng(assoc$pval, stat=assoc[[varstat]], df=df, filters, outfile)
    
    ## Manhattan plots
    outfile <- paste0(qqfname, "_manh_", type, ".png")
    manhattanPlotPng(assoc$pval, assoc$chromosome, filters[1:3], outfile)
}



