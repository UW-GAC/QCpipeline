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
pathprefix <- config["out_assoc_prefix"]
pathprefix
actions <- unlist(strsplit(config["gene_action"]," "))
actions
qqfname <- config["out_plot_prefix"]
qqfname
outcome <- unlist(strsplit(config["outcome"]," "))
outcome
covar.list <- getobj(config["covar.list"])
covar.list
model.type <- unlist(strsplit(config["model_type"]," "))
stopifnot(all(model.type %in% c("logistic", "linear")))
model.type
if (!is.na(config["plot_chroms"])) {
  plotchroms <- getobj(config["plot_chroms"])
  plotchroms
}

for (i in 1:length(actions))
{
  test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "-", model.type[i])
  print(test)
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  print(fname)
  combined <- getobj(fname)
  # only keep chroms in plotchroms
  if (!is.na(config["plot_chroms"]))
  {
    sub <- combined$snpID %in% snpID[getChromosome(snpAnnot) %in% plotchroms]
    combined <- combined[sub,]
  }

  maf.thresh <- switch(model.type[i],
                       linear=config["maf.linear.threshold"],
                       logistic=config["maf.logistic.threshold"])
  maf.text <- switch(config["maf.filter.type"],
                     absolute=paste("MAF >", config["maf.absolute.threshold"]),
                     snp.specific=paste("2*MAF*(1-MAF)*N >", maf.thresh))
  
  for (type in c("LR", "Wald")) { 
    # QQ plots
    png(paste(qqfname,"_model_", i, "_",actions[i],"_qq_",type,".png",sep=""), width=720, height=720)
    par(mfrow=c(2,2), mar=c(5,5,4,2)+0.1, lwd=1.5,
        cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
    varp <- paste(type,"pval", sep=".")
    pval <- combined[,varp]

    # QQ plots - unfiltered plot, subsetted with plotchroms
    pvaln <- pval[!is.na(pval)]
    title <- paste("no filter\n", length(pvaln), "SNPs")
    lambda <- median(-2*log(pvaln), na.rm=TRUE) / 1.39
    subtitle <- paste("lambda =", format(lambda, digits=4, nsmall=3))
    qqPlot(pvaln, trunc=FALSE, main=title, sub=subtitle)

    # QQ plots - filtered, subsetted with plotchroms
    pvaln <- pval[combined$quality.filter & (!is.na(pval))]
    title <- paste("quality filter\n", length(pvaln), "SNPs")
    lambda <- median(-2*log(pvaln), na.rm=TRUE) / 1.39
    subtitle <- paste("lambda =", format(lambda, digits=4, nsmall=3))
    qqPlot(pvaln, trunc=FALSE, main=title, sub=subtitle)

    # QQ plots - maf filtered, subsetted with plotchroms
    pvaln <- pval[combined$qual.maf.filter & (!is.na(pval))]
    title <- paste("quality filter +",maf.text,"\n", length(pvaln), "SNPs")
    lambda <- median(-2*log(pvaln), na.rm=TRUE) / 1.39
    subtitle <- paste("lambda =", format(lambda, digits=4, nsmall=3))
    qqPlot(pvaln, trunc=FALSE, main=title, sub=subtitle)

    # obs-exp plot
    pvalx <- -log10(sort(pvaln)) # sort() removes NAs
    n <- length(pvalx)
    x <- -log10((1:n)/n)
    plot(x, pvalx-x, xlab=substitute(paste(-log[10], "(expected P)")),
         ylab=substitute(paste(log[10], "(expected P)", -log[10], "(observed P)")),
         main=title)
    abline(h=0,col="red")
    dev.off()
    
    # Manhattan plots - no filter, subsetted with plotchroms
    png(paste(qqfname,"_model_", i, "_",actions[i],"_manh_",type,".png",sep=""), width=720, height=720)
    par(mfrow=c(3,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
    pvaln <- pval[(!is.na(pval))]
    chromosome <- getChromosome(snpAnnot, index=match(combined$snpID[(!is.na(pval))],snpID), char=TRUE)
    title <- paste("no filter\n", length(pvaln), "SNPs")
    manhattanPlot(p=pvaln,chromosome=chromosome,
                  main=title, signif=as.numeric(config["signif_line"]))

    # Manhattan plots - filtered, subsetted with plotchroms
    pvaln <- pval[combined$quality.filter & (!is.na(pval))]
    chromosome <- combined$chromosome[combined$quality.filter & (!is.na(pval))]
    title <- paste("quality filter\n", length(pvaln), "SNPs")
    manhattanPlot(p=pvaln,chromosome=chromosome,
                  main=title, signif=as.numeric(config["signif_line"]))

    # Manhattan plots - maf filtered, subsetted with plotchroms
    pvaln <- pval[combined$qual.maf.filter & (!is.na(pval))]
    chromosome <- combined$chromosome[combined$qual.maf.filter & (!is.na(pval))]
    title <- paste("quality filter +",maf.text,"\n", length(pvaln), "SNPs")
    manhattanPlot(p=pvaln,chromosome=chromosome,
                  main=title, signif=as.numeric(config["signif_line"]))
    dev.off()
  }

}


