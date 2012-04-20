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
required <- c("annot_snp_file", "assoc_output", "covar.list", "gene_action", "model_type", "outcome")
optional <- c("maf.filter", "plot_chroms", "plot_out", "quality.filter", "signif_line")
default <- c(0.02, NA, "assoc", "quality.filter", 5e-8)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# variables
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
pathprefix <- config["assoc_output"]
pathprefix
actions <-  config["gene_action"]
actions <- unlist(strsplit(actions," "))
actions
qqfname <- config["plot_out"]
qqfname
outcome <- config["outcome"]
outcome <- unlist(strsplit(outcome," "))
outcome
covar.list <- getobj(config["covar.list"])
covar.list
model.type <- config["model_type"]
model.type <- unlist(strsplit(model.type," "))
stopifnot(all(model.type %in% c("logistic", "linear", "Logistic", "Linear")))
model.type
qf <- config["quality.filter"]
qf
if (!is.na(config["plot_chroms"])) {
  plotchroms <- getobj(config["plot_chroms"])
  plotchroms
}
sub <- NULL

for (i in 1:length(actions))
{
  test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "-", model.type[i])
  print(test)
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.RData", sep="")
  print(fname)
  combined <- getobj(fname)
  # only keep chroms in plotchroms
  if (!is.na(config["plot_chroms"]))
  {
    sub <- combined$snpID %in% snpID[getChromosome(snpAnnot) %in% plotchroms]
  } else {
    sub <- rep(TRUE, nrow(combined))
  }
  
  # add filters, NOT subsetted with plotchroms
  stopifnot(all(combined$snpID %in% snpID))
  combined$quality.filter <- getVariable(snpAnnot, qf, index=match(combined$snpID,snpID))
  combined$qual.maf.filter <- combined$quality.filter & (!is.na(combined$minor.allele)) & combined$MAF>as.numeric(config["maf.filter"])
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  save(combined, file=fname)

  for (type in c("LR", "Wald")) { 
    # QQ plots
    png(paste(qqfname,"_model_", i, "_",actions[i],"_qq_",type,".png",sep=""), width=720, height=720)
    par(mfrow=c(2,2), mar=c(5,5,4,2)+0.1, lwd=1.5,
        cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
    varp <- paste("model", i, actions[i], type,"pval.G", sep=".")
    pval <- combined[,varp]

    # QQ plots - unfiltered plot, subsetted with plotchroms
    pvaln <- pval[!is.na(pval) & sub]
    title <- paste("no filter\nN =", length(pvaln))
    qqPlot(pvaln, trunc=F, main=title)

    # QQ plots - filtered, subsetted with plotchroms
    pvaln <- pval[combined$quality.filter & (!is.na(pval)) & sub]
    title <- paste("quality filter\nN =", length(pvaln))
    qqPlot(pvaln, trunc=F, main=title)

    # QQ plots - maf filtered, subsetted with plotchroms
    pvaln <- pval[combined$qual.maf.filter & (!is.na(pval)) & sub]
    title <- paste("quality filter + MAF >",config["maf.filter"],"\nN =", length(pvaln))
    qqPlot(pvaln, trunc=F, main=title)

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
    pvaln <- pval[(!is.na(pval)) & sub]
    chromosome <- getChromosome(snpAnnot, index=match(combined$snpID[(!is.na(pval)) & sub],snpID), char=TRUE)
    title <- paste("no filter\nN =", length(pvaln))
    manhattanPlot(p=pvaln,chromosome=chromosome,
                  main=title, signif=as.numeric(config["signif_line"]))

    # Manhattan plots - filtered, subsetted with plotchroms
    pvaln <- pval[combined$quality.filter & (!is.na(pval)) & sub]
    chromosome <- getChromosome(snpAnnot, index=match(combined$snpID[combined$quality.filter & (!is.na(pval)) & sub],snpID), char=TRUE)
    title <- paste("quality filter\nN =", length(pvaln))
    manhattanPlot(p=pvaln,chromosome=chromosome,
                  main=title, signif=as.numeric(config["signif_line"]))

    # Manhattan plots - maf filtered, subsetted with plotchroms
    chromosome <- getChromosome(snpAnnot, index=match(combined$snpID[combined$qual.maf.filter & (!is.na(pval)) & sub],snpID), char=TRUE)
    pvaln <- pval[combined$qual.maf.filter & (!is.na(pval)) & sub]
    title <- paste("quality filter + MAF >",config["maf.filter"],"\nN =", length(pvaln))
    manhattanPlot(p=pvaln,chromosome=chromosome,
                  main=title, signif=as.numeric(config["signif_line"]))
    dev.off()
  }

}


