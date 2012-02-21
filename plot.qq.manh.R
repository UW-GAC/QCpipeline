library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# variables
snpAnnot <- pData(getobj(config["annot_snp_file"]))
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
if (is.na(config["signif_line"])) config["signif_line"] <- 5e-8

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
    sub <- combined$snpID %in% snpAnnot$snpID[snpAnnot$chromosome %in% plotchroms]
  } else {
    sub <- rep(TRUE, nrow(combined))
  }
  
  # add filters, NOT subsetted with plotchroms
  stopifnot(all(combined$snpID %in% snpAnnot$snpID))
  combined$quality.filter <- snpAnnot[match(combined$snpID,snpAnnot$snpID), qf]
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
    chromosome <- snpAnnot$chromosome[match(combined$snpID[(!is.na(pval)) & sub],snpAnnot$snpID)]
    chroms <- 23:26
    names(chroms) <- c("X","XY","Y","M")
    idx <- match(chroms, unique(chromosome))
    chrom.labels <- unique(chromosome)
    chrom.labels[idx[!is.na(idx)]] <- names(chroms)[!is.na(idx)]
    title <- paste("no filter\nN =", length(pvaln))
    manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=chrom.labels,
                  main=title, signif=as.numeric(config["signif_line"]))

    # Manhattan plots - filtered, subsetted with plotchroms
    pvaln <- pval[combined$quality.filter & (!is.na(pval)) & sub]
    chromosome <- snpAnnot$chromosome[match(combined$snpID[combined$quality.filter & (!is.na(pval)) & sub],snpAnnot$snpID)]
    idx <- match(chroms, unique(chromosome))
    chrom.labels <- unique(chromosome)
    chrom.labels[idx[!is.na(idx)]] <- names(chroms)[!is.na(idx)]
    title <- paste("quality filter\nN =", length(pvaln))
    manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=chrom.labels,
                  main=title, signif=as.numeric(config["signif_line"]))

    # Manhattan plots - maf filtered, subsetted with plotchroms
    chromosome <- snpAnnot$chromosome[match(combined$snpID[combined$qual.maf.filter & (!is.na(pval)) & sub],snpAnnot$snpID)]
    pvaln <- pval[combined$qual.maf.filter & (!is.na(pval)) & sub]
    title <- paste("quality filter + MAF >",config["maf.filter"],"\nN =", length(pvaln))
    manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=chrom.labels,
                  main=title, signif=as.numeric(config["signif_line"]))
    dev.off()
  }

}


