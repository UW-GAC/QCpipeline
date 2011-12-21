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

for (i in 1:length(actions))
{
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
  
  png(paste(qqfname,".model.", i, ".",actions[i],".png",sep=""), width=720, height=720)
  par(mfrow=c(2,2), mar=c(5,5,4,2)+0.1, lwd=1.5,
      cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
  test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "\n", model.type[i])
  print(test)
  varp <- paste("model.",i,".", actions[i], ".pvalue.G", sep="")
  varl <- paste("model.",i,".", actions[i], ".chisq.G", sep="")
  pval <- combined[,varp]

  # unfiltered plot, subsetted with plotchroms
  pvaln <- pval[!is.na(pval) & sub]
  qqPlot(pvaln, trunc=F, main=paste(test, ", unfiltered", sep=""))

  # add filters, NOT subsetted with plotchroms
  stopifnot(all(combined$snpID %in% snpAnnot$snpID))
  combined$quality.filter <- snpAnnot[match(combined$snpID,snpAnnot$snpID), qf]
  combined$qual.maf.filter <- combined$quality.filter & (!is.na(combined$minor.allele)) & combined$MAF>0.05 & combined$MAF<0.95
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  save(combined, file=fname)

  # QQ plots - filtered, subsetted with plotchroms
  pvaln <- pval[combined$quality.filter & (!is.na(pval)) & sub]
  qqPlot(pvaln, trunc=F, main=paste(test, ", filtered", sep=""))

  # QQ plots - maf filtered, subsetted with plotchroms
  pvaln <- pval[combined$qual.maf.filter & (!is.na(pval)) & sub]
  qqPlot(pvaln, trunc=F, main=paste(test, ", MAF filtered", sep=""))
  dev.off()

  # Manhattan plots - no filter, subsetted with plotchroms
  png(paste(qqfname,".model.", i, ".",actions[i],".manh.png",sep=""), width=720, height=720)
  par(mfrow=c(3,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
  #png(paste(qqfname,".model.", i, ".",actions[i],".manh.no.filt.png",sep=""), width=1200, height=600)
  # chromosome <- snpAnnot$chromosome[match(combined$snpID,snpAnnot$snpID)][combined$quality.filter]
  pvaln <- pval[(!is.na(pval)) & sub]
  chromosome <- snpAnnot$chromosome[match(combined$snpID[(!is.na(pval)) & sub],snpAnnot$snpID)]
  chroms <- 23:26
  names(chroms) <- c("X","XY","Y","M")
  idx <- match(chroms, unique(chromosome))
  chrom.labels <- unique(chromosome)
  chrom.labels[idx[!is.na(idx)]] <- names(chroms)[!is.na(idx)]
  manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=chrom.labels,
                main=paste(test,"- not filtered"))
  #dev.off() 

  # Manhattan plots - filtered, subsetted with plotchroms
  #png(paste(qqfname,".model.", i, ".",actions[i],".manh.filt.png",sep=""), width=1200, height=600)
  # chromosome <- snpAnnot$chromosome[match(combined$snpID,snpAnnot$snpID)][combined$quality.filter]
  pvaln <- pval[combined$quality.filter & (!is.na(pval)) & sub]
  chromosome <- snpAnnot$chromosome[match(combined$snpID[combined$quality.filter & (!is.na(pval)) & sub],snpAnnot$snpID)]
  idx <- match(chroms, unique(chromosome))
  chrom.labels <- unique(chromosome)
  chrom.labels[idx[!is.na(idx)]] <- names(chroms)[!is.na(idx)]
  manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=chrom.labels,
                main=paste(test,"- filtered"))
  #dev.off() 


  # Manhattan plots - maf filtered, subsetted with plotchroms
  #png(paste(qqfname,".model.", i, ".",actions[i],".manh.maf.filt.png",sep=""), width=1200, height=600)
  chromosome <- snpAnnot$chromosome[match(combined$snpID[combined$qual.maf.filter & (!is.na(pval)) & sub],snpAnnot$snpID)]
  pvaln <- pval[combined$qual.maf.filter & (!is.na(pval)) & sub]
  manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=chrom.labels,
                main=paste(test,"- MAF filtered"))
  dev.off()
}


