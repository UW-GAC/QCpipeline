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


for (i in 1:length(actions))
{
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.RData", sep="")
  print(fname)
  combined <- getobj(fname)
  
  png(paste(qqfname,".model.", i, ".",actions[i],".png",sep=""), width=1200, height=600)
  par(mfrow=c(1,2))

  test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "\n", model.type[i])
  print(test)
  varp <- paste("model.",i,".", actions[i], ".pvalue.G", sep="")
  varl <- paste("model.",i,".", actions[i], ".chisq.G", sep="")
  pval <- combined[,varp]

  # unfiltered plot
  pvaln <- pval[!is.na(pval)]
  lambda <- median(-2*log(pval[!is.na(pval)]))/ 1.39 # change to new way to calculating lambda
  print(lambda)
  qqPlot(pvaln, trunc=F, main=paste(test, ",unfiltered", sep=""), cex.main = 0.9, sub=paste("lambda =",format(lambda,digits=4)))

  # add filters
  stopifnot(all(combined$snpID %in% snpAnnot$snpID))
  combined$quality.filter <- snpAnnot[match(combined$snpID,snpAnnot$snpID), qf]
  combined$qual.maf.filter <- combined$quality.filter & (!is.na(combined$minor.allele)) & combined$MAF>0.05 & combined$MAF<0.95
  fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
  save(combined, file=fname)

  pvaln <- pval[combined$quality.filter]
  lambda <- median(-2*log(pvaln[!is.na(pvaln)]))/ 1.39 # change to new way to calculating lambda
  print(lambda)
  qqPlot(pvaln, trunc=F, main=paste(test, ",\nfiltered", sep=""), cex.main = 0.9, sub=paste("lambda =",format(lambda,digits=4)))

  pvaln <- pval[combined$qual.maf.filter]
  lambda <- median(-2*log(pvaln[!is.na(pvaln)]))/ 1.39 # change to new way to calculating lambda
  print(lambda)
  qqPlot(pvaln, trunc=F, main=paste(test, ",\nMAF filtered", sep=""), cex.main = 0.9, sub=paste("lambda =",format(lambda,digits=4)))
  
  dev.off()

  # Manhattan plots - filtered
  png(paste(qqfname,".model.", i, ".",actions[i],".manh.filt.png",sep=""), width=1200, height=600)
  chromosome <- snpAnnot$chromosome[match(combined$snpID,snpAnnot$snpID)][combined$quality.filter]
  idx <- 1:26 %in% unique(chromosome)
  manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=c(1:22,"X","Y","XY","M")[idx],
                main=paste(test,"- filtered"), cex.main=1.5)
  dev.off() 

    # Manhattan plots - maf filtered
  png(paste(qqfname,".model.", i, ".",actions[i],".manh.maf.filt.png",sep=""), width=1200, height=600)
  chromosome <- snpAnnot$chromosome[match(combined$snpID,snpAnnot$snpID)][combined$qual.maf.filter]
  idx <- 1:26 %in% unique(chromosome)
  manhattanPlot(p=pvaln,chromosome=chromosome,chrom.labels=c(1:22,"X","Y","XY","M")[idx],
                main=paste(test,"- MAF filtered"), cex.main=1.5)
  dev.off()
}


