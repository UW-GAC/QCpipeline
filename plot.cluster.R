library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)


# variables
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
afs <- config["A.freq.study"]
afs

# make genotypedata and intensityData
scanAnnot <- getobj(config["annot_scan_file"])
snpAnnot <- getobj(config["annot_snp_file"])
nc <- NcdfGenotypeReader(config["nc_geno_file_samp"])
sid <- getScanID(nc)
stopifnot(all(scanAnnot$scanID==sid))
(genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)) 
xync <- NcdfIntensityReader(config["nc_xy_file_samp"])
(xyData <- IntensityData(xync,snpAnnot=snpAnnot))


for (i in 1:length(actions))
{
      fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
      print(fname)
      combined <- getobj(fname)
      combined <- combined[combined$quality.filter,]
      
      varp <- paste("model.",i,".", actions[i], ".pvalue.G", sep="")
      combined.intid <- combined[order(combined[,varp]),c("snpID",varp)]
      snp.intid <- combined.intid[1:27,]
      
      pdf(paste(qqfname,".model.", i, ".",actions[i],".lowP.hits.pdf",sep="")) 
      par(mfrow=c(3,3))        
      text <- paste(pData(snpAnnot)[match(snp.intid$snpID, snpAnnot$snpID), "rsID"],
                    "Chr",snpAnnot$chromosome[match(snp.intid$snpID, snpAnnot$snpID)])
      mtxt <- paste(text,"\np-value",sprintf("%.2e",
                    snp.intid[,paste("model.",i, ".", actions[i], ".pvalue.G", sep="")]))

      # plot
      test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "\n", model.type[i])
      genoClusterPlot(xyData,genoData, plot.type="RTheta", snp.intid$snpID, mtxt,cex.main=0.85)
      dev.off() 
}

