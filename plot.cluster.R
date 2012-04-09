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
if (!is.na(config["plot_chroms"])) {
  plotchroms <- getobj(config["plot_chroms"])
  plotchroms
}
sub <- NULL

# make genotypedata and intensityData
scanAnnot <- getobj(config["annot_scan_file"])
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
chrom <- getChromosome(snpAnnot)
rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"])
nc <- NcdfGenotypeReader(config["nc_geno_file_samp"])
(genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)) 
xync <- NcdfIntensityReader(config["nc_xy_file_samp"])
(xyData <- IntensityData(xync, scanAnnot=scanAnnot, snpAnnot=snpAnnot))


for (i in 1:length(actions))
{
      test <- paste(outcome[i],"~", paste(covar.list[[i]], collapse=" + "), "-", model.type[i])
      print(test)
      fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
      print(fname)
      combined <- getobj(fname)
      if (!is.na(config["plot_chroms"]))
      {
         sub <- combined$snpID %in% snpID[chrom %in% plotchroms]
      } else {
         sub <- rep(TRUE, nrow(combined))
      }
      combined <- combined[combined$quality.filter & sub,]
      
      varp <- paste("model.",i,".", actions[i], ".LR.pval.G", sep="")
      combined.intid <- combined[order(combined[,varp]),c("snpID",varp)]
      snp.intid <- combined.intid[1:27,]
      
      pdf(paste(qqfname,"_model_", i, "_",actions[i],"_lowP_hits.pdf",sep="")) 
      par(mfrow=c(3,3))
      ind <- match(snp.intid$snpID, snpID)
      text <- paste(rsID[ind], "Chr", chrom[ind])
      mtxt <- paste(text,"\np-value",sprintf("%.2e",
                    snp.intid[,paste("model.",i, ".", actions[i], ".LR.pval.G", sep="")]))

      # plot
      genoClusterPlot(xyData,genoData, plot.type="RTheta", snp.intid$snpID, mtxt)
      dev.off()

      # single page png for QC report
      png(paste(qqfname,"_model_", i, "_",actions[i],"_lowP_hits.png",sep=""), width=720, height=720) 
      par(mfrow=c(3,3), mar=c(5,5,4,2)+0.1, lwd=1.5,
          cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
      genoClusterPlot(xyData,genoData, plot.type="RTheta", snp.intid$snpID[1:9], mtxt[1:9])
      dev.off()
      
}

