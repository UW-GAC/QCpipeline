##########
# Cluster plots for association tests
# Usage: R --args config.file < plot.cluster.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "gene_action", "samp_geno_file", "samp_xy_file")
optional <- c("annot_snp_rsIDCol", "plot_chroms", "out_assoc_prefix", "out_plot_prefix")
default <- c("rsID", NA, "assoc", "assoc")
config <- setConfigDefaults(config, required, optional, default)
print(config)


# variables
(pathprefix <- config["out_assoc_prefix"])
(actions <- unlist(strsplit(config["gene_action"]," ")))
(qqfname <- config["out_plot_prefix"])
if (!is.na(config["plot_chroms"])) {
  (plotchroms <- getobj(config["plot_chroms"]))
}

# make genotypedata and intensityData
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
chrom <- getChromosome(snpAnnot)
rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"])
geno <- GenotypeReader(config["samp_geno_file"])
(genoData <- GenotypeData(geno,  snpAnnot=snpAnnot)) 
xy <- IntensityReader(config["samp_xy_file"])
(xyData <- IntensityData(xy, snpAnnot=snpAnnot))


for (i in 1:length(actions)) {
      fname <- paste(pathprefix, ".model.", i, ".",actions[i], ".combined.qual.filt.RData", sep="")
      print(fname)
      combined <- getobj(fname)
      if (!is.na(config["plot_chroms"])) {
         combined <- combined[combined$chromosome %in% plotchroms,]
      }
      combined <- combined[combined$composite.filter,]
      
      varp <- "LR.pval"
      ## no LR test for models with interactions
      if (!(varp %in% names(combined))) next
      
      combined.intid <- combined[order(combined[,varp]),c("snpID",varp)]
      snp.intid <- combined.intid[1:27,]
      
      pdf(paste(qqfname,"_model_", i, "_",actions[i],"_lowP_hits.pdf",sep="")) 
      par(mfrow=c(3,3))
      ind <- match(snp.intid$snpID, snpID)
      text <- paste(rsID[ind], "Chr", chrom[ind])
      mtxt <- paste(text,"\np-value",sprintf("%.2e",snp.intid[,varp]))

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

