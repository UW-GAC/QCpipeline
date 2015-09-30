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
required <- c("annot_snp_file", "samp_geno_file", "samp_xy_file")
optional <- c("annot_snp_rsIDCol", "out_assoc_prefix", "out_plot_prefix")
default <- c("rsID", "assoc", "assoc")
config <- setConfigDefaults(config, required, optional, default)
print(config)


# variables
(pathprefix <- config["out_assoc_prefix"])
(qqfname <- config["out_plot_prefix"])

# make genotypedata and intensityData
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
chrom <- getChromosome(snpAnnot)
rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"])
geno <- GenotypeReader(config["samp_geno_file"])
(genoData <- GenotypeData(geno,  snpAnnot=snpAnnot)) 
xy <- IntensityReader(config["samp_xy_file"])
(xyData <- IntensityData(xy, snpAnnot=snpAnnot))


fname <- paste0(pathprefix, "_combined_qual_filt.RData")
combined <- getobj(fname)
combined <- combined[combined$composite.filter,]

# select pvalue to use for plots
(varp <- intersect(paste0(c("LR", "Wald", "z"), ".pval"), names(combined))[1])

combined.intid <- combined[order(combined[,varp]),c("snpID",varp)]
snp.intid <- combined.intid[1:27,]

pdf(paste0(qqfname, "_lowP_hits.pdf")) 
par(mfrow=c(3,3))
ind <- match(snp.intid$snpID, snpID)
text <- paste(rsID[ind], "Chr", chrom[ind])
mtxt <- paste(text,"\np-value",sprintf("%.2e",snp.intid[,varp]))

## plot
genoClusterPlot(xyData,genoData, plot.type="RTheta", snp.intid$snpID, mtxt)
dev.off()

## single page png for QC report
png(paste0(qqfname, "_lowP_hits.png"), width=720, height=720) 
par(mfrow=c(3,3), mar=c(5,5,4,2)+0.1, lwd=1.5,
    cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
genoClusterPlot(xyData,genoData, plot.type="RTheta", snp.intid$snpID[1:9], mtxt[1:9])
dev.off()     

