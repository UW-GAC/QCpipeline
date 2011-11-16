##########
# HWE plots
# Usage: R --args config.file < hwe_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

hwe <- getobj(paste(config["out_hwe_prefix"], "RData", sep="."))

# qq plots
png(config["out_qq_plot"], width=720, height=720)
par(mfrow=c(2,2))
qqPlot(hwe$p.value[hwe$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(hwe$p.value[hwe$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")
qqPlot(hwe$p.value[hwe$chromosome == 23], trunc=FALSE, main="X chromosome, all")
qqPlot(hwe$p.value[hwe$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
dev.off()


# inbreeding coefficient
aut <- hwe$chromosome < 23
summary(hwe$f)
summary(hwe$f[aut])

png(config["out_inbrd_plot"], width=600, height=600)
hist(hwe$f[aut], xlab="Inbreeding coefficient estimate", ylab="number of SNPs",
     main="Autosomes, all")
abline(v=mean(hwe$f, na.rm=TRUE), lty=2, col="gray")
dev.off()


# MAF vs Pvalue
png(config["out_maf_plot"], width=720, height=720)
plot(hwe$MAF[aut], -log10(hwe$p.value[aut]), xlab="MAF", ylab="-log10(p-value)",
     main="MAF vs P-value - autosomal")
dev.off()


# cluster plots
# To decide what p-value to use for the hwe filter, we look at cluster plots
# for a sample of SNPs, binned by their p-value
# get missing call rate
snpAnnot <- getobj(config["annot_snp_file"])
snp <- getVariable(snpAnnot, c("snpID", config["annot_snp_rsIDCol"], config["annot_snp_missingCol"]))
names(snp) <- c("snpID", "rsID", "missing")
hwe <- merge(hwe, snp)
# alog: autosomal snps with non-NA p values and an MCR < .05
alog <- hwe$chromosome<23 & !is.na(hwe$p.value) & hwe$missing<0.05
table(alog)

x <- c(0, 1e-24, 1e-12, 1e-6, 1e-4, 1e-2, 1e-1, 0.25, 0.50, 0.75, 1)
bins <- rep(NA, (length(x)-1))
ids <- list()
for(i in 1:(length(x)-1)){
  ids[[i]] <- hwe$snpID[ alog & hwe$p.value>=x[i] &
                             hwe$p.value<x[i+1] ]
  
  # put indices of snps with pvalue = 1 in the last bin
  if(i==(length(x)-1)) {
    ids[[i]] <- c(ids[[i]], hwe$snpID[alog & hwe$p.value==1])
  }

  # number of snps in each bin
  bins[i] <- length(ids[[i]]) 
}
sum(bins) #  - covers all snps in alog
sum(unlist(lapply(ids,length)))
bins

xyNC <- NcdfIntensityReader(config["nc_qxy_file"])
xyData <- IntensityData(xyNC, snpAnnot=snpAnnot)
genoNC <- NcdfGenotypeReader(config["nc_geno_file"])
genoData <- GenotypeData(genoNC, snpAnnot=snpAnnot)

for(j in 1:3){  # independent sets of samples, each in a different png file
  # sample n from each bin
  n <- 9
  sids <- lapply(ids, function(x) {if (length(x) > n) sample(x, n) else x}) 
  sids <- unlist(sids)
  # retrieve all sampled data
  dat <- hwe[hwe$snpID %in% sids,
                         c("snpID", "chromosome", "rsID", "p.value")]
  dat <- dat[order(dat$p.value),] # ordered by pvalues
  
  mtxt <- paste("Chr",dat$chromosome,dat$rsID, "\np-value", format(dat$p.value,digits=3))
  png(paste(config["out_clust_prefix"], "_",
            j,"%03d.png",sep=""),width=720,height=720)
  par(mfrow=c(3,3))
  genoClusterPlot(xyData, genoData, snpID=dat$snpID, main.txt=mtxt)
  dev.off()
}

## get snp counts based upon various thresholds
sum(hwe$p.value < 1e-6, na.rm=TRUE)
sum(hwe$p.value < 1e-5, na.rm=TRUE)
sum(hwe$p.value < 1e-4, na.rm=TRUE)
sum(hwe$p.value < 1e-3, na.rm=TRUE)

#compare 1e-2 to 1e-4 versus 1e-4 to 1e-6
#these are bins 4-5
pbin <- c("1e-6_1e-4", "1e-4_1e-2")
for (i in 4:5) {
  # sample n from each bin
  n <- 9*7
  sids <- unlist(if (length(ids[[i]]) > n) sample(ids[[i]], n) else ids[[i]])
  if (length(sids) > 0) {
  # retrieve all sampled data
    dat <- hwe[hwe$snpID %in% sids,
                         c("snpID", "chromosome", "rsID", "p.value")]
    dat <- dat[order(dat$p.value),] # ordered by pvalues

    mtxt <- paste("Chr",dat$chromosome,dat$rsID, "\np-value", format(dat$p.value,digits=3))
    png(paste(config["out_clust_prefix"], "_",
            pbin[i-3],"%03d.png",sep=""),width=720,height=720)
    par(mfrow=c(3,3))
    genoClusterPlot(xyData, genoData, snpID=dat$snpID, main.txt=mtxt)
    dev.off()
  }
}

