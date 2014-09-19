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

# check config and set defaults
required <- c("annot_snp_file", "samp_geno_file", "samp_xy_file", "out_hwe_prefix")
optional <- c("annot_snp_missingCol", "annot_snp_rsIDCol", "out_clust_prefix", "out_inbrd_plot", "out_maf_plot", "out_qq_plot", "out_sim_prefix")
default <- c("missing.n1", "rsID", "hwe_clust", "hwe_inbrd.pdf", "hwe_maf.png", "hwe_qq.png", "hwe_sim")
config <- setConfigDefaults(config, required, optional, default)
print(config)

hwe <- getobj(paste(config["out_hwe_prefix"], "RData", sep="."))

# qq plots
# check if any X chrom p-values are valid (in case of all-male study)
plotX <- sum(!is.na(hwe$p.value[hwe$chromosome == 23])) > 0
if (plotX) nrow <- 2 else nrow <- 1
png(config["out_qq_plot"], width=720, height=(360*nrow))
par(mfrow=c(nrow,2), mar=c(5,5,4,2)+0.1, lwd=1.5,
    cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)
qqPlot(hwe$p.value[hwe$chromosome < 23],  trunc=FALSE, main="Autosomes, all")
qqPlot(hwe$p.value[hwe$chromosome < 23], trunc=TRUE, main="Autosomes, truncated")
if (plotX) {
  qqPlot(hwe$p.value[hwe$chromosome == 23], trunc=FALSE, main="X chromosome, all")
  qqPlot(hwe$p.value[hwe$chromosome == 23], trunc=TRUE, main="X chromosome, truncated")
}
dev.off()


# inbreeding coefficient
aut <- hwe$chromosome < 23
summary(hwe$f)
summary(hwe$f[aut])

pdf(config["out_inbrd_plot"], width=6, height=6)
hist(hwe$f[aut], xlab="Inbreeding coefficient estimate",
     main="Autosomal SNPs", breaks=40)
abline(v=mean(hwe$f[aut], na.rm=TRUE), lty=2, col="gray")
abline(v=0, col="red")
dev.off()

# simulated vs observed inbreeding coefficient
sim <- getobj(paste(config["out_sim_prefix"], ".RData", sep=""))
sim <- sim[!is.na(sim$f.sim), ]
d.obs <- density(sim$f.obs)
d.sim <- density(sim$f.sim)
#xlim <- c(min(d.obs$x, d.sim$x), max(d.obs$x, d.sim$x))
xlim <- c(-0.2, 0.2)
ylim <- c(0, max(d.obs$y, d.sim$y))
main.txt <- paste("Simulated F -", nrow(sim), "SNPs - ", max(sim$N), "Samples")
pdf(paste(config["out_sim_prefix"], "_f.pdf", sep=""))
plot(density(sim$f.obs), col="black", lwd=2, ylim=ylim, xlim=xlim, main=main.txt, xlab="F", type="n")
abline(v=0, h=0, col="gray", lty=2)
lines(density(sim$f.obs), col="black", lwd=2)
lines(density(sim$f.sim), col="red")
legend("topright", c("Data", "Simulation"), col=c("black", "red"), lwd=c(2,2))
dev.off()

# MAF vs Pvalue
png(config["out_maf_plot"], width=600, height=600)
par(mar=c(5,5,4,2)+0.1, lwd=1.5,
    cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
plot(hwe$MAF[aut], -log10(hwe$p.value[aut]), xlab="MAF", ylab="-log10(p-value)",
     main="Autosomal SNPs")
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

xy <- IntensityReader(config["samp_xy_file"])
xyData <- IntensityData(xy, snpAnnot=snpAnnot)
geno <- GenotypeReader(config["samp_geno_file"])
genoData <- GenotypeData(geno, snpAnnot=snpAnnot)

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
  par(mfrow=c(3,3), mar=c(5,5,4,2)+0.1, lwd=1.5,
      cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
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
            pbin[i-3],"_%03d.png",sep=""),width=720,height=720)
    par(mfrow=c(3,3), mar=c(5,5,4,2)+0.1, lwd=1.5,
        cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    genoClusterPlot(xyData, genoData, snpID=dat$snpID, main.txt=mtxt)
    dev.off()
  }
}

