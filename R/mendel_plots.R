##########
# Cluster plots for Mendelian errors
# Usage: R --args config.file < mendel_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "mend.bin.end", "mend.bin.start", "geno_file", "xy_file")
optional <- c("annot_snp_rsIDCol", "out_mend_clust_prefix", "out_mend_file")
default <- c("rsID", "mendel_clust", "mendel_err.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

men <- getobj(config["out_mend_file"])
men <- men$snp
stopifnot(names(men) == c("snpID", "error.cnt", "check.cnt"))
names(men) <- c("snpID", "err", "chk")
men$err[men$chk == 0] <- NA

snpAnnot <- getobj(config["annot_snp_file"])
snp <- getVariable(snpAnnot, c("snpID", config["annot_snp_rsIDCol"], "chromosome"))
names(snp) <- c("snpID", "rsID", "chromosome")
snp <- merge(snp, men)

bin.start <- config["mend.bin.start"]
bin.start <- as.integer(unlist(strsplit(bin.start, " ", fixed=TRUE), use.names=FALSE))
bin.end <- config["mend.bin.end"]
bin.end <- as.integer(unlist(strsplit(bin.end, " ", fixed=TRUE), use.names=FALSE))
stopifnot(length(bin.start) == length(bin.end))

bins <- rep(NA, (length(bin.start)))
ids <- list()
for(i in 1:(length(bin.start))) {
  ids[[i]] <- snp$snpID[snp$err %in% bin.start[i]:bin.end[i]]
  
  # number of snps in each bin
  bins[i] <- length(ids[[i]]) 
}
sum(bins) 
sum(unlist(lapply(ids,length)))
bins

xy <- IntensityReader(config["xy_file"])
xyData <- IntensityData(xy, snpAnnot=snpAnnot)
geno <- GenotypeReader(config["geno_file"])
genoData <- GenotypeData(geno, snpAnnot=snpAnnot)

for(j in 1:3){  # independent sets of samples, each in a different png file
  # sample n from each bin
  n <- 9
  sids <- lapply(ids, function(x) {if (length(x) > n) sample(x, n) else x}) 
  sids <- unlist(sids)
  # retrieve all sampled data
  dat <- snp[snp$snpID %in% sids,]
  dat <- dat[order(dat$err),]
  
  mtxt <- paste("Chr",dat$chromosome,dat$rsID, "\n", dat$err,"Mendelian error(s)")
  png(paste(config["out_mend_clust_prefix"], "_",
            j,"%03d.png",sep=""),width=720,height=720)
  par(mfrow=c(3,3), mar=c(5,5,4,2)+0.1, lwd=1.5,
      cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
  genoClusterPlot(xyData, genoData, snpID=dat$snpID, main.txt=mtxt)
  dev.off()
}
