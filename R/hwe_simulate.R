##########
# HWE simulation of inbreeding coefficient
# Usage: R --args config.file < hwe_simulate.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "scan_include_file")
optional <- c("out_hwe_prefix", "nsnp_simulate", "out_sim_prefix")
default <- c("hwe", 50000, "hwe_sim")
config <- setConfigDefaults(config, required, optional, default)
print(config)

hwe <- getobj(paste(config["out_hwe_prefix"], "RData", sep="."))

# autosomes only that are not monomorphic
aut <- hwe[hwe$chr %in% 1:22 & !is.na(hwe$f), ]
dim(aut)

(nsnp <- min(as.integer(config["nsnp_simulate"]), nrow(aut)))
(nsamp <- max(aut$nAA + aut$nAB + aut$nBB))

# random sample of snps
aut.random <- aut[sample(1:nsnp), ]
aut.random$afreq <- ifelse(aut.random$minor.allele == "A", aut.random$MAF, 1-aut.random$MAF)

geno <- hweSimulateGenotypeMatrix(nsnp, nsamp, aut.random$afreq)

res <- hweSimulate(geno, aut.random$afreq)

res$snpID <- aut.random$snpID
res$pval.obs <- aut.random$pval
res$f.obs <- aut.random$f
res$N <- nsamp

save(res, file=paste0(config["out_sim_prefix"], ".RData"))
