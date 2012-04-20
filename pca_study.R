##########
# PCA - unrelated study samples
# Usage: R --args config.file < pca_study.R
##########

library(GWASTools)
library(QCpipeline)
library(SNPRelate)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("gds_geno_file", "out_pruned_file", "study_unrelated_file")
optional <- c("num_evs_to_plot", "out_corr_file", "out_pca_file")
default <- c(12, "pca_corr.RData", "pca.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

snp.ids <- getobj(config["out_pruned_file"])
length(snp.ids)

scan.ids <- getobj(config["study_unrelated_file"])
length(scan.ids)

gdsobj <- openfn.gds(config["gds_geno_file"])
pca <- snpgdsPCA(gdsobj, sample.id=scan.ids, snp.id=snp.ids)
save(pca, file=config["out_pca_file"])

nev <- as.integer(config["num_evs_to_plot"])
pca.corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:nev)
save(pca.corr, file=config["out_corr_file"])

closefn.gds(gdsobj)
