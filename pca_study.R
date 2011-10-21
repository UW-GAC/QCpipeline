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
print(config)

snp.ids <- getobj(config["out_pruned_file"])
length(snp.ids)

scan.ids <- getobj(config["study_unrelated_file"])
length(scan.ids)

gdsobj <- openfn.gds(config["gds_geno_file"])
pca <- snpgdsPCA(gdsobj, sample.id=scan.ids, snp.id=snp.ids)
save(pca, file=config["out_pca_file"])

pca.corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:8)
save(pca.corr, file=config["out_corr_file"])

closefn.gds(gdsobj)
