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

# which SNPs for the correlation?
corr.snp <- NULL
if (!is.na(config["corr_only_pruned_snps"])) {
  if (as.logical(config["corr_only_pruned_snps"])) {
    corr.snp <- snp.ids
  }
}

nev <- as.integer(config["num_evs_to_plot"])
pca.corr <- snpgdsPCACorr(pca, gdsobj, snp.id=corr.snp, eig.which=1:nev)
save(pca.corr, file=config["out_corr_file"])

closefn.gds(gdsobj)
