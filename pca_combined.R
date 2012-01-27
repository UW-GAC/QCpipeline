##########
# PCA - combined with external data
# Usage: R --args config.file < pca_combined.R
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

# pruned SNPs
snps.pruned <- getobj(config["out_pruned_file"])
snpAnnot <- getobj(config["annot_snp_file"])
snpID <- getSnpID(snpAnnot)
rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"])
snp.ids <- rsID[snpID %in% snps.pruned]
length(snp.ids)

# select unduplicated scans from study
study.scan.ids <- getobj(config["study_unduplicated_file"])
length(study.scan.ids)

# select unrelated scans from external
ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
ext.unrel <- getVariable(ext.scanAnnot, config["ext_annot_scan_unrelCol"])
stopifnot(is.logical(ext.unrel))
ext.scan.ids <- getScanID(ext.scanAnnot)[ext.unrel]
length(ext.scan.ids)

comb.scan.ids <- c(study.scan.ids, ext.scan.ids)

# remove duplicate scans between study and external
discord <- getobj(config["out_disc_file"])
disc <- discord$discordance.by.subject
dups <- vector()
for (i in 1:length(disc)) {
  all.scans <- unlist(dimnames(disc[[i]]))
  sel.scan <- intersect(study.scan.ids, all.scans)
  if (length(sel.scan) == 0) sel.scan <- all.scans[1]
  dups <- append(dups, setdiff(all.scans, sel.scan))
}
dups <- as.integer(dups)

undup.scans <- setdiff(comb.scan.ids, dups)
length(undup.scans)

gdsobj <- openfn.gds(config["out_comb_gds_geno_file"])
pca <- snpgdsPCA(gdsobj, sample.id=undup.scans, snp.id=snp.ids)
save(pca, file=config["out_pca_file"])

pca.corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:8)
save(pca.corr, file=config["out_corr_file"])

closefn.gds(gdsobj)
