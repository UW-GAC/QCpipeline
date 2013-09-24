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

# check config and set defaults
required <- c("annot_scan_file", "ext_annot_scan_file", "out_comb_annot_scan_file",
              "out_comb_gds_geno_file", "out_pruned_file",
              "study_unduplicated_file")
optional <- c("annot_scan_hapmapCol", "annot_scan_subjectCol", "annot_scan_unrelCol",
              "ext_annot_scan_unrelCol", "num_evs_to_plot", "out_corr_file",
              "out_disc_file", "out_pca_file")
default <- c("geno.cntl", "subjectID", "unrelated",
             "unrelated", 12, "pca_combined_corr.RData", NA, "pca_combined.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# multithreading on pearson?
nSlots <- Sys.getenv("NSLOTS")
nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
print(paste("Running with", nThreads,"thread(s)."))

# pruned SNPs
snp.ids <- getobj(config["out_pruned_file"])
length(snp.ids)

# excluded scans?
comb.scanAnnot <- getobj(config["out_comb_annot_scan_file"])
scanAnnot <- getobj(config["annot_scan_file"])
scanAnnot <- scanAnnot[scanAnnot$scanID %in% comb.scanAnnot$scanID,]
ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
ext.scanAnnot <- ext.scanAnnot[ext.scanAnnot$scanID %in% comb.scanAnnot$scanID,]

# select unduplicated scans from study
study.scan.ids <- getobj(config["study_unduplicated_file"])
length(study.scan.ids)

# check
scanID <- getScanID(scanAnnot)
stopifnot(all(study.scan.ids %in% scanID))
subjID <- getVariable(scanAnnot, config["annot_scan_subjectCol"])
if (any(duplicated(subjID[scanID %in% study.scan.ids]))) {
  stop("study_unduplicated_file contains some duplicate scans")
}

# select unrelated hapmaps
unrel <- getVariable(scanAnnot, config["annot_scan_unrelCol"])
stopifnot(is.logical(unrel))
hapmap <- getVariable(scanAnnot, config["annot_scan_hapmapCol"])
hm.unrel.ids <- scanID[hapmap == 1 & unrel]
length(hm.unrel.ids)

# select unrelated scans from external
ext.unrel <- getVariable(ext.scanAnnot, config["ext_annot_scan_unrelCol"])
stopifnot(is.logical(ext.unrel))
ext.scan.ids <- getScanID(ext.scanAnnot)[ext.unrel]
length(ext.scan.ids)

comb.scan.ids <- unique(c(study.scan.ids, hm.unrel.ids, ext.scan.ids))
length(comb.scan.ids)

# remove duplicate scans between study and external (if any)
if (!is.na(config["out_disc_file"])) {
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
} else {
  undup.scans <- comb.scan.ids
}

gdsobj <- openfn.gds(config["out_comb_gds_geno_file"])
pca <- snpgdsPCA(gdsobj, sample.id=undup.scans, snp.id=snp.ids, num.thread=nThreads)
save(pca, file=config["out_pca_file"])

nev <- as.integer(config["num_evs_to_plot"])
pca.corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:nev, num.thread=nThreads)
save(pca.corr, file=config["out_corr_file"])

closefn.gds(gdsobj)
