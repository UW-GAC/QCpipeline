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
required <- c("annot_scan_file", "ext_annot_scan_file", "out_comb_prefix",
              "out_pruned_file", "study_unduplicated_file")
optional <- c("annot_scan_hapmapCol", "annot_scan_subjectCol", "annot_scan_unrelCol",
              "ext_annot_scan_unrelCol", "num_evs_to_plot", "out_corr_file",
              "out_pca_file", "include_study_hapmaps")
default <- c("geno.cntl", "subjectID", "unrelated",
             "unrelated", 12, "pca_combined_corr.RData", "pca_combined.RData",
             TRUE)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# multithreading on pearson?
nSlots <- Sys.getenv("NSLOTS")
nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
## allow running outside the cluster environment
if (nThreads == 0) nThreads <- 1
print(paste("Running with", nThreads,"thread(s)."))

# pruned SNPs
snp.ids <- getobj(config["out_pruned_file"])
length(snp.ids)

# excluded scans?
comb.scanAnnot <- getobj(paste0(config["out_comb_prefix"], "_scanAnnot.RData"))
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

## might not want study HapMaps if they are related to external samples
if (as.logical(config["include_study_hapmaps"])) {
    comb.scan.ids <- unique(c(study.scan.ids, hm.unrel.ids, ext.scan.ids))
} else {
    comb.scan.ids <- unique(c(study.scan.ids, ext.scan.ids))
}
length(comb.scan.ids)

# remove duplicate scans between study and external (if any)
dups <- comb.scanAnnot$scanID[duplicated(comb.scanAnnot$subjectID)]
undup.scans <- setdiff(comb.scan.ids, dups)
length(undup.scans)

gdsobj <- snpgdsOpen(paste0(config["out_comb_prefix"], ".gds"))
pca <- snpgdsPCA(gdsobj, sample.id=undup.scans, snp.id=snp.ids, num.thread=nThreads)
save(pca, file=config["out_pca_file"])

nev <- as.integer(config["num_evs_to_plot"])
pca.corr <- snpgdsPCACorr(pca, gdsobj, eig.which=1:nev, num.thread=nThreads)
save(pca.corr, file=config["out_corr_file"])

snpgdsClose(gdsobj)
