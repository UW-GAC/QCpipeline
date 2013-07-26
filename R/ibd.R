##########
# IBD
# Usage: R --args config.file < ibd.R
##########

library(GWASTools)
#library(QCpipeline)
library(QCpipeline)
library(SNPRelate)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "gds_geno_file", "out_snp_file")
optional <- c("annot_scan_familyCol", "ibd_method", "out_ibd_file",
              "out_ibd_kc32_file", "scan_ibd_include_file")
default <- c(NA, "MoM", "ibd.RData", "ibd_kc32.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# multithreading on pearson?
nSlots <- Sys.getenv("NSLOTS")
nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
print(paste("Running with", nThreads,"thread(s)."))


snp.ids <- getobj(config["out_snp_file"])
length(snp.ids)

if (!is.na(config["scan_ibd_include_file"])) {
  scan.ids <- getobj(config["scan_ibd_include_file"])
} else {
  # if scan.ids=NULL for IBD functions, all scans are used
  scan.ids <- NULL
}
length(scan.ids)

if (!is.na(config["annot_scan_familyCol"])) {
  scanAnnot <- getobj(config["annot_scan_file"])
  family <- getVariable(scanAnnot, config["annot_scan_familyCol"])
  if (!is.null(scan.ids)) {
    family <- family[getScanID(scanAnnot) %in% scan.ids]
  }
} else {
  family <- NULL
}


gdsobj <- openfn.gds(config["gds_geno_file"])
if (config["ibd_method"] == "MoM") {
  ibd <- snpgdsIBDMoM(gdsobj, sample.id=scan.ids, snp.id=snp.ids, num.thread=nThreads)
} else if (config["ibd_method"] == "MLE") {
  ibd <- snpgdsIBDMLE(gdsobj, sample.id=scan.ids, snp.id=snp.ids, num.thread=nThreads,
                      method="EM")
} else if (config["ibd_method"] == "KING") {
  ibd <- snpgdsIBDKING(gdsobj, sample.id=scan.ids, snp.id=snp.ids, num.thread=nThreads,
                       family.id=family)
} else {
  stop("ibd method not recognized")
}
closefn.gds(gdsobj)

save(ibd, file=config["out_ibd_file"])

# save only pairs where KC > 1/32
ibd <- snpgdsIBDSelection(ibd, kinship.cutoff=1/32)
save(ibd, file=config["out_ibd_kc32_file"])
