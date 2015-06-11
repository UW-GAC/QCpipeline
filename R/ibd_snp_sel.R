##########
# SNP selection for IBD
# Usage: R --args config.file < ibd_snp_sel.R
##########

library(GWASTools)
library(SNPRelate)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("gds_geno_file")
optional <- c("ld_r_threshold", "ld_win_size", "maf_threshold",
              "out_snp_file", "scan_pruning_include_file",
              "snp_pruning_include_file")
default <- c(0.32, 10, 0.05, "ibd_snp_sel.RData", NA, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)


# multithreading on pearson?
# num.thread is not implemented in LD pruning, but uncomment these when it is
##nSlots <- Sys.getenv("NSLOTS")
##nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
##print(paste("Running with", nThreads,"thread(s)."))
nThreads <- 1


if (!is.na(config["scan_pruning_include_file"])) {
  scan.ids <- getobj(config["scan_pruning_include_file"])
} else {
  # if scan.ids=NULL for IBD functions, all scans are used
  scan.ids <- NULL
}
length(scan.ids)

if (!is.na(config["snp_pruning_include_file"])) {
  snp.ids <- getobj(config["snp_pruning_include_file"])
} else {
  # if snp.ids=NULL for IBD functions, all snps are used
  snp.ids <- NULL
}
length(snp.ids)

maf <- as.numeric(config["maf_threshold"])
r <- as.numeric(config["ld_r_threshold"])
win <- as.numeric(config["ld_win_size"]) * 1e6

gdsobj <- snpgdsOpen(config["gds_geno_file"])
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.ids, snp.id=snp.ids,
                          autosome.only=TRUE, maf=maf, missing.rate=0.02,
                          method="corr", slide.max.bp=win, ld.threshold=r, 
                          num.thread=nThreads)
snpgdsClose(gdsobj)

snps.ibd <- unlist(snpset, use.names=FALSE)
save(snps.ibd, file=config["out_snp_file"])
