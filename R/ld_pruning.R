##########
# LD pruning
# Usage: R --args config.file pca.type < ld_pruning.R
##########

library(GWASTools)
library(QCpipeline)
library(SNPRelate)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check for type
if (length(args) < 2) stop("missing pca type (study or combined)")
type <- args[2]

# check config and set defaults
if (type == "study") {
    required <- c("annot_scan_file", "annot_snp_file", "build", "gds_geno_file", "study_unrelated_file")
    scanfile <- config["annot_scan_file"]
    snpfile <- config["annot_snp_file"]
    genofile <- config["gds_geno_file"]
} else if (type == "combined") {
    required <- c("out_comb_prefix", "build", "study_unrelated_file")
    scanfile <- paste0(config["out_comb_prefix"], "_scanAnnot.RData")
    snpfile <- paste0(config["out_comb_prefix"], "_snpAnnot.RData")
    genofile <- paste0(config["out_comb_prefix"], ".gds")
}
optional <- c("ld_r_threshold", "ld_win_size",
              "snp_pruning_include_file", "out_pruned_file")
default <- c(0.32, 10, NA, "snps_pruned.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)


# multithreading on pearson?
# num.thread is not implemented in LD pruning, but uncomment these when it is
##nSlots <- Sys.getenv("NSLOTS")
##nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
##print(paste("Running with", nThreads,"thread(s)."))
nThreads <- 1


# select scans
(scanAnnot <- getobj(scanfile))
scanID <- getScanID(scanAnnot)
scan.sel <- getobj(config["study_unrelated_file"])
stopifnot(all(scan.sel %in% scanID))
length(scan.sel)

# select initial SNPs
(snpAnnot <- getobj(snpfile))
snpID <- getSnpID(snpAnnot)
if (!is.na(config["snp_pruning_include_file"])) {
  snp.sel <- intersect(snpID, getobj(config["snp_pruning_include_file"]))
} else {
  snp.sel <- snpID
}
length(snp.sel)

# remove SNPs in spike regions
filt <- get(data(list=paste("pcaSnpFilters", config["build"], sep=".")))
chrom <- getChromosome(snpAnnot)
pos <- getPosition(snpAnnot)
snp.filt <- rep(TRUE, length(snpID))
for (f in 1:nrow(filt)) {
  snp.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos & pos < filt$end.base[f]] <- FALSE
}
snp.sel <- intersect(snp.sel, snpID[snp.filt])
length(snp.sel)

r <- as.numeric(config["ld_r_threshold"])
win <- as.numeric(config["ld_win_size"]) * 1e6

gdsobj <- snpgdsOpen(genofile)
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.02,
                          method="corr", slide.max.bp=win, ld.threshold=r,
                          num.thread=nThreads)
snpgdsClose(gdsobj)

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned)
save(snp.pruned, file=config["out_pruned_file"])
