##########
# LD pruning
# Usage: R --args config.file < ld_pruning.R
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
required <- c("annot_scan_file", "annot_snp_file", "build", "gds_geno_file", "study_unrelated_file")
optional <- c("annot_snp_rsIDCol", "ld_r_threshold", "ld_win_size", "out_disc_file", "out_pruned_file")
default <- c("rsID", 0.32, 10, NA, "snps_pruned.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# select scans
(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)
scan.sel <- getobj(config["study_unrelated_file"])
stopifnot(all(scan.sel %in% scanID))
length(scan.sel)

# select initial SNPs
(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

# remove discordant SNPs if duplicate discordance has been run
if (!is.na(config["out_disc_file"])) {
  discord <- getobj(config["out_disc_file"])
  disc <- discord$discordance.by.snp
  snp.conc.rsid <- row.names(disc[disc$discordant == 0,])
  rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"])
  snp.sel <- snpID[rsID %in% snp.conc.rsid]
} else {
  message("no discordance file specified; using all SNPs")
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

gdsobj <- openfn.gds(config["gds_geno_file"])
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=win, ld.threshold=r)
closefn.gds(gdsobj)

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned)
save(snp.pruned, file=config["out_pruned_file"])
