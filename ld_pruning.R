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

discord <- getobj(config["out_disc_file"])
disc <- discord$discordance.by.snp
snp.conc.rsid <- row.names(disc[disc$discordant == 0,])
rsID <- getVariable(snpAnnot, config["annot_snp_rsIDCol"])
snp.sel <- snpID[rsID %in% snp.conc.rsid]
length(snp.sel)

gdsobj <- openfn.gds(config["gds_geno_file"])
snpset <- snpgdsLDpruning(gdsobj, sample.id=scan.sel, snp.id=snp.sel,
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10e6,
                          ld.threshold=as.numeric(config["ld_r_threshold"]))
closefn.gds(gdsobj)

snp.pruned <- unlist(snpset, use.names=FALSE)
length(snp.pruned)
save(snp.pruned, file=config["out_pruned_file"])
