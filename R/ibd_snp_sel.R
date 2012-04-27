##########
# SNP selection for IBD
# Usage: R --args config.file < ibd_snp_sel.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file", "out_afreq_file")
optional <- c("annot_snp_missingCol", "out_snp_file")
default <- c("missing.n1", "ibd_snp_sel.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

afreq <- getobj(config["out_afreq_file"])
stopifnot(allequal(rownames(afreq), snpID))

chrom <- getChromosome(snpAnnot)
pos <- getPosition(snpAnnot)
missing <- getVariable(snpAnnot, config["annot_snp_missingCol"])
maf.filt <- !is.na(afreq[,"all"]) & afreq[,"all"] > 0 & afreq[,"all"] < 1
pool <- chrom < 23 & missing < 0.05 & maf.filt
table(pool)

rsnp <- apartSnpSelection(chrom, pos, min.dist=15000, init.sel=pool)
table(rsnp)

snps.ibd <- snpID[rsnp]
save(snps.ibd, file=config["out_snp_file"])
