##########
# SNP selection for IBD
# Usage: R --args config.file < ibd_snp_sel.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]

(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)

afreq <- getobj(config["out_afreq_file"])
stopifnot(allequal(rownames(afreq), snpID))

chrom <- getChromosome(snpAnnot)
pos <- getPosition(snpAnnot)
missing <- getVariable(snpAnnot, config["annot_snp_missingCol"])
pool <- chrom < 23 & missing < 0.05 & afreq[,"all"] > 0 & afreq[,"all"] < 1
table(pool)

rsnp <- apartSnpSelection(chrom, pos, min.dist=15000, init.sel=pool)
table(rsnp)

snps.ibd <- snpID[rsnp]
save(snps.ibd, file=config["out_snp_file"])
