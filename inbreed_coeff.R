##########
# Individual inbreeding coefficient
# Usage: R --args config.file < inbreed_coeff.R
##########

library(GWASTools)
library(SNPRelate)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]
print(config.table)

# SNP filter - contain the condition of every two SNPs being at least 15kb apart
# pre-selected SNPs: 15kb apart, autosomal, MAF > 0%, and MCR < 5%
# also select SNPs with MAF > 5%
snp.ids <- getobj(config["out_snp_file"])
length(snp.ids)
afreq <- getobj(config["out_afreq_file"])
snpID <- rownames(afreq)
maf.filt <- !is.na(afreq[,"all"]) & afreq[,"all"] < 0.95 & afreq[,"all"] > 0.05
snp.sel <- intersect(snp.ids, snpID[maf.filt])
length(snp.sel)
afreq.sel <- afreq[as.character(snp.sel), "all"]
stopifnot(length(snp.sel) == length(afreq.sel))
stopifnot(allequal(names(afreq.sel), snp.sel))

gdsobj <- openfn.gds(config["gds_geno_file"])
inbrd.coeff <- snpgdsIndInb(gdsobj, snp.id=snp.sel, allele.freq=afreq.sel)
closefn.gds(gdsobj)

save(inbrd.coeff, file=config["out_inbrd_file"])
