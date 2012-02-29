##########
# Individual inbreeding coefficient
# Usage: R --args config.file < inbreed_coeff.R
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

# SNP filter - contain the condition of every two SNPs being at least 15kb apart
# pre-selected SNPs: 15kb apart, autosomal, MAF > 0%, and MCR < 5%
snp.sel <- getobj(config["out_snp_file"])
length(snp.sel)
afreq <- getobj(config["out_afreq_file"])
afreq.sel <- afreq[as.character(snp.sel), "all"]
stopifnot(length(snp.sel) == length(afreq.sel))
stopifnot(allequal(names(afreq.sel), snp.sel))

gdsobj <- openfn.gds(config["gds_geno_file"])
inbrd.coeff <- snpgdsIndInb(gdsobj, snp.id=snp.sel, allele.freq=afreq.sel)
closefn.gds(gdsobj)

save(inbrd.coeff, file=config["out_inbrd_file"])
