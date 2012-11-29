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

# check config and set defaults
required <- c("gds_geno_file", "out_snp_file")
optional <- c("out_inbrd_file")
default <- c("inbreed_coeff.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

snp.sel <- getobj(config["out_snp_file"])
length(snp.sel)

gdsobj <- openfn.gds(config["gds_geno_file"])
inbrd.coeff <- snpgdsIndInb(gdsobj, snp.id=snp.sel)
closefn.gds(gdsobj)

save(inbrd.coeff, file=config["out_inbrd_file"])
