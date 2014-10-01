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
optional <- c("out_inbrd_file", "scan_ibd_include_file", "ibd_method")
default <- c("inbreed_coeff.RData", NA, "MoM")
config <- setConfigDefaults(config, required, optional, default)
print(config)

snp.sel <- getobj(config["out_snp_file"])
length(snp.sel)

if (!is.na(config["scan_ibd_include_file"])) {
  scan.sel <- getobj(config["scan_ibd_include_file"])
} else {
  scan.sel <- NULL
}
length(scan.sel)

if (config["ibd_method"] == "MoM") {
  method <- "mom.weir"
} else if (config["ibd_method"] == "MLE") {
  method <- "mle"
}
method

gdsobj <- snpgdsOpen(config["gds_geno_file"])
inbrd.coeff <- snpgdsIndInb(gdsobj, snp.id=snp.sel, sample.id=scan.sel, method=method)
snpgdsClose(gdsobj)

save(inbrd.coeff, file=config["out_inbrd_file"])
