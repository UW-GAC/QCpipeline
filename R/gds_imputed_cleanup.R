##########
# Imputed dosage GDS file
# Usage: R --args config.file < gds_imputed_cleanup.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()


# read configuration
args <- commandArgs(trailingOnly=TRUE)
print(args)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("orig_annot_scan_file",
              "impute2_sample_prefix",
              "impute2_geno_prefix",
              "metrics_path",
              "snp_segment_mapping_file",
              "out_gds_dir",
              "out_gds_prefix")
optional <- c()
default <- c()
config <- setConfigDefaults(config, required, optional, default)
print(config)

files <- list.files(config["out_gds_dir"], full.names=TRUE)

# check the scan annotations
(files_scan <- files[grepl("scanAnnot.tmp.RData", files)])

# read in the first one
scanAnnot <- getobj(files_scan[1])

## changes to imputedDosageFile make this not feasible
# # loop through and read in the others
# for (f in files_scan[2:length(files_scan)]){
#   print(paste("Checking file:", f))
#   scanAnnot.chk <- getobj(f)
#   for (var in varLabels(scanAnnot)[varLabels(scanAnnot) != "added"]) stopifnot(allequal(scanAnnot[[var]], scanAnnot.chk[[var]]))
# }

scanAnnot$added <- NULL

scan.file <- file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_scanAnnot.RData", sep=""))
save(scanAnnot, file=scan.file)

# remove the temporary files
file.remove(files_scan)


# check the snp annotations against the mapping file
files_snp <- files[grepl("snpAnnot.RData", files)]

snp.segment.map <- read.csv(config["snp_segment_mapping_file"], stringsAsFactors=FALSE, header=TRUE, as.is=TRUE)

keeplist <- list()

for (i in 1:length(files_snp)){
  snpAnnot <- getobj(files_snp[i])
  keeplist[[i]] <- snp.segment.map[snp.segment.map$chrom %in% snpAnnot$chromosome & snp.segment.map$segment %in% snpAnnot$segment, ]
}

snp.segment.keep <- do.call(rbind, keeplist)

# copy the snp-segment mapping file so we have a record.. I think we will need it later.
write.csv(snp.segment.keep, file=file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_snp_segment_map.csv", sep="")), quote=FALSE, row.names=FALSE)


