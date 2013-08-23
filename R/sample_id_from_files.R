##########
## Get sample ID in each FinalReport file
## Usage: R --args config <outfile> < sample_id_from_files.R
##########

library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

basedir <- config["raw_path"]
sep <- config["raw_sepType"]
skip <- as.integer(config["raw_skipNum"])
samp.col <- as.integer(config["raw_sampleCol"])

if (length(args) > 1) {
  outfile <- args[2]
} else {
  outfile <- "sample_file_map.RData"
}

files <- list.files(basedir)
length(files)

## for each file, extract sample ID
samples <- character(length(files))
for (i in 1:length(files)) {
  dat <- read.table(file.path(basedir, files[i]),
                    sep=sep, skip=skip, header=FALSE,
                    nrow=1, as.is=TRUE)
  samples[i] <- dat[,samp.col]
} 

file.df <- data.frame(samples, file=files, stringsAsFactors=FALSE)
names(file.df)[1] <- config["annot_scan_nameCol"]
save(file.df, file=outfile)
