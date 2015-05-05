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

# check config and set defaults
required <- c("raw_path")
optional <- c("annot_scan_fileCol", "annot_scan_nameCol",
              "raw_scanNameInFile", "raw_sampleRow",
              "raw_sepType", "raw_skipNum", "raw_sampleCol")
default <- c("file", "Sample.Name", 0, 9, ",", 11, 9)
config <- setConfigDefaults(config, required, optional, default)
print(config)

basedir <- config["raw_path"]
sep <- config["raw_sepType"]
skip <- as.integer(config["raw_skipNum"])
samp.col <- as.integer(config["raw_sampleCol"])

if (length(args) > 1) {
  outfile <- args[2]
} else {
  outfile <- "sample_file_map.RData"
}

files <- list.files(basedir, full.names=TRUE)
length(files)

## for each file, extract sample ID
samples <- character(length(files))
for (i in 1:length(files)) {
    if (config["raw_scanNameInFile"] == 1) {
        dat <- read.table(files[i],
                          sep=sep, skip=skip, header=FALSE,
                          nrow=1, as.is=TRUE)
        samples[i] <- dat[,samp.col]
    } else {
        skip <- as.integer(config["raw_sampleRow"]) - 1
        
        f <- file(files[i], "r")
        x <- readLines(f, n=skip)
        x <- readLines(f, n=1)
        close(f)
        samples[i] <- strsplit(x, ",")[[1]][2]
    }
} 

file.df <- data.frame(samples, file=files, stringsAsFactors=FALSE)
names(file.df)[1] <- config["annot_scan_nameCol"]
save(file.df, file=outfile)
