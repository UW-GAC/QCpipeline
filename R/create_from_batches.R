##########
# Combine GDS files created in batches
# Usage: R --args config.file nbatches type < create_from_batches.R
##########

library(GWASTools)
library(QCpipeline)
library(gdsfmt)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("geno_file", "xy_file", "bl_file")
optional <- c()
default <- c()
config <- setConfigDefaults(config, required, optional, default)
print(config)

nbatches <- as.integer(args[2])

type <- args[3]
stopifnot(type %in% c("geno", "xy", "bl"))

batchfile <- function(x, batch) {
    file.path(dirname(x), paste0("batch", batch, ".", basename(x)))
}

## select files based on type
filename <- function(type, batch) {
    batchfile(config[paste0(type, "_file")], batch)
}

## open first file
file1 <- filename(type, 1)
gds <- openfn.gds(file1, readonly=FALSE)
n1 <- objdesp.gdsn(index.gdsn(gds, "sample.id"))$dim

## variables to copy
vars <- ls.gdsn(gds)
vars <- vars[!grepl("^sample", vars) & !grepl("^snp", vars)]

## decompress data for writing
vars.compress <- c("sample.id",  setdiff(vars, "genotype"))
for (v in vars.compress) compression.gdsn(index.gdsn(gds, v), "")
## must close and reopen
closefn.gds(gds)
gds <- openfn.gds(file1, readonly=FALSE)

## for each subsequent file, read and append to first file
for (b in 2:nbatches) {
    bfile <- filename(type, b)
    bgds <- GdsReader(bfile)

    # get data
    nsamp <- getDimension(bgds, "sample.id")
    for (s in 1:nsamp) {
        sample.id <- getVariable(bgds, "sample.id", start=s, count=1)
        dat <- lapply(setNames(vars, vars),
                      function(v) getVariable(bgds, v, start=c(1,s), count=c(-1,1)))
        GWASTools:::.addData.gds.class(gds, vars, dat, sample.id)
    }
    close(bgds)
    sync.gds(gds)
}

## compress and close
for (v in vars.compress) compression.gdsn(index.gdsn(gds, v), "ZIP.max")
GWASTools:::.close.gds.class(gds, verbose=TRUE)


## for each batch (starting with second), compare both files
gds <- GdsReader(file1)
sample.id <- getVariable(gds, "sample.id")

i <- n1
for (b in 2:nbatches) {
    bfile <- filename(type, b)
    bgds <- GdsReader(bfile)

    bsamp <- getVariable(bgds, "sample.id")
    nsamp <- length(bsamp)
    stopifnot(allequal(bsamp, sample.id[(i+1):(i+nsamp)]))
    for (s in 1:nsamp) {
        for (v in vars) {
            d1 <- getVariable(gds, v, start=c(1,i+s), count=c(-1,1))
            db <- getVariable(bgds, v, start=c(1,s), count=c(-1,1))
            stopifnot(allequal(d1, db))
        }
    }
    i <- n1 + nsamp
    close(bgds)
}

close(gds)

file.rename(file1, config[paste0(type, "_file")])
for (b in 2:nbatches) file.remove(filename(type, b))
