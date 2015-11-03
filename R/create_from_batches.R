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
bgds <- openfn.gds(file1)
nsnp <- objdesp.gdsn(index.gdsn(bgds, "snp.id"))$dim

## variables to copy
vars <- ls.gdsn(bgds)
vars.snp <- vars[grepl("^snp", vars)]
vars <- vars[!grepl("^sample", vars) & !grepl("^snp", vars)]

## create new gds file
newfile <- config[paste0(type, "_file")]
gds <- createfn.gds(newfile)

## sample ID (fill in later)
add.gdsn(gds, "sample.id", storage="integer", valdim=0, compress="ZIP_RA")

## add snp variables
for (v in vars.snp) {
    node <- index.gdsn(bgds, v)
    add.gdsn(gds, v, storage=objdesp.gdsn(node)$storage, valdim=0, compress="ZIP_RA")
    newnode <- index.gdsn(gds, v)
    append.gdsn(newnode, node)
    readmode.gdsn(newnode)
}
sync.gds(gds)

## add other variables (but not data yet)
if ("genotype" %in% vars) {
    geno.node <- add.gdsn(gds, "genotype", storage="bit2", valdim=c(nsnp, 0))
    put.attr.gdsn(geno.node, "snp.order")
}
for (v in setdiff(vars, "genotype")) {
    node <- index.gdsn(bgds, v)
    add.gdsn(gds, v, storage=objdesp.gdsn(node)$storage, valdim=c(nsnp, 0), compress="ZIP_RA.max:8M")
}
sync.gds(gds)

closefn.gds(bgds)


## for each batch, read and append
for (b in 1:nbatches) {
    message("Copying batch ", b)
    bfile <- filename(type, b)
    bgds <- openfn.gds(bfile)

    ## add sample.id
    samp.node <- index.gdsn(bgds, "sample.id")
    append.gdsn(index.gdsn(gds, "sample.id"), samp.node)
    
    ## add other variables
    for (v in vars) {
        append.gdsn(index.gdsn(gds, v), index.gdsn(bgds, v))
    }
    closefn.gds(bgds)
    sync.gds(gds)
}

## close
GWASTools:::.close.gds.class(gds, verbose=TRUE)


## for each batch, compare both files
gds <- GdsReader(newfile)
sample.id <- getVariable(gds, "sample.id")

i <- 0
for (b in 1:nbatches) {
    message("Checking batch ", b)
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
    i <- i + nsamp
    close(bgds)
}

close(gds)

message("remove the following files after successful completion:\n",
        paste(filename(type, 1:nbatches), collapse="\n"))
