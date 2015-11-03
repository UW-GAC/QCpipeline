##########
# Imputed dosage GDS file
# Usage: R --args config.file chromosome <test> < gds_create_imputed.R
##########

library(GWASTools)
library(QCpipeline)
library(gdsfmt)
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
optional <- c("quality_name", "quality_minimum",
              "scan_annot_familyCol", "scan_annot_individCol",
              "logfile_prefix", "gds_type")
default <- c("info", 0,
             "family", "subjectID",
             "log", "dosage")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# check for test
if (length(args) < 2) stop ("missing chromosome ranges")
chromosomes <- eval(parse(text=args[2]))
  
# check for chromosome
if (length(args) < 3) stop("missing chromosome")
chromosome <- as.integer(args[3])
print(chromosome)

# make sure it's a valid chromosome - only 1-23 are imputed.
stopifnot(chromosome %in% 1:23)

scanAnnot.orig <- getobj(config["orig_annot_scan_file"])

# check for test
#if (length(args) > 3) {
#  if ((args[3]) == "test") {
#    message("testing first 5 scans")
#    scanAnnot.orig <- scanAnnot.orig[1:5,] # testing
#  }
#}
nsamp <- nrow(scanAnnot.orig)

# check for impute files
samp.file <- paste(config["impute2_sample_prefix"], "_chr", chromosome, ".sample.gz", sep="")
stopifnot(file.exists(samp.file))

impute.file <- paste(config["impute2_geno_prefix"], "_chr", chromosome, ".gprobs.gz", sep="")
stopifnot(file.exists(impute.file))

# set up input scan data frame
samples <- read.table(samp.file, header=TRUE, as.is=TRUE)
samples <- samples[-1, ] # first row is useless.

# subset for samples included in imputation
scanAnnot.orig <- scanAnnot.orig[scanAnnot.orig$subj.plink & scanAnnot.orig$geno.cntl == 0, ]

# construct the scan.df data frame
x <- apply(getVariable(scanAnnot.orig, c(config["scan_annot_familyCol"], config["scan_annot_individCol"])), 1, paste, collapse=" ")
y <- getScanID(scanAnnot.orig)
scan.df <- data.frame(scanID=y, sampleID=x, stringsAsFactors=FALSE)

# check for quality metric threshold
filename <- paste(config["metrics_path"], "_chr", chromosome, ".metrics.gz", sep="")
stopifnot(file.exists(filename))

# read in the metric for that chromosome
metrics <- read.table(filename, header=TRUE, comment.char="")
snp.exclude.logical <- metrics[[config["quality_name"]]] < config["quality_minimum"] # logical version is needed later.
snp.exclude <- which(snp.exclude.logical)
print(length(snp.exclude))

# starting snp ID -- start at 1 but increment by previous chromosomes
snp.id.start <- 1

# add # of snps included on previous chromosomes to snp.id.start

for (chr in chromosomes[chromosomes < chromosome]){
  fname <- paste(config["metrics_path"], "_chr", chr, ".metrics.gz", sep="")
  ## for troubleshooting
  #message("\treading in: ", fname)
  #if (!file.exists(fname)) next
  tmp <- read.table(fname, header=TRUE, comment.char="")
  snp.id.start <- snp.id.start + sum(tmp[[config["quality_name"]]] >= config["quality_minimum"])
}
message("Starting at snpID ", snp.id.start)

# set up input for gdsImputedDosage
input.files <- c(impute.file, samp.file)

gds.filename <- file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, ".gds", sep=""))
snp.filename.tmp <- file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, "_snpAnnot_tmp.RData", sep=""))
scan.filename <- file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, "_scanAnnot.tmp.RData", sep=""))

genotypeDim <- "scan,snp"

# write to the scratch disk of each node
gds.filename.tmp <- tempfile()

message("will write ",config["gds_type"], " gds file")
message("gds temporarily located at ", gds.filename.tmp)

imputedDosageFile(input.files, gds.filename.tmp, chromosome, input.type="IMPUTE2", input.dosage=FALSE, output.type=config["gds_type"],
                  file.type="gds",
                  snp.annot.filename=snp.filename.tmp, scan.annot.filename=scan.filename,
                  verbose=TRUE, genotypeDim=genotypeDim, scan.df=scan.df, snp.exclude=snp.exclude,
                  snp.id.start=snp.id.start)

# compress the genotypes
gds <- openfn.gds(gds.filename.tmp, readonly=F)
compression.gdsn(index.gdsn(gds, "genotype"), compress="ZIP_RA.max:8M")
closefn.gds(gds)

cleanup.gds(gds.filename.tmp, verbose=T)

# copy it
file.copy(gds.filename.tmp, gds.filename)
# remove the tmp file
file.remove(gds.filename.tmp)

# open it
scanAnnot <- getobj(scan.filename)
snpAnnot <- getobj(snp.filename.tmp)
geno.gds <- GdsGenotypeReader(gds.filename)
genoData <- GenotypeData(geno.gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)

# logfile for snp and sample IDs of all missing data in gds
if(!is.na(config["logfile_prefix"])) {
	message("writing out snp and sample IDs with missing ", config["gds_type"], " data") 
	logfile <- file.path(paste(config["logfile_prefix"], "_chr-", chromosome, "_log.txt", sep=""))
} else if (is.na(config["logfile_prefix"])) {
	logfile <- NULL
	message("NOT writing out log file of all snp and sample IDs with missing ", config["gds_type"], "data")
}
	
res <- checkImputedDosageFile(genoData, snpAnnot, scanAnnot, 
                       input.files=input.files, chromosome,
                       input.type="IMPUTE2", 
					   output.type=config["gds_type"],
                       input.dosage=FALSE,
                       snp.exclude=snp.exclude,
                       snp.id.start=snp.id.start,
                       na.logfile=logfile)
stopifnot(res)
close(genoData)

# add snp-chunk mapping
snp.segment.map <- read.csv(config["snp_segment_mapping_file"], stringsAsFactors=FALSE, header=TRUE, as.is=TRUE)

# make them unique
#snp.segment.map$unique.segment <- 1:nrow(snp.segment.map)

segments <- snp.segment.map[snp.segment.map$chrom %in% chromosome, ]

map <- segments
#snp <- pData(snpAnnot)

#stopifnot(all(map$mb.start[2:nrow(map)] == map$mb.end[1:(nrow(map)-1)]))
levels <- c(map$bp.start, map$bp.end[nrow(map)])
fact <- cut(snpAnnot$position, breaks=levels, include.lowest=TRUE, right=FALSE)
#levels(fact) <- map$unique.segment
levels(fact) <- map$segment
snpAnnot$segment <- as.numeric(levels(fact))[fact]


# add metrics - add all of them for now
vars <- names(metrics)[!(names(metrics) %in% c("position", "rs_id"))]
for (metric.name in vars){
  snpAnnot[[metric.name]] <- metrics[[metric.name]][!snp.exclude.logical]
}

snp.filename <- file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, "_snpAnnot.RData", sep=""))
# save the imputation segments
save(snpAnnot, file=snp.filename)

# open the gds file once more
geno.gds <- GdsGenotypeReader(gds.filename)
genoData <- GenotypeData(geno.gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
close(genoData)


# remove the tmp snp annotation file
file.remove(snp.filename.tmp)

