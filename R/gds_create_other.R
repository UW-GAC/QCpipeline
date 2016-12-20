library(GWASTools)
library(OLGApipeline)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("orig_annot_scan_file",
              "orig_annot_snp_file",
              "orig_gds_file",
              "snp_segment_mapping_file",
              "out_gds_prefix",
              "out_gds_dir")
optional <- c("quality_name", "quality_minimum",
              "genotype_dim",
              "orig_snp_annot_filterCol",
              "scan_annot_familyCol", "scan_annot_individCol")
default <- c("info", 0,
             "scan,snp",
             "composite.filter",
             "family", "subjectID")
config <- setConfigDefaults(config, required, optional, default)
print(config)

# check for chromosome
if (length(args) < 2) stop("missing chromosome")
chromosome <- args[2]

# make sure it's a valid chromosome
stopifnot(chromosome %in% c("other", "failed"))

scanAnnot <- getobj(config["orig_annot_scan_file"])

# for "misc", get chromsomes 24-28, !failed
# for "misc", get all failed=TRUE
snpAnnot <- getobj(config["orig_annot_snp_file"])

if (chromosome %in% "other"){  
  snp.include <- getobj(config["snp_include_other"])
} else if (chromosome %in% "failed"){
  snp.include <- getobj(config["snp_include_failed"])
}
print(length(snp.include))

# # check for test
# if (length(args) > 3) {
#   if ((args[3]) == "test") {
#     message("testing first 5 scans")
#     scanAnnot <- scanAnnot[1:5,] # testing
#   }
# }
# nsamp <- nrow(scanAnnot)

# subset for samples included in imputation
# scanAnnot.tmp <- scanAnnot[scanAnnot$subj.plink & scanAnnot$geno.cntl == 0, ]
# construct the scan.df data frame
scan.df <- data.frame(scanID=scanAnnot$scanID[scanAnnot$subj.plink & scanAnnot$geno.cntl == 0])

# determine starting snpID - order 1:23, other, failed
olgaData <- OlgaGenotypeData(directory=config["out_gds_dir"], base=config["out_gds_prefix"])
if (chromosome %in% "other") {
  snpAnnot.tmp <- getSnpAnnotation(olgaData, 23)
  snp.id.start <- snpAnnot.tmp$snpID[nrow(snpAnnot.tmp)]
} else if (chromosome %in% "failed") {
  snpAnnot.tmp <- getSnpAnnotation(olgaData, "other")
  snp.id.start <- snpAnnot.tmp$snpID[nrow(snpAnnot.tmp)]
}
snpID <- (1:length(snp.include)) + snp.id.start
message(paste("Starting at snpID", snpID[1]))


# genoData object for original data
gds.orig <- GdsGenotypeReader(config["orig_gds_file"])
genoData.orig <- GenotypeData(gds.orig) # do not attach the scan annotation - scan annotation might not match. need to use getScanID to match them later.


# check genotypedim
stopifnot(config["genotype_dim"] %in% c("scan,snp", "snp,scan"))

# set up gds file
chrom <- getChromosome(gds.orig)[snp.include]
pos <- getPosition(gds.orig)[snp.include]

new.file <- file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, ".gds", sep=""))
gfile <- createfn.gds(new.file)

add.gdsn(gfile, "sample.id", scan.df$scanID, compress="LZMA_RA", closezip=TRUE)
add.gdsn(gfile, "snp.id", snpID, compress="LZMA_RA", closezip=TRUE)
if (hasVariable(gds.orig, "snp.allele")) {
  add.gdsn(gfile, "snp.allele", getVariable(gds.orig, "snp.allele")[snp.include],
           compress="LZMA_RA", closezip=TRUE)
}
add.gdsn(gfile, "snp.position", pos, compress="LZMA_RA", closezip=TRUE)
add.gdsn(gfile, "snp.chromosome", chrom, storage="uint8",
         compress="LZMA_RA", closezip=TRUE)
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.start", min(autosomeCode(gds.orig)))
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.end", max(autosomeCode(gds.orig)))
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "X", XchromCode(gds.orig))
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "XY", XYchromCode(gds.orig))
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "Y", YchromCode(gds.orig))
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "M", MchromCode(gds.orig))
put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "MT", MchromCode(gds.orig))     
sync.gds(gfile)

if (config["genotype_dim"] == "snp,scan"){
  valdim <- c(length(snpID), nrow(scan.df))
  geno.ord <- "snp.order"
} else if (config["genotype_dim"] == "scan,snp"){
  valdim <- c(nrow(scan.df), length(snpID))
  geno.ord <- "scan.ord"
}
gGeno <- add.gdsn(gfile, "genotype", storage="bit2", valdim=valdim)
put.attr.gdsn(gGeno, geno.ord)


# match scanAnnot (gds file) to scan.df (new gds file)
scan.df$ind.old <- match(scan.df$scanID, getScanID(gds.orig))

# loop over samples
for (ind.new in 1:nrow(scan.df)) {
  if (ind.new %% 10 == 0) message(paste("sample", ind.new, "of", nrow(scan.df)))
  stopifnot(getScanID(gds.orig)[scan.df$ind.old[ind.new]] == scan.df$scanID[ind.new])
  geno <- getGenotype(gds.orig, snp=c(1,-1), scan=c(scan.df$ind.old[ind.new], 1))
  # subset by requested genotypes
  geno <- geno[snp.include]
  # set missing values
  geno[!(geno %in% c(0,1,2))] <- 3
  
  # choose start and count based on genotype dimension requested
  if (config["genotype_dim"] == "snp,scan") {
    start <- c(1, ind.new)
    count <- c(-1, 1)
  } else if (config["genotype_dim"] == "scan,snp") {
    start <- c(ind.new, 1)
    count <- c(1, -1)    
  }
  write.gdsn(gGeno, geno, start=start, count=count)
}


closefn.gds(gfile)

# make a new snp annotation
close(gds.orig)


## write code to check against original netcdf check here.


scanAnnot.new <- ScanAnnotationDataFrame(scan.df)
save(scanAnnot.new, file=file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, "_scanAnnot.tmp.RData", sep="")))


# snp annotation
## do we want to add segments? yes, probably.
snp.df <- data.frame(snpID=snpID, snpID.old=snpAnnot$snpID[snp.include], chromosome=snpAnnot$chromosome[snp.include],
                     position=snpAnnot$position[snp.include])
snpAnnot.new <- SnpAnnotationDataFrame(snp.df)
save(snpAnnot.new, file=file.path(config["out_gds_dir"], paste(config["out_gds_prefix"], "_chr-", chromosome, "_snpAnnot.RData", sep="")))



# check - open the file
genoData <- getGenoData(olgaData, "other", snpAnnot=TRUE, scanAnnot=TRUE)
close(genoData)


