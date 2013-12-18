##########
# Combine two netCDF files into a single GDS file
# Usage: R --args config.file < combine_gds.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file", "ext_annot_scan_file",
              "ext_annot_snp_file", "ext_nc_geno_file", "gds_geno_file")
optional <- c("annot_scan_subjectCol", "annot_snp_rsIDCol",
              "comb_scan_exclude_file", "comb_snp_exclude_file",
              "ext_annot_scan_subjectCol", "ext_annot_snp_rsIDCol",
              "out_comb_annot_scan_file", "out_comb_annot_snp_file", "out_comb_gds_geno_file")
default <- c("subjectID", "rsID", NA, NA, "subjectID", "rsID",
             "comb_scan_annot.RData", "comb_snp_annot.RData", "comb_geno.gds")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["annot_scan_file"]); dim(scanAnnot)
# this might be a subject-level netCDF file, so subset annotation
nc <- GenotypeReader(config["gds_geno_file"])
scanAnnot <- scanAnnot[match(getScanID(nc), getScanID(scanAnnot)),]; dim(scanAnnot)
close(nc)

ext.scanAnnot <- getobj(config["ext_annot_scan_file"]); dim(ext.scanAnnot)
# this might be a subject-level netCDF file, so subset annotation
nc <- GenotypeReader(config["ext_nc_geno_file"])
ext.scanAnnot <- ext.scanAnnot[match(getScanID(nc), getScanID(ext.scanAnnot)),]; dim(ext.scanAnnot)
close(nc)
dupids <- intersect(scanAnnot$scanID, ext.scanAnnot$scanID)
if (length(dupids) > 0) {
  stop("scan IDs are not unique - files cannot be combined")
}

# any scans to exclude?
if (!is.na(config["comb_scan_exclude_file"])) {
  scan.exclude <- getobj(config["comb_scan_exclude_file"])
  scanAnnot <- scanAnnot[!(scanAnnot$scanID %in% scan.exclude),]
  ext.scanAnnot <- ext.scanAnnot[!(ext.scanAnnot$scanID %in% scan.exclude),]
}
dim(scanAnnot)
dim(ext.scanAnnot)

proj <- c("study", "ext")

# SNP information
snpAnnot <- getobj(config["annot_snp_file"])
snp.info1 <- getVariable(snpAnnot, c("snpID", "chromosome", "position", config["annot_snp_rsIDCol"],
                                     "alleleA", "alleleB"))
names(snp.info1)[4] <- "rsID"
ext.snpAnnot <- getobj(config["ext_annot_snp_file"])
snp.info2 <- getVariable(ext.snpAnnot, c("snpID", "chromosome", "position", config["ext_annot_snp_rsIDCol"]))
names(snp.info2)[4] <- "rsID"
snp.info.list <- list(snp.info1, snp.info2)
names(snp.info.list) <- proj

# match SNPs
snp.include <- merge(snp.info1[,c("rsID", "chromosome", "position")],
                     snp.info2[,c("rsID", "chromosome", "position")])
snp.info <- snp.info1[snp.info1$rsID %in% snp.include$rsID,]
nrow(snp.info)

# any snps to exclude?
if (!is.na(config["comb_snp_exclude_file"])) {
  snp.exclude <- getobj(config["comb_snp_exclude_file"])
  snp.info <- snp.info[!(snp.info$snpID %in% snp.exclude),]
}
nrow(snp.info)

# sample information
samp.info1 <- getVariable(scanAnnot, c("scanID", config["annot_scan_subjectCol"], "sex"))
names(samp.info1)[2] <- "subjectID"
samp.info2 <- getVariable(ext.scanAnnot, c("scanID", config["ext_annot_scan_subjectCol"], "sex"))
names(samp.info2)[2] <- "subjectID"
samp.info <- list(samp.info1, samp.info2)
names(samp.info) <- proj

FileList <- list(config["gds_geno_file"], config["ext_nc_geno_file"])
names(FileList) <- proj

# save annotation files
samp.info.all <- ScanAnnotationDataFrame(rbind(samp.info1, samp.info2))
save(samp.info.all, file=config["out_comb_annot_scan_file"])
snp.info.all <- SnpAnnotationDataFrame(snp.info[,c("snpID", "rsID", "chromosome", "position", "alleleA", "alleleB")])
save(snp.info.all, file=config["out_comb_annot_snp_file"])


#############################################################
# create GDS files
#
#
# prepare variables
gfile1 <- createfn.gds(config["out_comb_gds_geno_file"])

# add "sample.id"
sample1.id <- NULL
for (pro in proj)
{
	s1 <- with(samp.info[[pro]], scanID)
	sample1.id <- c(sample1.id, s1)
}

length(sample1.id)

add.gdsn(gfile1, "sample.id", sample1.id, compress="ZIP.max", closezip=TRUE) # scanID

# add "snp.id"
nSNP <- dim(snp.info)[1] 

add.gdsn(gfile1, "snp.id", snp.info$snpID, compress="ZIP.max", closezip=TRUE)

# add "rs.id"
add.gdsn(gfile1, "snp.rs.id", snp.info$rsID, compress="ZIP.max", closezip=TRUE) 

# add "position"
add.gdsn(gfile1, "snp.position", as.integer(snp.info$position), compress="ZIP.max", closezip=TRUE)

# add "chromosome"
add.gdsn(gfile1, "snp.chromosome", snp.info$chromosome, storage="uint8", compress="ZIP.max", closezip=TRUE)

# add "allele"
add.gdsn(gfile1, "snp.allele", paste(snp.info$alleleA, snp.info$alleleB, sep="/"), compress="ZIP.max", closezip=TRUE)

# sync GDS files
sync.gds(gfile1)

# add "genotype", 2 bits to store one genotype - 0, 1, 2
gGeno1 <- add.gdsn(gfile1, "genotype", storage="bit2", valdim=c(nSNP, length(sample1.id)))
put.attr.gdsn(gGeno1, "snp.order")

# for-loop each project
# samples in GDS are ordered by project, then by samp.info file

# sample index
samp1.idx <- 1

for (pro in names(FileList))
{
	pro.fn <- FileList[[pro]]

	# start
	cat("\n")
	cat(date(), "\tProject Name,", pro, ":", pro, "\n")

	# open netCDF
## 	nc <- open.ncdf(pro.fn); print(nc)
## 	sampleID <- get.var.ncdf(nc, "sampleID") # all scanIDs in the netCDF for pro
        nc <- GenotypeReader(pro.fn); print(nc)
        sampleID <- getScanID(nc)

	# match snps
	snp.flag <- match(snp.info$rsID, snp.info.list[[pro]]$rsID)
                         
	# sample selection
       	s1 <- which(sampleID %in% samp.info[[pro]]$scanID) # s1: all scanIDs in samp.info[[pro]]

	# for-loop for each sample
	#cat(pro, ":", length(samp1.flag), "samples for set1\n")
        for (i in s1)
	{
## 		v <- get.var.ncdf(nc, "genotype", start=c(1, i), count=c(-1, 1))[snp.flag]
		v <- getGenotype(nc, snp=c(1,-1), scan=c(i,1))[snp.flag]
                v[is.na(v)] <- 3  ## 3 is missing value for GDS genotype file
		write.gdsn(gGeno1, v, start=c(1, samp1.idx), count=c(-1, 1))
		if (samp1.idx %% 100 == 0) cat(date(), "\tset1:", samp1.idx, "\n")
		samp1.idx <- samp1.idx + 1
	}
	# close the ncdf file
## 	close.ncdf(nc)
	close(nc)
}

# close GDS files
closefn.gds(gfile1)

# check
## check annotation
gds <- GdsGenotypeReader(config["out_comb_gds_geno_file"])
(gData <- GenotypeData(gds, snpAnnot=snp.info.all, scanAnnot=samp.info.all))

## check genotypes
nc1 <- GenotypeReader(config["gds_geno_file"])
nc2 <- GenotypeReader(config["ext_nc_geno_file"])

scanID <- getScanID(gData)
scan1 <- getScanID(nc1)
scan2 <- getScanID(nc2)
snp.flag1 <- match(snp.info$rsID, snp.info.list[["study"]]$rsID)
snp.flag2 <- match(snp.info$rsID, snp.info.list[["ext"]]$rsID)
for (i in 1:length(scanID)) {
  geno <- getGenotype(gData, snp=c(1,-1), scan=c(i,1))
  if (scanID[i] %in% scan1) {
    geno.orig <- getGenotype(nc1, snp=c(1,-1), scan=c(which(scan1 == scanID[i]), 1))[snp.flag1]
  } else {
    geno.orig <- getGenotype(nc2, snp=c(1,-1), scan=c(which(scan2 == scanID[i]), 1))[snp.flag2]
  }
  stopifnot(allequal(geno, geno.orig))
}

close(gData)
close(nc1)
close(nc2)
