##########
# Combine GDS files
# Usage: R --args config.file < combine_gds.R
##########

library(GWASTools)
library(QCpipeline)
library(SNPRelate)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

scanAnnot <- getobj(config["annot_scan_file"])
ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
dupids <- intersect(scanAnnot$scanID, ext.scanAnnot$scanID)
if (length(dupids) > 0) {
  stop("scan IDs are not unique - files cannot be combined")
}

#snpgdsCombineGeno(c(config["gds_geno_file"], config["ext_gds_geno_file"]),
#                  config["out_gds_geno_file"])

proj <- c("study", "ext")

# SNP information
snpAnnot <- getobj(config["annot_snp_file"])
snp.info1 <- getVariable(snpAnnot, c("snpID", "chromosome", "position", config["annot_snp_rsIDCol"]))
names(snp.info1)[4] <- "rsID"
ext.snpAnnot <- getobj(config["ext_annot_snp_file"])
snp.info2 <- getVariable(ext.snpAnnot, c("snpID", "chromosome", "position", config["ext_annot_snp_rsIDCol"]))
names(snp.info2)[4] <- "rsID"
snp.info.list <- list(snp.info1, snp.info2)
names(snp.info.list) <- proj

# match SNPs
snp.include <- intersect(snp.info1$rsID, snp.info2$rsID)
snp.info <- snp.info1[snp.info1$rsID %in% snp.include,]
nrow(snp.info)

# sample information
samp.info1 <- getVariable(scanAnnot, c("scanID", config["annot_scan_subjectCol"]))
names(samp.info1)[2] <- "subjectID"
samp.info2 <- getVariable(ext.scanAnnot, c("scanID", config["ext_annot_scan_subjectCol"]))
names(samp.info2)[2] <- "subjectID"
samp.info <- list(samp.info1, samp.info2)
names(samp.info) <- proj

FileList <- list(config["nc_geno_file"], config["ext_nc_geno_file"])
names(FileList) <- proj


#############################################################
# create GDS files
#
#
# prepare variables
gfile1 <- createfn.gds(config["comb_gds_geno_file"])

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

add.gdsn(gfile1, "snp.id", snp.info$rsID, compress="ZIP.max", closezip=TRUE) 

# add "position"
add.gdsn(gfile1, "snp.position", as.integer(snp.info$position), compress="ZIP.max", closezip=TRUE)

# add "chromosome"
add.gdsn(gfile1, "snp.chromosome", snp.info$chromosome, storage="uint8", compress="ZIP.max", closezip=TRUE)

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
	nc <- open.ncdf(pro.fn); print(nc)
	sampleID <- get.var.ncdf(nc, "sampleID") # all scanIDs in the netCDF for pro

	# match snps
	snp.flag <- match(snp.info$rsID, snp.info.list[[pro]]$rsID)
                         
	# sample selection
       	s1 <- samp.info[[pro]]$scanID # s1: all scanIDs in samp.info[[pro]]

	# for-loop for each sample
	#cat(pro, ":", length(samp1.flag), "samples for set1\n")
        for (i in 1:length(s1))
	{
		v <- get.var.ncdf(nc, "genotype", start=c(1, i), count=c(-1, 1))[snp.flag]
		write.gdsn(gGeno1, v, start=c(1, samp1.idx), count=c(-1, 1))
		if (samp1.idx %% 100 == 0) cat(date(), "\tset1:", samp1.idx, "\n")
		samp1.idx <- samp1.idx + 1
	}
	# close the ncdf file
	close.ncdf(nc)
}

# close GDS files
closefn.gds(gfile1)
