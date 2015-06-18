###################
# Plotting script for locus zoom plots (currently assumes results are stored chromosome-level from gmmat mixed models)
###################

library(GWASTools)
library(QCpipeline)
library(GWASbyChr)
library(OLGAtools)

args <- commandArgs(trailingOnly=TRUE)
print(args)

if (length(args) < 1) stop("missing config file")
(configFile <- args[1])

if (length(args) < 2) stop("missing LD population")
(pop <- args[2])

if (length(args) < 3) stop("missing snpID")
(snpID <- args[3])

flankingRegion <- 500

config <- readConfig(configFile)
required <- c("geno_file", "out_assoc_prefix", "scan_include_file", "out_plot_prefix")
optional <- c()
default <- c()
config <- setConfigDefaults(config, required, optional, default)

annotByChr <- AnnotationByChr(config["geno_file"], suffix="_snpAnnot.RData", chromSep="_chr-")

assocByChr <- AssocResultsByChr(dirname(config["out_assoc_prefix"]),
	base=basename(config["out_assoc_prefix"]),
	chromSep="_chr-")

assocAnnot <- locusZoomPrepareResults(annotByChr, assocByChr, snpID, flankingRegion, macThreshold=50)
# for the GMMAT imputation, rsID is whatever was in impute2 output -- should probably change this in the code that makes the SNP annotations files.
# since there are no rsIDs here, just use the locus zoom name
# you may end up with some duplicated lz.names here, so we should modify locusZoomPrepareResults to optionally use rsID (specify a column with an argument) or just use chr:pos inidcators
assocAnnot$lz.name <- assocAnnot$lz.pos

refsnp <- assocAnnot[assocAnnot$refsnp, ]

if (pop %in% c("analysis")){

	# gentoype data
	genoByChr <- GenotypeDataByChr(config["geno_file"])
	genoData <- getGenoData(genoByChr, chromosome=chromToInt(refsnp$chromosome), snpAnnot=TRUE)

	# included scans
	scan.include <- getobj(config["scan_include_file"])

	ldFile <- tempfile()
	ld.cmd <- locusZoomLD(genoData, assocAnnot, scan.include=scan.include, filename=ldFile)
	ld.title <- "analysis"

} else {
	ld.cmd <- paste("--pop", pop, "--source 1000G_March2012")
	ld.title <- paste("1000G", pop)
}

pval <- "Score.pval"

print(paste0("plotting ", pval, "..."))

## write association test results
assoc <- assocAnnot[, c("lz.name", pval, "lz.obs")]
names(assoc) <- c("MarkerName", "P-value", "annotation")
assoc.filename <- tempfile()
write.table(assoc, file=assoc.filename, row.names=F, quote=F, sep="\t")


prefix <- paste0(basename(config["out_plot_prefix"]), "_", "locuszoom_snp", snpID)

title <- paste(refsnp$lz.name, "- LD:", ld.title, "- MAF:", formatC(refsnp$MAF, digits=3))

setwd(dirname(config["out_plot_prefix"]))

command <- paste("/projects/geneva/geneva_sata/apps/locuszoom/bin/locuszoom",
               "theme=publication",
               "--cache None",
               "--no-date",
               "--chr", refsnp$chromosome,
               "--plotonly",
               "--metal", assoc.filename,
               ld.cmd,
               paste0("--refsnp \"", refsnp$lz.name, "\""),
               "--build hg19",
               paste0("--flank ", flankingRegion, "kb"),
               "--prefix ", prefix,
               "annotCol='annotation' annotOrder='refsnp_genotyped,refsnp_imputed,genotyped,imputed' annotPch='23,25,21,4'",
               paste0("title=\"", title, "\""),
               paste0("signifLine=\"", -log10(5e-8), "\" signifLineColor=\"gray\" signifLineWidth=\"2\""))

cat(paste(command, "\n"))

system(command)