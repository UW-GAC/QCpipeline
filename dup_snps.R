##########
# Duplicate snp discordance
# Usage: R --args config.file < dup_snps.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)

# are there any scans to exclude?
if (!is.na(config["dupsnp_scan_exclude_file"])) {
  scan.exclude <- getobj(config["dupsnp_scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
length(scan.exclude)
scan.sel <- setdiff(scanID, scan.exclude)

# select snps
(snpAnnot <- getobj(config["annot_snp_file"]))
snpID <- getSnpID(snpAnnot)
snp.sel <- which(!is.na(getVariable(snpAnnot, config["annot_snp_dupSnpCol"])))
length(snp.sel)

# match up SNPs into pairs
snpset <- pData(snpAnnot)[snp.sel, c("snpID", config["annot_snp_rsIDCol"], config["annot_snp_dupSnpCol"])]
names(snpset) <- c("snpID", "rsID", "dupID")
table(table(snpset$dupID))
snpset <- snpset[order(snpset$dupID),]

dupsnp1 <- snpset[seq(1, length(snp.sel), by=2),]
names(dupsnp1)[1:3] <- c("snpID.1", "rsID.1","dupID.1")
dupsnp2 <- snpset[seq(2, length(snp.sel), by=2),]
names(dupsnp2)[1:3] <- c("snpID.2", "rsID.2","dupID.2")
snpset2 <- cbind(dupsnp1, dupsnp2)
stopifnot(allequal(snpset2[,"dupID.1"], snpset2[,"dupID.2"]))
snpset2$dupID <- snpset2$dupID.1
snpset2$dupID.1 <- NULL
snpset2$dupID.2 <- NULL

ncfile <- config["nc_geno_file"]
nc <- NcdfGenotypeReader(ncfile)
genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)


# for each pair of SNPs, find discordance
n <- nrow(snpset2)
disc.count <- rep(NA, n)
disc.sampsize <- rep(NA, n)
for(i in 1:n){
	indx1 <- which(is.element(snpID, snpset2$snpID.1[i]))
	indx2 <- which(is.element(snpID, snpset2$snpID.2[i]))
	genox1 <- getGenotype(genoData, snp=c(indx1,1), scan=c(1,-1))[scan.sel]
	genox2 <- getGenotype(genoData, snp=c(indx2,1), scan=c(1,-1))[scan.sel]
	x <- !is.na(genox1) & !is.na(genox2)
	disc.count[i] <- sum(genox1[x]!=genox2[x])
	disc.sampsize[i] <- length(genox1[x])
	if(i %%100 ==0) print(i)
}
snpset2$disc.count <- disc.count
snpset2$disc.sampsize <- disc.sampsize
snpset2$disc.fraction <- snpset2$disc.count/snpset2$disc.sampsize

# probability of discordance for various error rates
N <- max(snpset2$disc.sampsize)
prob.disc <- duplicateDiscordanceProbability(N)

# find out how  many snps fall into each category of discordance
num <- rep(NA, 8)
discordant <- snpset2$disc.count
for(i in 1:8) num[i] <- length(discordant[!is.na(discordant) & discordant>(i-1)])
prob.tbl <- cbind(prob.disc, num)

disc <- list("disc"=snpset2, "probability"=prob.tbl)
save(disc, file=config["out_dupsnp_file"])
