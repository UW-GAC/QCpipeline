##########
# Duplicate discordance across two datasets
# Usage: R --args config.file type < dup_disc_2sets.R
##########

library(GWASTools)
library(QCpipeline)
library(tools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file_1", "annot_snp_file_2", "out_prefix")
optional <- c("annot_snp_clustSepCol_2", "annot_snp_mafCol_1",
              "annot_snp_snpCol_1", "annot_snp_snpCol_2")
default <- c("cluster.sep", "MAF",
             "rsID", "rsID")
config <- setConfigDefaults(config, required, optional, default)
print(config)

if (length(args) < 2) stop("missing type")
discfile <- paste(config["out_prefix"], "_", args[2], ".RData", sep="")
outfile <- paste(file_path_sans_ext(discfile), "_summary.RData", sep="")

disc <- getobj(discfile)
snp.disc <- disc$discordance.by.snp
snp.id <- row.names(snp.disc)

snpAnnot1 <- getobj(config["annot_snp_file_1"])
snp.annot1 <- pData(snpAnnot1)[snpAnnot1[[config["annot_snp_snpCol_1"]]] %in% snp.id,
                               c(config["annot_snp_snpCol_1"], config["annot_snp_mafCol_1"])]
names(snp.annot1) <- c("snp.id", "MAF")

snpAnnot2 <- getobj(config["annot_snp_file_2"])
snp.annot2 <- pData(snpAnnot2)[snpAnnot2[[config["annot_snp_snpCol_1"]]] %in% snp.id,
                               c(config["annot_snp_snpCol_1"], config["annot_snp_clustSepCol_2"])]
names(snp.annot2) <- c("snp.id", "cluster.sep")

snp.annot <- merge(snp.annot1, snp.annot2)
row.names(snp.annot) <- snp.annot$snp.id
snp.annot <- snp.annot[snp.id, c("MAF", "cluster.sep")]

stopifnot(allequal(row.names(snp.disc), row.names(snp.annot)))
snp <- cbind(snp.disc, snp.annot)

# bin by MAF and cluster.sep
maf.bins <- c(0, 0.01, 0.5)
sep.bins <- c(0, 0.4, 0.8, 1)
res <- matrix(NA, nrow=4, ncol=2, dimnames=list(
                c("clust.sep 0-0.4", "clust.sep 0.4-0.8", "clust.sep 0.8-1", "all"),
                                    c("MAF <= 0.01", "MAF > 0.01")))
nsnp <- res
for (i in 1:(length(maf.bins)-1)) {
  maf.sel <- !is.na(snp$MAF) & maf.bins[i] < snp$MAF & snp$MAF <= maf.bins[i+1]
  if (i == 1) maf.sel <- maf.sel | (!is.na(snp$MAF) & snp$MAF == maf.bins[1])
  for (j in 1:(length(sep.bins)-1)) {
    sep.sel <- !is.na(snp$cluster.sep) & sep.bins[j] < snp$cluster.sep & snp$cluster.sep <= sep.bins[j+1]
    if (j == 1) sep.sel <- sep.sel | (!is.na(snp$cluster.sep) & snp$cluster.sep == sep.bins[1])
    res[j,i] <- 1 - mean(snp$discord.rate[maf.sel & sep.sel], na.rm=TRUE)
    nsnp[j,i] <- sum(maf.sel & sep.sel)
  }
  res["all",i] <- 1 - mean(snp$discord.rate[maf.sel], na.rm=TRUE)
  nsnp["all",i] <- sum(maf.sel)
}

overall <- 1 - mean(snp$discord.rate, na.rm=TRUE)

out <- list("overall"=overall, "binned"=res, "counts"=nsnp)
save(out, file=outfile)
