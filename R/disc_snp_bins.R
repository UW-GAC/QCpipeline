##########
# Duplicate discordance across two datasets
# Usage: R --args config.file type < disc_snp_bins.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_snp_file_1", "annot_snp_file_2", "out_prefix")
optional <- c("annot_snp_clustSepCol_2", "annot_snp_mafCol_1",
              "annot_snp_snpCol_1", "annot_snp_snpCol_2",
              "bins_maf", "bins_clustSep",
              "out_summary_prefix", "summary_include_file")
default <- c("cluster.sep", "MAF",
             "rsID", "rsID",
             "0 0.01 0.05 0.5", "0 0.4 0.8 1",
             NA, NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

if (length(args) < 2) stop("missing type")
type <- args[2]

discfile <- paste(config["out_prefix"], "_", type, ".RData", sep="")
disc <- getobj(discfile)

if (is.na(config["out_summary_prefix"])) {
  outfile <- paste(config["out_prefix"], "_", type, "_summary.RData", sep="")
} else {
  outfile <- paste(config["out_summary_prefix"], "_", type, "_summary.RData", sep="")
}

if (type == "senspec") {
  snp.disc <- disc
  metric <- names(disc)[-1]
} else {
  snp.disc <- disc$discordance.by.snp
  # metric to summarize = concordance rate
  snp.disc$conc.rate <- 1 - snp.disc$discord.rate
  metric <- "conc.rate"
}
snp.id <- row.names(snp.disc)

snpAnnot1 <- getobj(config["annot_snp_file_1"])
snp.annot1 <- pData(snpAnnot1)[snpAnnot1[[config["annot_snp_snpCol_1"]]] %in% snp.id,
                               c(config["annot_snp_snpCol_1"], config["annot_snp_mafCol_1"])]
names(snp.annot1) <- c("snp.id", "MAF")

snpAnnot2 <- getobj(config["annot_snp_file_2"])
snp.annot2 <- pData(snpAnnot2)[snpAnnot2[[config["annot_snp_snpCol_2"]]] %in% snp.id,
                               c(config["annot_snp_snpCol_2"], config["annot_snp_clustSepCol_2"])]
names(snp.annot2) <- c("snp.id", "cluster.sep")

snp.annot <- merge(snp.annot1, snp.annot2)
row.names(snp.annot) <- snp.annot$snp.id
snp.annot <- snp.annot[snp.id, c("MAF", "cluster.sep")]

stopifnot(allequal(row.names(snp.disc), row.names(snp.annot)))
snp <- cbind(snp.disc, snp.annot)

if (!is.na(config["summary_include_file"])) {
  snp.include <- getobj(config["summary_include_file"])
  snp <- snp[snp.include,]
}

# bin by MAF and cluster.sep
maf.bins <- as.numeric(unlist(strsplit(config["bins_maf"], " ", fixed=TRUE), use.names=FALSE))
sep.bins <- as.numeric(unlist(strsplit(config["bins_clustSep"], " ", fixed=TRUE), use.names=FALSE))

# loop over metrics in file
out <- list()
for (m in metric) {
  res <- matrix(NA, nrow=length(sep.bins), ncol=length(maf.bins), dimnames=list(
                  c(paste("clust.sep", sep.bins[1:(length(sep.bins)-1)], "-",
                          sep.bins[2:length(sep.bins)]), "all"),
                  c(paste("MAF", maf.bins[1:(length(maf.bins)-1)], "-",
                          maf.bins[2:length(maf.bins)]), "all")))
  nsnp <- res
  for (i in 1:(length(maf.bins)-1)) {
    maf.sel <- !is.na(snp$MAF) & maf.bins[i] < snp$MAF & snp$MAF <= maf.bins[i+1]
    if (i == 1) maf.sel <- maf.sel | (!is.na(snp$MAF) & snp$MAF == maf.bins[1])
    for (j in 1:(length(sep.bins)-1)) {
      sep.sel <- !is.na(snp$cluster.sep) & sep.bins[j] < snp$cluster.sep & snp$cluster.sep <= sep.bins[j+1]
      if (j == 1) sep.sel <- sep.sel | (!is.na(snp$cluster.sep) & snp$cluster.sep == sep.bins[1])
      res[j,i] <- mean(snp[maf.sel & sep.sel, m], na.rm=TRUE)
      nsnp[j,i] <- sum(maf.sel & sep.sel & snp$npair > 0)
      res[j,"all"] <- mean(snp[sep.sel, m], na.rm=TRUE)
      nsnp[j,"all"] <- sum(sep.sel & snp$npair > 0)
    }
    res["all",i] <- mean(snp[maf.sel, m], na.rm=TRUE)
    nsnp["all",i] <- sum(maf.sel & snp$npair > 0)
  }

  res["all","all"] <- mean(snp[,m], na.rm=TRUE)
  nsnp["all","all"] <- sum(snp$npair > 0)

  out[[m]] <- list("metric"=res, "counts"=nsnp)
}
if (length(out) == 1) out <- out[[1]]
save(out, file=outfile)
