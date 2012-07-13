##########
# Duplicate discordance with external data
# Usage: R --args config.file < dup_disc_ext.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("out_comb_annot_scan_file", "out_comb_annot_snp_file",
              "out_comb_gds_geno_file")
optional <- c("out_disc_file", "out_disc_plot")
default <- c("dup_disc_ext.RData", "dup_disc_ext.pdf")
config <- setConfigDefaults(config, required, optional, default)
print(config)

scanAnnot <- getobj(config["out_comb_annot_scan_file"]); dim(scanAnnot)
snpAnnot <- getobj(config["out_comb_annot_snp_file"]); dim(snpAnnot)
gds <- GdsGenotypeReader(config["out_comb_gds_geno_file"])
(genoData <- GenotypeData(gds, scanAnnot=scanAnnot, snpAnnot=snpAnnot))

discord <- duplicateDiscordance(genoData, subjName.col="subjectID",
  one.pair.per.subj = FALSE)

save(discord, file=config["out_disc_file"])

# print results
nsubj <- length(discord$discordance.by.subject)
allpairs <- unlist(discord$discordance.by.subject)
npairs <- length(allpairs)
subjsum <- paste(npairs, "duplicate pairs among", nsubj, "subjects")
print(subjsum)

ndisc <- sum(discord$discordance.by.snp$discordant > 0)
print(paste(ndisc, "SNPs with > 0 discordant calls"))

# plot
disc <- allpairs[order(allpairs)]
rank <- 1:length(disc)
pdf(config["out_disc_plot"], width=6, height=6)
plot(disc, rank, xlab="discordance rate", ylab="rank", main=subjsum)
dev.off()

close(genoData)
