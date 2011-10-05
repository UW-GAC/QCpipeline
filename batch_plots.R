##########
# Plots to check genotyping bach
# Usage: R --args config.file < batch_plots.R
##########

library(GWASTools)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config.table <- read.table(args[1], as.is=TRUE)
config <- config.table[,2]
names(config) <- config.table[,1]
print(config.table)

(scanAnnot <- getobj(config["annot_scan_file"]))

batch.chisq <- getobj(paste(config["out_chisq_file"], "RData", sep="."))

batches <- names(batch.chisq$mean.chisq)
n <- length(batches)
batch <- getVariable(scanAnnot, config["annot_scan_batchCol"])
race <- getVariable(scanAnnot, config["annot_scan_raceCol"])
racetbl <- table(race)
majority <- names(racetbl)[racetbl == max(racetbl)]
racefrac <- rep(NA, n); names(racefrac) <- batches
for (i in 1:n) {
  nsamp <- sum(batch %in% batches[i] & !is.na(race))
  nmaj <- sum(batch %in% batches[i] & race %in% majority)
  racefrac[i] <- nmaj / nsamp
}

pdf(config["out_meanchisq_race_plot"], width=6, height=6)
plot(racefrac, batch.chisq$mean.chisq, ylab=expression(paste("mean ", chi^2, " test statistic")), xlab=paste("fraction of", majority, "samples per batch"))
abline(v=mean(racefrac), lty=2) # mean over all plates
dev.off()

pdf(config["out_lambda_race_plot"], width=6, height=6)
plot(racefrac, batch.chisq$lambda, ylab=expression(paste("genomic inflation factor ", lambda)), xlab=paste("fraction of", majority, "samples per batch"))
abline(v=mean(racefrac), lty=2) # mean over all plates
dev.off()


# distribution of number of samples per batch
pdf(config["out_hist_plot"], width=6, height=6)
hist(table(batch), xlab="number of samples per batch", ylab="number of batches", main="")
dev.off()


# mean autosomal missing call rate per batch
missing <- getVariable(scanAnnot, config["annot_scan_missAutoCol"])
batches <- unique(na.omit(batch))
n <- length(batches)
bmiss <- rep(NA,n); names(bmiss) <- batches
bn <- rep(NA,n); names(bn) <- batches
for(i in 1:n) {
  x <- missing[is.element(batch, batches[i])]
  bmiss[i] <- mean(x)
  bn[i] <- length(x)
}
pdf(config["out_meanmcr_nscan_plot"], width=6, height=6)
plot(bn, bmiss, xlab="number of samples per batch", ylab="mean autosomal missing call rate")
y <- lm(bmiss ~ bn)
abline(y$coefficients)
anova(y)
dev.off()

pdf(config["out_meanmcr_meanchisq_plot"], width=6, height=6)
tmp <- batch.chisq$mean.chisq[match(names(bmiss), names(batch.chisq$mean.chisq))]
plot(tmp, bmiss, xlab=expression(paste("mean ", chi^2, " test statistic")), ylab="mean autosomal missing call rate")
dev.off()


# study only
# are there any scans to exclude?
scanID <- getScanID(scanAnnot)
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanID))
} else {
  scan.exclude <- NULL
}
incl.sel <- !(scanID %in% scan.exclude)
study <- scanAnnot[incl.sel,]

# batch effect on missing call rate
# make plate names shorter for plotting
batch <- getVariable(study, config["annot_scan_batchCol"])
batchLabel <- vapply(strsplit(as.character(batch), "-"), function(x) x[[1]][1], "a")
missing <- getVariable(study, config["annot_scan_missAutoCol"])
model <- log10(missing) ~ as.factor(batchLabel)
lm.result <- lm(model)
anova(lm.result)
pdf(config["out_mcr_plot"], width=6, height=6)
boxplot(model, varwidth=TRUE, las=2, ylab="log10(autosomal missing call rate)", main=config["annot_scan_batchCol"])
dev.off()


# chromosome 1 intensity by batch
mninten <- getobj(config["inten_file"])
mninten <- mninten[[1]]
pdf(config["out_inten_plot"], width=6, height=6)
boxplot(mninten[incl.sel,"1"] ~ as.factor(batchLabel), varwidth=TRUE, las=2, ylab="mean chromosome 1 intensity", main=config["annot_scan_batchCol"])
dev.off()
