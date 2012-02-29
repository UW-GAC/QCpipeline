##########
# Plots to check genotyping bach
# Usage: R --args config.file test.type < batch_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])
print(config)

# check for test
if (length(args) < 2) stop("missing test type (chisq or fisher)")
type <- args[2]
  
(scanAnnot <- getobj(config["annot_scan_file"]))

if (type == "chisq") {
  batch.res <- getobj(paste(config["out_chisq_file"], "RData", sep="."))
} else if (type == "fisher") {
  batch.res <- getobj(paste(config["out_fisher_file"], "RData", sep="."))
} else {
  stop("test type must be chisq or fisher")
}

batches <- names(batch.res$lambda)
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

pch <- rep(1, length(batches))
redo <- getVariable(scanAnnot, config["annot_scan_redoCol"])
if (!is.null(redo)) {
  redobatches <- unique(batch[redo %in% "Y"])
  pch[batches %in% redobatches] <- 2
}


if (type == "chisq") {
  pdf(config["out_meanchisq_race_plot"], width=6, height=6)
  plot(racefrac, batch.res$mean.chisq, ylab=expression(paste("mean ", chi^2, " test statistic")), xlab=paste("fraction of", majority, "samples per batch"), pch=pch)
  abline(v=mean(racefrac), lty=2) # mean over all plates
  if (!is.null(redo)) {
    legend(bestLegendPos(racefrac, batch.res$mean.chisq), c("redo","mean"), pch=c(2,-1), lty=c(0,2))
  } else {
    legend(bestLegendPos(racefrac, batch.res$mean.chisq), c("mean"), pch=c(-1), lty=c(2))
  }
  dev.off()
} else if (type == "fisher") {
  pdf(config["out_meanor_race_plot"], width=6, height=6)
  plot(racefrac, batch.res$mean.or, ylab="mean Fisher's OR", xlab=paste("fraction of", majority, "samples per batch"), pch=pch)
  abline(v=mean(racefrac), lty=2) # mean over all plates
  if (!is.null(redo)) {
    legend(bestLegendPos(racefrac, batch.res$mean.or), c("redo","mean"), pch=c(2,-1), lty=c(0,2))
  } else {
    legend(bestLegendPos(racefrac, batch.res$mean.or), c("mean"), pch=c(-1), lty=c(2))
  }
  dev.off()
}

pdf(config["out_lambda_race_plot"], width=6, height=6)
plot(racefrac, batch.res$lambda, ylab=expression(paste("genomic inflation factor ", lambda)), xlab=paste("fraction of", majority, "samples per batch"), pch=pch)
abline(v=mean(racefrac), lty=2) # mean over all plates
if (!is.null(redo)) {
  legend(bestLegendPos(racefrac, batch.res$lambda), c("redo","mean"), pch=c(2,-1), lty=c(0,2))
} else {
  legend(bestLegendPos(racefrac, batch.res$lambda), c("mean"), pch=c(-1), lty=c(2))
}
dev.off()


# distribution of number of samples per batch
pdf(config["out_hist_plot"], width=6, height=6)
hist(table(batch), xlab="number of samples per batch", ylab="number of batches", main="")
dev.off()


# mean autosomal missing call rate per batch
missing <- getVariable(scanAnnot, config["annot_scan_missAutoCol"])
bmiss <- rep(NA,n); names(bmiss) <- batches
bn <- rep(NA,n); names(bn) <- batches
for(i in 1:n) {
  x <- missing[is.element(batch, batches[i])]
  bmiss[i] <- mean(x)
  bn[i] <- length(x)
}
pdf(config["out_meanmcr_nscan_plot"], width=6, height=6)
plot(bn, bmiss, xlab="number of samples per batch", ylab="mean autosomal missing call rate", pch=pch)
y <- lm(bmiss ~ bn)
#abline(y$coefficients)
anova(y)
if (!is.null(redo)) legend(bestLegendPos(bn, bmiss), "redo", pch=2)
dev.off()

if (type == "chisq") {
  pdf(config["out_meanmcr_meanchisq_plot"], width=6, height=6)
  tmp <- batch.res$mean.chisq[match(names(bmiss), names(batch.res$mean.chisq))]
  plot(tmp, bmiss, xlab=expression(paste("mean ", chi^2, " test statistic")), ylab="mean autosomal missing call rate", pch=pch)
  if (!is.null(redo)) legend(bestLegendPos(tmp, bmiss), "redo", pch=2)
  dev.off()
} else if (type == "fisher") {
  pdf(config["out_meanmcr_meanor_plot"], width=6, height=6)
  tmp <- batch.res$mean.or[match(names(bmiss), names(batch.res$mean.or))]
  plot(tmp, bmiss, xlab="mean Fisher's OR", ylab="mean autosomal missing call rate", pch=pch)
  if (!is.null(redo)) legend(bestLegendPos(tmp, bmiss), "redo", pch=2)
  dev.off()
}

if (type == "chisq") {
  pdf(config["out_meanchisq_ncsan_plot"], width=6, height=6)
  tmp <- batch.res$mean.chisq[match(names(bn), names(batch.res$mean.chisq))]
  plot(bn, tmp, ylab=expression(paste("mean ", chi^2, " test statistic")), xlab="number of samples per batch", pch=pch)
  if (!is.null(redo)) legend(bestLegendPos(bn, tmp), "redo", pch=2)
  dev.off()
} else if (type == "fisher") {
  pdf(config["out_meanor_ncsan_plot"], width=6, height=6)
  tmp <- batch.res$mean.or[match(names(bn), names(batch.res$mean.or))]
  plot(bn, tmp, ylab="mean Fisher's OR", xlab="number of samples per batch", pch=pch)
  if (!is.null(redo)) legend(bestLegendPos(bn, tmp), "redo", pch=2)
  dev.off()
}


# batch effect on missing call rate
# make plate names shorter for plotting
batch <- getVariable(scanAnnot, config["annot_scan_batchCol"])
batchLabel <- vapply(strsplit(as.character(batch), "-"), function(x) x[[1]][1], "a")
missing <- getVariable(scanAnnot, config["annot_scan_missAutoCol"])
model <- log10(missing) ~ as.factor(batchLabel)
lm.result <- lm(model)
anova(lm.result)
pdf(config["out_mcr_plot"], width=6, height=6)
par(mar=c(6, 4, 4, 2) + 0.1)
boxplot(model, varwidth=TRUE, las=2, ylab="log10(autosomal missing call rate)", main=config["annot_scan_batchCol"])
dev.off()


# chromosome 1 intensity by batch
mninten <- getobj(config["inten_file"])
mninten <- mninten[[1]]
pdf(config["out_inten_plot"], width=6, height=6)
par(mar=c(6, 4, 4, 2) + 0.1)
boxplot(mninten[,"1"] ~ as.factor(batchLabel), varwidth=TRUE, las=2, ylab="mean chromosome 1 intensity", main=config["annot_scan_batchCol"])
dev.off()
