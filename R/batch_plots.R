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

# check config and set defaults
required <- c("annot_scan_file", "annot_scan_raceCol", "inten_file")
optional <- c("annot_scan_batchCol", "annot_scan_missAutoCol", "annot_scan_redoCol",
              "out_chisq_file", "out_fisher_file", "out_hist_plot", "out_inten_plot",
              "out_lambda_race_plot", "out_mcr_plot", "out_meanchisq_nscan_plot",
              "out_meanchisq_race_plot", "out_meanmcr_meanchisq_plot",
              "out_meanmcr_meanor_plot", "out_meanmcr_nscan_plot",
              "out_meanor_nscan_plot", "out_meanor_race_plot")
default <- c("Sample.Plate", "miss.e1.auto", "Redo.Processing.Plate", "batch_chisq",
             "batch_fisher", "batch_nscan_hist.pdf", "batch_chr1inten.pdf",
             "batch_lambda_race.pdf", "batch_mcr.pdf", "batch_meanchisq_nscan.pdf",
             "batch_meanchisq_race.pdf", "batch_meanmcr_meanchisq.pdf",
             "batch_meanmcr_meanor.pdf", "batch_meanmcr_nscan.pdf",
             "batch_meanor_nscan.pdf", "batch_meanor_race.pdf")
config <- setConfigDefaults(config, required, optional, default)
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

pcol <- rep("black", length(batches))
bpcol <- rep("black", length(batches))
redo <- getVariable(scanAnnot, config["annot_scan_redoCol"])
if (!is.null(redo)) {
  redobatches <- unique(batch[redo %in% c("Y", "Yes", "yes", "YES", TRUE)])
  pcol[batches %in% redobatches] <- "red"
  # boxplot labels will always be sorted
  bpcol[sort(batches) %in% redobatches] <- "red"
}


if (type == "chisq") {
  pdf(config["out_meanchisq_race_plot"], width=6, height=6)
  plot(racefrac, batch.res$mean.chisq, ylab=expression(paste("mean ", chi^2, " test statistic")), xlab=paste("fraction of", majority, "samples per batch"), col=pcol)
  abline(v=mean(racefrac), lty=2) # mean over all plates
  if (!is.null(redo)) {
    legend(bestLegendPos(racefrac, batch.res$mean.chisq), c("redo","mean"), col=c("red","black"), pch=c(1,-1), lty=c(0,2))
  } else {
    legend(bestLegendPos(racefrac, batch.res$mean.chisq), c("mean"), pch=c(-1), lty=c(2))
  }
  dev.off()
} else if (type == "fisher") {
  pdf(config["out_meanor_race_plot"], width=6, height=6)
  plot(racefrac, batch.res$mean.or, ylab="mean Fisher's OR", xlab=paste("fraction of", majority, "samples per batch"), col=pcol)
  abline(v=mean(racefrac), lty=2) # mean over all plates
  if (!is.null(redo)) {
    legend(bestLegendPos(racefrac, batch.res$mean.or), c("redo","mean"),  col=c("red","black"), pch=c(1,-1), lty=c(0,2))
  } else {
    legend(bestLegendPos(racefrac, batch.res$mean.or), c("mean"), pch=c(-1), lty=c(2))
  }
  dev.off()
}

pdf(config["out_lambda_race_plot"], width=6, height=6)
plot(racefrac, batch.res$lambda, ylab=expression(paste("genomic inflation factor ", lambda)), xlab=paste("fraction of", majority, "samples per batch"), col=pcol)
abline(v=mean(racefrac), lty=2) # mean over all plates
if (!is.null(redo)) {
  legend(bestLegendPos(racefrac, batch.res$lambda), c("redo","mean"),  col=c("red","black"), pch=c(1,-1), lty=c(0,2))
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
plot(bn, bmiss, xlab="number of samples per batch", ylab="mean autosomal missing call rate", col=pcol)
y <- lm(bmiss ~ bn)
#abline(y$coefficients)
anova(y)
if (!is.null(redo)) legend(bestLegendPos(bn, bmiss), "redo", col="red", pch=1)
dev.off()

if (type == "chisq") {
  pdf(config["out_meanmcr_meanchisq_plot"], width=6, height=6)
  tmp <- batch.res$mean.chisq[match(names(bmiss), names(batch.res$mean.chisq))]
  plot(tmp, bmiss, xlab=expression(paste("mean ", chi^2, " test statistic")), ylab="mean autosomal missing call rate", col=pcol)
  if (!is.null(redo)) legend(bestLegendPos(tmp, bmiss), "redo", col="red", pch=1)
  dev.off()
} else if (type == "fisher") {
  pdf(config["out_meanmcr_meanor_plot"], width=6, height=6)
  tmp <- batch.res$mean.or[match(names(bmiss), names(batch.res$mean.or))]
  plot(tmp, bmiss, xlab="mean Fisher's OR", ylab="mean autosomal missing call rate", col=pcol)
  if (!is.null(redo)) legend(bestLegendPos(tmp, bmiss), "redo", col="red", pch=1)
  dev.off()
}

if (type == "chisq") {
  pdf(config["out_meanchisq_nscan_plot"], width=6, height=6)
  tmp <- batch.res$mean.chisq[match(names(bn), names(batch.res$mean.chisq))]
  plot(bn, tmp, ylab=expression(paste("mean ", chi^2, " test statistic")), xlab="number of samples per batch", col=pcol)
  if (!is.null(redo)) legend(bestLegendPos(bn, tmp), "redo", col="red", pch=1)
  dev.off()
} else if (type == "fisher") {
  pdf(config["out_meanor_nscan_plot"], width=6, height=6)
  tmp <- batch.res$mean.or[match(names(bn), names(batch.res$mean.or))]
  plot(bn, tmp, ylab="mean Fisher's OR", xlab="number of samples per batch", col=pcol)
  if (!is.null(redo)) legend(bestLegendPos(bn, tmp), "redo", col="red", pch=1)
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
boxplot(model, varwidth=TRUE, las=2, ylab="log10(autosomal missing call rate)", main=config["annot_scan_batchCol"], border=bpcol)
dev.off()


# chromosome 1 intensity by batch
mninten <- getobj(config["inten_file"])
mninten <- mninten[[1]]
pdf(config["out_inten_plot"], width=6, height=6)
par(mar=c(6, 4, 4, 2) + 0.1)
boxplot(mninten[,"1"] ~ as.factor(batchLabel), varwidth=TRUE, las=2, ylab="mean chromosome 1 intensity", main=config["annot_scan_batchCol"], border=bpcol)
dev.off()
