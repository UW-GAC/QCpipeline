##########
# Plots to check genotyping bach
# Usage: R --args config.file test.type < batch_plots.R
##########

library(GWASTools)
library(QCpipeline)
library(plyr)
library(ggplot2)
sessionInfo()

theme_set(theme_bw() + theme(legend.position="top"))

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_scan_raceCol", "inten_file")
optional <- c("annot_scan_batchCol", "annot_scan_missAutoCol", "annot_scan_redoCol",
              "out_batch_prefix", "out_hist_plot", "out_inten_plot",
              "out_lambda_race_plot", "out_mcr_plot", "out_meanchisq_nscan_plot",
              "out_meanchisq_race_plot", "out_meanmcr_meanchisq_plot",
              "out_meanmcr_meanor_plot", "out_meanmcr_nscan_plot",
              "out_meanor_nscan_plot", "out_meanor_race_plot",
              "scan_exclude_file")
default <- c("Sample.Plate", "miss.e1.auto", "Redo.Processing.Plate",
             "batch_test", "batch_nscan_hist.pdf", "batch_chr1inten.pdf",
             "batch_lambda_race.pdf", "batch_mcr.pdf", "batch_meanchisq_nscan.pdf",
             "batch_meanchisq_race.pdf", "batch_meanmcr_meanchisq.pdf",
             "batch_meanmcr_meanor.pdf", "batch_meanmcr_nscan.pdf",
             "batch_meanor_nscan.pdf", "batch_meanor_race.pdf",
             NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

# check for test
if (length(args) < 2) stop("missing test type (chisq or fisher)")
type <- args[2]
  
(scanAnnot <- getobj(config["annot_scan_file"]))

# check for scans to exclude
if (!is.na(config["scan_exclude_file"])) {
  scan.exclude <- getobj(config["scan_exclude_file"])
  stopifnot(all(scan.exclude %in% scanAnnot$scanID))
} else {
	scan.exclude <- NULL
}
length(scan.exclude)

# exclude hapmaps from plots
if (!is.na(config["annot_scan_hapmapCol"])) {
  hapmap <- getVariable(scanAnnot, config["annot_scan_hapmapCol"])
  hm.ids <- getScanID(scanAnnot)[hapmap == 1]
  scan.exclude <- union(scan.exclude, hm.ids)
}
length(scan.exclude)

batch.res <- getobj(paste(config["out_batch_prefix"], "RData", sep="."))

scan.index <- !(scanAnnot$scanID %in% scan.exclude)
stopifnot(all(!(getScanID(scanAnnot, index=scan.index) %in% scan.exclude)))
sum(scan.index)

scan <- getVariable(scanAnnot, c("scanID",
                                 config["annot_scan_batchCol"],
                                 config["annot_scan_raceCol"],
                                 config["annot_scan_redoCol"],
                                 config["annot_scan_missAutoCol"]), index=scan.index)
names(scan) <- c("scanID", "plate", "race", "redo", "missing")
racetbl <- table(scan$race)
majority <- names(racetbl)[racetbl == max(racetbl)]
batches <- ddply(scan, .(plate), summarise, nsamp=sum(!is.na(race)), nmaj=sum(race %in% majority), racefrac=nmaj/nsamp, meanmiss=mean(missing))
redoPlates <- unique(scan$plate[scan$redo %in% c("Y", "Yes", "yes", "YES", TRUE)])
batches$redo <- batches$plate %in% redoPlates

color_redo <- scale_color_manual(values=c("black", "red"), breaks=c("FALSE", "TRUE"))

if (type == "chisq"){
  batches$stat <- batch.res$mean.chisq[match(batches$plate, names(batch.res$mean.chisq))]
  ylab <- expression(paste("mean ", chi^2, " test statistic"))
  outfile <- config["out_meanchisq_race_plot"]
} else if (type == "fisher") {
  batches$stat <- batch.res$mean.or[match(batches$plate, names(batch.res$mean.or))]
  ylab <- "mean Fisher's OR"
  outfile <- config["out_meanor_race_plot"]
}
p <- ggplot(batches, aes(x=racefrac, y=stat, color=redo)) +
  geom_vline(x=mean(batches$racefrac), linetype='dashed') +  # mean over all plates
  geom_point() +
  color_redo +
  ylab(ylab) +
  xlab(paste("fraction of", majority, "samples per batch"))
ggsave(outfile, plot=p, width=6, height=6)





batches$lambda <- batch.res$lambda[match(batches$plate, names(batch.res$lambda))]
p <- ggplot(batches, aes(x=racefrac, y=lambda, color=redo)) +
  geom_vline(x=mean(batches$racefrac), linetype='dashed') + # mean over all plates
  geom_point() +
  color_redo +
  ylab(expression(paste("genomic inflation factor ", lambda))) +
  xlab(paste("fraction of", majority, "samples per batch"))
ggsave(config["out_lambda_race_plot"], plot=p, width=6, height=6)


# distribution of number of samples per batch
p <- ggplot(batches, aes(x=nsamp)) + geom_histogram(binwidth=1)
ggsave(config["out_hist_plot"], plot=p, width=6, height=6)


# mean autosomal missing call rate per batch
p <- ggplot(batches, aes(x=nsamp, y=meanmiss, color=redo)) +
  geom_point() +
  color_redo +
  ylab("mean autosomal missing call rate") +
  xlab("number of samples per batch")
ggsave(config["out_meanmcr_nscan_plot"], plot=p, width=6, height=6)

lmr <- lm(meanmiss ~ nsamp, data=batches)
anova(lmr)


# missing call rate vs test stat
if (type == "chisq") {
  xlab <- expression(paste("mean ", chi^2, " test statistic"))
  outfile <- config["out_meanmcr_meanchisq_plot"]
} else if (type == "fisher") {
  xlab <- "mean Fisher's OR"
  outfile <- config["out_meanmcr_meanor_plot"]
}
p <- ggplot(batches, aes(x=stat, y=meanmiss, color=redo)) +
  geom_point() +
  xlab(xlab) +
  ylab("mean autosomal missing call rate") +
  color_redo
ggsave(file=outfile, plot=p, widt=6, height=6)


# test stat vs number of scans
if (type == "chisq") {
  ylab <- expression(paste("mean ", chi^2, " test statistic"))
  outfile <- config["out_meanchisq_nscan_plot"]
} else if (type == "fisher") {
  ylab <- "mean Fisher's OR"
  outfile <- config["out_meanor_nscan_plot"]
}
p <- ggplot(batches, aes(x=nsamp, y=stat, color=redo)) +
  geom_point() +
  xlab("number of samples per batch") +
  ylab(ylab) +
  color_redo
ggsave(file=outfile, plot=p, widt=6, height=6)



# batch effect on missing call rate
# make plate names shorter for plotting
scan$redo <- scan$plate %in% batches$plate[batches$redo]
scan$batchLabel <- vapply(strsplit(as.character(scan$plate), "-"), function(x) x[[1]][1], "a")
model <- log10(missing) ~ as.factor(batchLabel)
lm.result <- lm(model, data=scan)
anova(lm.result)

p <- ggplot(scan, aes(x=batchLabel, y=missing, color=redo)) +
  geom_boxplot(varwidth=TRUE) +
  ylab("log10(autosomal missing call rate)") +
  xlab("") +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  color_redo
ggsave(file=config["out_mcr_plot"], plot=p, width=6, height=6)


# chromosome 1 intensity by batch
mninten <- getobj(config["inten_file"])
mninten <- mninten[[1]]
stopifnot(all(!(names(mninten[scan.index, "1"]) %in% scan.exclude)))
scan$mninten <- mninten[(match(scan$scanID, rownames(mninten))), "1"]

p <- ggplot(scan, aes(x=batchLabel, y=mninten, color=redo)) +
  geom_boxplot(varwidth=TRUE) +
  ylab("mean chromosome 1 intensity") +
  xlab("") +
  theme(axis.text.x=element_text(angle = 90, hjust = 0)) +
  color_redo
ggsave(file=config["out_inten_plot"], plot=p, width=6, height=6)
