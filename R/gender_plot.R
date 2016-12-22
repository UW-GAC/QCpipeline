##########
# Plot X and Y intensity and X heterozygosity for gender check
# Usage: R --args config.file < gender_plot.R
##########

library(GWASTools)
library(QCpipeline)
library(ggplot2)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "annot_snp_file")
optional <- c("annot_scan_sexCol", "out_pdf", "out_autosome_prefix",
              "out_het_file", "out_inten_file")
default <- c("sex", "sex_check.pdf", "autosome",
             "het_by_scan.RData", "mean_inten.RData")
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))

# mean intensity
mninten <- getobj(config["out_inten_file"])
mninten <- mninten[[1]] # only want sum of X and Y intensities
stopifnot(allequal(getScanID(scanAnnot), rownames(mninten)))
colnames(mninten) <- paste0("mninten.", colnames(mninten))

# heterozygosity
het <- getobj(config["out_het_file"])
stopifnot(allequal(getScanID(scanAnnot), rownames(het)))
colnames(het) <- paste0("het.", colnames(het))

dat <- getVariable(scanAnnot, c("scanID", config["annot_scan_sexCol"]))
names(dat)[2] <- "sex.annot"

dat <- cbind(dat, mninten)
dat <- cbind(dat, het)

# best guess at real sex
dat$sex.genotype <- ifelse(dat$mninten.Y > 0.5 & dat$mninten.X < 1.0, "M", "F")

# order so that mismatches are on top
dat$mismatch <- ifelse(dat$sex.genotype == dat$sex.annot, 0, 1)
dat$mismatch[is.na(dat$mismatch)] <- 1
dat <- dat[order(dat$mismatch),]

# replace so it shows up on the axis label
dat$sex.annot[is.na(dat$sex.annot)] <- "NA"

# plot labels - X and Y probes
(snpAnnot <- getobj(config["annot_snp_file"]))
chr <- getChromosome(snpAnnot, char=TRUE)
nx <- sum(chr == "X")
ny <- sum(chr == "Y")
xlab <- paste("X intensity (",nx," probes)", sep="")
ylab <- paste("Y intensity (",ny," probes)", sep="")


plots <- list()
theme_set(theme_bw())

ggcol <- scale_color_manual(breaks=c("F", "M", "NA"), values=c("#FF4D4D", "#0099FF", "black"))

plots[[1]] <- ggplot(dat, aes(x=mninten.X, y=mninten.Y, color=sex.annot)) +
  geom_point() +
  xlab(xlab) +
  ylab(ylab) +
  ggcol + 
  guides(alpha=FALSE) +
  theme(legend.position=c(1,1), legend.justification=c(1,1))

plots[[2]] <- ggplot(dat, aes(x=mninten.X, y=het.X, color=sex.annot)) +
  geom_point() +
  xlab(xlab) +
  ylab("X heterozygosity") +
  ggcol + 
  guides(color=FALSE, alpha=FALSE)

plots[[3]] <- ggplot(dat, aes(x=mninten.Y, y=het.X, color=sex.annot)) +
  geom_point() +
  xlab(ylab) +
  ylab("X heterozygosity") +
  ggcol +
  guides(color=FALSE, alpha=FALSE)

plots[[4]] <- ggplot(dat[dat$sex.annot %in% "F", ], aes(x=het.A, y=het.X, color=sex.annot)) +
  geom_point() +
  xlab("Autosomal heterozygosity") +
  ylab("X heterozygosity") +
  ggcol +
  guides(color=FALSE, alpha=FALSE)

pdf(file=config["out_pdf"], width=8, height=8)
multiplot(plotlist=plots, cols=2, byrow=TRUE)
dev.off()


# plot autosome intensities
autoPrefix <- config["out_autosome_prefix"]

autosomes <- intersect(unique(chr), as.character(1:22))
for (autosome in autosomes){
  ni <- sum(chr %in% autosome)
  ylab <- paste("Chr ", autosome, " intensity (", ni, " probes)", sep="")
  
  png(paste(autoPrefix, "_", autosome, ".png", sep=""), width=600, height=600)
  
  p <- ggplot(dat, aes_string(x="mninten.X", y=paste0("mninten.", autosome), color="sex.annot")) +
    geom_point() +
    xlab(xlab) + ylab(ylab) +
    ggcol +
    theme(legend.position=c(1,1), legend.justification=c(1,1))
  print(p)
  ggsave(paste(autoPrefix, "_", autosome, ".png", sep=""), width=8, height=8)
}
