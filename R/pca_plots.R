##########
# PCA plots
# Usage: R --args config.file pca.type < pca_plots.R
##########

library(GWASTools)
library(QCpipeline)
library(RColorBrewer)
library(MASS)
library(ggplot2)
library(GGally)
library(ggExtra)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check for type
if (length(args) < 2) stop("missing pca type (study or combined)")
type <- args[2]

theme_set(theme_bw())


# check config and set defaults
if (type == "study") {
    required <- c("annot_scan_file", "annot_scan_raceCol", "annot_snp_file")
    optional <- c("annot_scan_ethnCol", "annot_snp_rsIDCol",
                  "num_evs_to_plot", "out_corr_file", "out_pca_file",
                  "out_corr_plot_prefix", "out_corr_pruned_plot_prefix",
                  "out_dens_plot", "out_ev12_plot", "out_pairs_plot", "out_scree_plot",
                  "out_parcoord_plot",
                  "parcoord_vars", "out_parcoord_var_prefix")
    default <- c(NA, "rsID", 12, "pca_corr.RData", "pca.RData", 
                 "pca_corr", NA, "pca_dens.pdf",
                 "pca_ev12.pdf", "pca_pairs.png", "pca_scree.pdf",
                 "pca_parcoord.png",
                 "", "pca_parcoord")
    snpfile <- config["annot_snp_file"]
} else if (type == "combined"){
    required <- c("annot_scan_file", "annot_scan_raceCol", "out_comb_prefix")
    optional <- c("annot_scan_ethnCol", "annot_snp_rsIDCol", "ext_annot_scan_file",
                  "ext_annot_scan_raceCol", 
                  "num_evs_to_plot", "out_corr_file", "out_pca_file",
                  "out_corr_plot_prefix", "out_corr_pruned_plot_prefix",
                  "out_dens_plot", "out_ev12_plot", "out_pairs_plot", "out_scree_plot",
                  "out_parcoord_plot",
                  "out_ev12_plot_hapmap", "out_ev12_plot_study",
                  "parcoord_vars", "out_parcoord_var_prefix")
    default <- c(NA, "rsID", NA, "pop.group", 12, "pca_combined_corr.RData",
                 "pca_combined.RData", "pca_corr", NA, "pca_dens.pdf",
                 "pca_ev12.pdf", "pca_pairs.png", "pca_scree.pdf",
                 "pca_parcoord.png",
                 "pca_ev12_hapmap.pdf", "pca_ev12_study.pdf",
                 "", "pca_parcoord")
    snpfile <- paste0(config["out_comb_prefix"], "_snpAnnot.RData")
}
config <- setConfigDefaults(config, required, optional, default)
print(config)


# functions for parallel coordinate plots later
.getN <- function(samp, var){
  ntab <- table(samp[[var]])
  n <- ntab[as.character(samp[[var]])]
  n[is.na(n)] <- sum(is.na(samp[[var]]))
  n
}
# transparency based on number of samples in a group
.getParcoordAlpha <- function(samp, var) {
  n <- .getN(samp, var)
  return(ifelse(n < 10, 1,
                ifelse(n < 100, 0.5,
                       ifelse(n < 1000, 0.3, 0.1))) * 255) 
}

# parallel coordinates plot variables
vars <- unlist(strsplit(config["parcoord_vars"], " "))

# scan annotation
if (type == "study") {
  scanAnnot <- getobj(config["annot_scan_file"])
  samp <- getVariable(scanAnnot, c("scanID", c(config["annot_scan_raceCol"], vars)))
  names(samp) <- c("scanID", "race", vars)
  if (!is.na(config["annot_scan_ethnCol"])) {
    samp$ethnicity <- getVariable(scanAnnot, config["annot_scan_ethnCol"])
  } else samp$ethnicity <- NA
} else if (type == "combined") {
  scanAnnot <- getobj(config["annot_scan_file"])
  scan1 <- getVariable(scanAnnot, c("scanID", config["annot_scan_raceCol"], config["annot_scan_hapmapCol"]))
  names(scan1) <- c("scanID", "race", "geno.cntl")
  if (!is.na(config["annot_scan_ethnCol"])) {
    scan1$ethnicity <- getVariable(scanAnnot, config["annot_scan_ethnCol"])
  } else scan1$ethnicity <- NA
  if (sum(is.na(scan1$race)) > 0 & hasVariable(scanAnnot, config["ext_annot_scan_raceCol"])) {
    scan1$race2 <- getVariable(scanAnnot, config["ext_annot_scan_raceCol"])
    scan1$race[is.na(scan1$race)] <- scan1$race2[is.na(scan1$race)]
    scan1$race2 <- NULL
  }
  ext.scanAnnot <- getobj(config["ext_annot_scan_file"])
  scan2 <- getVariable(ext.scanAnnot, c("scanID", config["ext_annot_scan_raceCol"]))
  names(scan2) <- c("scanID", "race")
  scan2$geno.cntl <- 1
  scan2$ethnicity <- NA
  samp <- rbind(scan1, scan2)
} else {
  stop("pca type must be study or combined")
}

# get PCA results
pca <- getobj(config["out_pca_file"])
samp <- samp[match(pca$sample.id, samp$scanID),]
stopifnot(allequal(pca$sample.id, samp$scanID))
table(samp$race, samp$ethnicity, useNA="ifany")

# why are we doing this? (sort capital letters before lower case)
Sys.setlocale("LC_COLLATE", "C")

# color by race
race <- as.character(sort(unique(samp$race)))
if (length(race) > 0) {
  #stopifnot(all(race %in% names(config)))
  
  cmapRace <- setNames(config[race], race)
  chk <- which(is.na(cmapRace))
  if (length(chk) > 0) {
    message(sprintf("Using default colors for %s races: %s", length(chk), paste(names(cmapRace[chk]), collapse=", ")))
    defaultColors <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
    cmapRace[chk] <- defaultColors[1:length(chk)]
  }
  
  colorScale <- scale_color_manual("race", values=cmapRace, breaks=names(cmapRace), na.value="grey")
} else {
  colorScale <- scale_color_manual("race", values="black", breaks="hack", na.value="black")
}
rm(race)

# plot symbol by ethnicity
ethn <- as.character(sort(unique(samp$ethnicity)))
if (length(ethn) > 0){
  stopifnot(all(ethn %in% names(config)))
  symbolMap <- config[ethn]
  mode(symbolMap) <- "integer"
  symbolScale <- scale_shape_manual("ethnicity", values=symbolMap, breaks=names(symbolMap), na.value=16)
} else {
  symbolScale <- scale_shape_manual("ethnicity", values=1, breaks="hack", na.value=16)
}
rm(ethn)


# labels
## recent change in SNPRelate - pca$eigenval only returns first 32 values
#(x <- pca$eigenval[1:4]/sum(pca$eigenval))
x <- pca$varprop[1:4]
lbls <- paste("EV", 1:4, " (", format(100*x,digits=2), "%)", sep="")

#samp$nrace <- table(samp$race, useNA="ifany")[samp$race]
#samp$nrace[is.na(samp$nrace)] <- sum(is.na(samp$nrace))
#zorder <- order(-samp$nrace)

pcs <- pca$eigenvect
colnames(pcs) <- paste0("EV", 1:ncol(pcs))
pcs <- as.data.frame(pcs)
pcs$scanID <- pca$sample.id

dat <- merge(pcs, samp)

# plot the first four pcs
nev <- 4
pairs <- ggpairs(dat,
                 columns=which(names(dat) %in% sprintf("EV%s", 1:nev)),
                 upper=list(continuous="points"),
                 color="race",
                 pch="ethnicity",
                 columnLabels=lbls[1:nev],
                 axisLabels="internal",
                 params=c(alpha=0))
for (i in 1:length(pairs$columns)){
  for (j in 1:length(pairs$columns)){
    subplot <- getPlot(pairs, i, j)
    if (i != j){
      subplot <- subplot + geom_point(alpha=0.7, size=2, aes(pch=ethnicity))
    } else {
    }
    subplot <- subplot + colorScale + symbolScale
    #subplot <- subplot + theme()
    #subplot <- subplot + guides(colour = guide_legend(override.aes = list(size=3)))
    #if (i == 2 & j == 1) leg <- g_legend(subplot)
    pairs <- putPlot(pairs, subplot, i, j)
  }
}
png(config["out_pairs_plot"], width=720, height=720)
print(pairs)
dev.off()


# plot EV1 vs EV2 with density plots
p <- ggplot(dat, aes(x=EV1, y=EV2, color=race, pch=ethnicity)) +
  geom_point(alpha=0.7) +
  colorScale +
  symbolScale +
  theme(legend.position="none") +
  xlab(lbls[1]) + ylab(lbls[2])
pdf(config["out_dens_plot"], width=6, height=6)
ggMarginal(p, type="density")
dev.off()


# plot EV1 vs EV2
p <- p + theme(legend.position="right")
ggsave(config["out_ev12_plot"], plot=p, width=6, height=6)


ggParcoordTheme <- theme(axis.title.x=element_blank(),
                         axis.ticks.y=element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.y=element_blank(),
                         axis.text.x=element_text(colour="black"),
                         panel.grid.major.y=element_blank(),
                         panel.grid.minor.y=element_blank(),
                         legend.position="top")


# parallel coordinates plot
ev.ind <- which(names(dat) %in% sprintf("EV%s", 1:12))
p <- ggparcoord(dat, columns=ev.ind, groupColumn="race", scale="uniminmax", alphaLines=0.5) +
  colorScale + ggParcoordTheme
ggsave(config["out_parcoord_plot"], plot=p, width=10, height=5)


## other variables for parallel coordinate, specified by user
if (type == "study" & length(vars) > 0){
  for (var in vars){
    stopifnot(var %in% names(samp))
    
    # auto filename
    fname <- paste(config["out_parcoord_var_prefix"], "_", var, ".png", sep="")
    dat[["fvar"]] <- as.factor(dat[[var]])
    
    p <- ggparcoord(dat, columns=ev.ind, groupColumn="fvar", scale="uniminmax", alphaLines=0.5) +
      scale_color_brewer(var, palette="Set1", na.value="grey") + ggParcoordTheme
    ggsave(fname, plot=p, width=10, height=5)
    
  }
}

if (type == "combined"){
  
  xlim <- range(dat$EV1)
  ylim <- range(dat$EV2)
  
  # hapmap plot
  dat$plotcol <- dat$race
	dat$plotcol[dat$geno.cntl %in% 0] <- NA
  p <- ggplot(dat, aes(x=EV1, y=EV2, color=plotcol, pch=ethnicity)) +
    geom_point() +
    colorScale +
    symbolScale +
    xlab(lbls[1]) + ylab(lbls[2]) +
    xlim(xlim)
  ggsave(config["out_ev12_plot_hapmap"], plot=p, width=6, height=6)
  	
  p <- ggplot(dat[dat$geno.cntl %in% 0, ], aes(x=EV1, y=EV2, color=race, pch=ethnicity)) +
    geom_point() +
    colorScale +
    symbolScale +
    xlab(lbls[1]) + ylab(lbls[2]) +
    xlim(xlim)
  ggsave(config["out_ev12_plot_study"], plot=p, width=6, height=6)

}


#plot SNP-PC correlation
snpAnnot <- getobj(snpfile)
corr <- getobj(config["out_corr_file"])
snp <- snpAnnot[match(corr$snp.id, getSnpID(snpAnnot)),]
chrom <- getChromosome(snp, char=TRUE)

nev <- as.integer(config["num_evs_to_plot"])

png(paste(config["out_corr_plot_prefix"], "_%03d.png", sep=""), height=720, width=720)
par(mfrow=c(4,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
for(i in 1:nev){
  snpCorrelationPlot(abs(corr$snpcorr[i,]), chrom,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()

if (!is.na(config["out_corr_pruned_plot_prefix"])) {
  snps.pruned <- getobj(config["out_pruned_file"])
  ind <- getSnpID(snp) %in% snps.pruned
  
  png(paste(config["out_corr_pruned_plot_prefix"], "_%03d.png", sep=""), height=720, width=720)
  par(mfrow=c(4,1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
  for(i in 1:nev){
    snpCorrelationPlot(abs(corr$snpcorr[i,ind]), chrom[ind],
                       main=paste("Eigenvector",i), ylim=c(0,1))
  }
  dev.off()
}
  
# scree plot
dat <- data.frame(ev=1:nev, varprop=pca$varprop[1:nev])

p <- ggplot(dat, aes(x=factor(ev), y=100*varprop)) +
  geom_point() +
  xlab("Eigenvector") + ylab("Percent of variance accounted for")
ggsave(config["out_scree_plot"], plot=p, width=6, height=6)

