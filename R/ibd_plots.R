##########
# IBD plots
# Usage: R --args config.file < ibd_plots.R
##########

library(GWASTools)
library(QCpipeline)
library(RColorBrewer)
library(ggplot2)
sessionInfo()


## read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

## check config and set defaults
required <- c("annot_scan_file")
optional <- c("annot_scan_subjectCol", "exp_rel_file", "ibd_method",
              "out_ibd_kc32_file", "out_ibd_con_file",
              "out_ibd_con_plot", "out_ibd_exp_plot", "out_ibd_obs_plot",
              "out_ibd_rel_file", "out_ibd_unexp_plot",
              "out_ibd_unobs_dup_file", "out_ibd_unobs_rel_file",
              "scan_ibd_include_file")
default <- c("subjectID", NA, "KING", "ibd_kc32.RData",
             "ibd_connectivity.RData", "ibd_connectivity.pdf",
             "ibd_expected.pdf", "ibd_observed.pdf", "ibd_obsrel.RData",
             "ibd_unexpected.pdf",
             "ibd_unobs_dup.RData", "ibd_unobs_rel.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

theme_set(theme_bw() + theme(legend.position=c(1, 1), legend.justification=c(1,1), legend.background = element_rect(colour = "black")))

ibd <- getobj(config["out_ibd_kc32_file"])
dim(ibd)

## if ibd has > 5000 rows, make png instead of pdf
if (nrow(ibd) > 5000) {
  plotfile <- function(filename, ...) {
    png(plotname(filename), width=650, height=650, ...)
    par(cex=1.5)
  }
  plotname <- function(filename){
    sub("pdf$", "png", filename)
  }
  
} else {
  plotfile <- function(filename, ...) {
    pdf(plotname(filename), width=6, height=6)
  }
  plotname <- function(filename) {
    filename
  }
}


(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)
samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_subjectCol"]))
names(samp) <- c("scanID", "Individ")
## keep only scans used in IBD
if (!is.na(config["scan_ibd_include_file"])) {
  scan.ids <- getobj(config["scan_ibd_include_file"])
  samp <- samp[samp$scanID %in% scan.ids,]
}
dim(samp)

## add subjectID so that we can bring in expected relation
ibd <- merge(ibd, samp, by.x="ID1", by.y="scanID")
ibd <- merge(ibd, samp, by.x="ID2", by.y="scanID", suffixes=c("1","2"))
if (config["ibd_method"] == "KING") {
  cols <- c("IBS0", "kinship")
} else {
  cols <- c("k0", "k1", "kinship")
}
ibd <- ibd[,c("ID1", "ID2", "Individ1", "Individ2", cols)]
ibd$ii <- pasteSorted(ibd$Individ1, ibd$Individ2)

if (!is.na(config["exp_rel_file"])) {
  relprs <- getobj(config["exp_rel_file"])
  message("expected relative pairs (subjects)")
  print(table(relprs$relation))
  print(table(relprs$exp.rel))
  relprs$ii <- pasteSorted(relprs$Individ1, relprs$Individ2)

  ibd <- merge(ibd, relprs[,c("ii", "relation", "exp.rel")], all.x=TRUE)
  ## names(ibd)[names(ibd) == "relation"] <- "exp.rel"
  message("expected relative pairs (samples)")
  print(table(ibd$exp.rel, useNA="ifany"))
} else {
  ibd$exp.rel <- NA
}

ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"
ibd$exp.rel[is.na(ibd$exp.rel)] <- "U"
message("expected relative pairs (including dups)")
table(ibd$exp.rel)


## make our own color code so we can modify the plot legend
rels <- c("Dup", "PO", "FS", "Deg1", "Deg2", "Deg3", "Q", "U")
cols <- c(brewer.pal(length(rels)-1, "Dark2")[c(1, 2, 3, 6, 5, 4, 7)], "black")
cmap <- setNames(cols, rels)


if (config["ibd_method"] == "KING") {

  ## thresholds for assigning relationships using kinship coefficients in table 1 of Manichaikul (2010) - KING paper
  cut.dup <- 1/(2^(3/2))
  cut.deg1 <- 1/(2^(5/2))
  cut.deg2 <- 1/(2^(7/2))
  cut.deg3 <- 1/(2^(9/2))
  cut.ibs <- 0.003 # should be 0 for PO, but sometimes is greater due to genotyping error. 0.003 works for OLGA, Kittner, and HRS2
  
  alpha <- 0.7
  
  p <- ggplot(ibd, aes(x=IBS0, y=kinship, color=exp.rel)) +
    geom_hline(y=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), linetype='dashed', color="grey") +
    geom_vline(x=cut.ibs, linetype='dashed', color="grey") +
    geom_point(alpha=alpha) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    guides(colour=guide_legend(override.aes=list(alpha=1))) +
    xlab("Fraction of IBS=0") + ylab("Kinship coefficient") +
    theme(legend.position=c(1, 1), legend.justification=c(1,1)) +
    ggtitle("IBD - expected")
  ggsave(plotname(config["out_ibd_exp_plot"]), plot=p, width=6, height=6)
  
  ## assign observed relationships (duplicates)
  ibd$obs.rel <- ibdAssignRelatednessKing(ibd$IBS0, ibd$kinship)
  message("observed relative pairs")
  print(table(ibd$obs.rel, useNA="ifany"))
  
  p <- ggplot(ibd, aes(x=IBS0, y=kinship, color=obs.rel)) +
    geom_hline(y=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), linetype='dashed', color="grey") +
    geom_vline(x=cut.ibs, linetype='dashed', color="grey") +
    geom_point(alpha=alpha) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    guides(colour=guide_legend(override.aes=list(alpha=1))) +
    xlab("Fraction of IBS=0") + ylab("Kinship coefficient") +
    ggtitle("IBD - observed")
  ggsave(plotname(config["out_ibd_obs_plot"]), plot=p, width=6, height=6)
  
  
  ibd$unexp <- ibd$exp.rel != ibd$obs.rel & ibd$kinship > cut.deg2 # use degree 2 cutoff in KING paper - ~0.088
  ## ## check for Deg2 and Deg3
  ## deg2 <- ibd$exp.rel %in% deg2.rel & ibd$obs.rel == "Deg2"
  ## deg3 <- ibd$exp.rel %in% deg3.rel & ibd$obs.rel == "Deg3"
  ## unexp <- unexp & !deg2 & !deg3
  message("unexpected relative pairs")
  print(table(ibd$obs.rel[ibd$unexp]))
  
  ibd<-ibd[order(ibd$unexp),] # order the file before plotting to highlight the unexpected relationships
  p <- ggplot(ibd, aes(x=IBS0, y=kinship, color=exp.rel, pch=unexp)) +
    geom_hline(y=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), linetype='dashed', color="grey") +
    geom_vline(x=cut.ibs, linetype='dashed', color="grey") +
    geom_point(alpha=alpha) +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    guides(colour=guide_legend(override.aes=list(alpha=1))) +
    xlab("Fraction of IBS=0") + ylab("Kinship coefficient") +
    ggtitle("IBD - unexpected")
  ggsave(plotname(config["out_ibd_unexp_plot"]), plot=p, width=6, height=6)
  
} else {
  plotfile(config["out_ibd_exp_plot"])
  ibdPlot(ibd$k0, ibd$k1, relation=ibd$exp.rel, main="IBD - expected")
  dev.off()

  ## assign observed relationships
  ibd$obs.rel <- ibdAssignRelatedness(ibd$k0, ibd$k1)
  message("observed relative pairs")
  print(table(ibd$obs.rel))

  plotfile(config["out_ibd_obs_plot"])
  ibdPlot(ibd$k0, ibd$k1, relation=ibd$obs.rel, main="IBD - observed")
  dev.off()

  ## plot of unexpected relationships (kinship > 0.1)
  unexp <- ibd$exp.rel != ibd$obs.rel & ibd$kinship > 0.09833927
  ## ## check for Deg2 and Deg3
  ## deg2 <- ibd$exp.rel %in% deg2.rel & ibd$obs.rel == "Deg2"
  ## deg3 <- ibd$exp.rel %in% deg3.rel & ibd$obs.rel == "Deg3"
  ## unexp <- unexp & !deg2 & !deg3
  message("unexpected relative pairs")
  print(table(ibd$obs.rel[unexp]))

  
  plotfile(config["out_ibd_unexp_plot"])
  psym <- rep(1, nrow(ibd))
  psym[unexp] <- 2
  cmap[c("Deg3", "Q", "U")] <- "black"
  pcol <- cmap[ibd$exp.rel]
  ibdPlot(ibd$k0, ibd$k1, color=pcol, pch=psym)
  points(ibd$k0[unexp], ibd$k1[unexp], col=pcol[unexp], pch=psym[unexp])
  rel <- rels[rels %in% unique(ibd$exp.rel)]
  col <- cmap[rel]
  sym <- c(rep(1, length(col)), 2)
  rel[rel == "U"] <- "Unrel"
  rel <- paste("Exp", rel)
  rel <- c(rel, "Unexp")
  col <- c(col, "black")
  legend("topright", legend=rel, col=col, pch=sym)
  dev.off()
}


## check for expected relationships not observed
if (!is.na(config["exp_rel_file"])) {
  relprs <- relprs[relprs$exp.rel != "U",]
  unobs.sel <- !(relprs$ii %in% ibd$ii)
  if (sum(unobs.sel) > 0) {
    message(paste(sum(unobs.sel), "relative pairs not observed"))
    unobs.rel <- relprs[unobs.sel,]
    unobs.rel$ii <- NULL
    save(unobs.rel, file=config["out_ibd_unobs_rel_file"])
  } else {
    message("all expected relatives observed")
  }
}

## check for expected duplicates not observed
dupsubj <- unique(samp$Individ[duplicated(samp$Individ)])
unobs.dup <- list()
for (d in dupsubj) {
  exp.scans <- samp$scanID[samp$Individ == d]
  this.subj <- ibd$Individ1 == d & ibd$exp.rel == "Dup"
  obs.set <- c(ibd$ID1[this.subj], ibd$ID2[this.subj])
  unobs.scans <- setdiff(exp.scans, obs.set)
  if (length(unobs.scans) > 0) {
    unobs.dup[[as.character(d)]] <- unobs.scans
  }
}
if (length(unobs.dup) > 0) {
  message(paste(length(unobs.dup), "duplicate pairs not observed"))
  save(unobs.dup, file=config["out_ibd_unobs_dup_file"])
} else {
  message("all expected duplicates observed")
}

ibd$ii <- NULL
save(ibd, file=config["out_ibd_rel_file"])


## IBD connectivity
dupids <- ibd$ID2[ibd$obs.rel %in% "Dup"]
length(dupids)
ibd.nodup <- ibd[!(ibd$ID1 %in% dupids) & !(ibd$ID2 %in% dupids),]
(npr <- nrow(ibd.nodup))

allsamp <- c(ibd.nodup$ID1, ibd.nodup$ID2)
uniqsamp <- unique(allsamp)
(n <- length(uniqsamp))
con <- rep(NA, n)
for(i in 1:n) con[i] <- sum(allsamp %in% uniqsamp[i])

dat <- data.frame(con=sort(con), rank=1:length(con))
p <- ggplot(dat, aes(x=rank, y=con)) +
  geom_point() +
  ylab("sample connectivity") +
  ggtitle(paste(npr, "pairs with KC > 1/32"))
ggsave(plotname(config["out_ibd_con_plot"]), plot=p, width=6, height=6)


con.df <- data.frame("scanID"=uniqsamp, "connectivity"=con)
save(con.df, file=config["out_ibd_con_file"])
