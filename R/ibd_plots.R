##########
# IBD plots
# Usage: R --args config.file < ibd_plots.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "out_ibd_kc32_file")
optional <- c("annot_scan_subjectCol", "exp_rel_file", "out_ibd_con_file",
              "out_ibd_con_plot", "out_ibd_exp_plot", "out_ibd_obs_plot",
              "out_ibd_rel_file", "out_ibd_unexp_plot",
              "out_ibd_unobs_dup_file", "out_ibd_unobs_rel_file",
              "scan_ibd_include_file")
default <- c("subjectID", NA, "ibd_connectivity.RData", "ibd_connectivity.pdf",
             "ibd_expected.pdf", "ibd_observed.pdf", "ibd_obsrel.RData",
             "ibd_unexpected.pdf",
             "ibd_unobs_dup.RData", "ibd_unobs_rel.RData", NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

ibd <- getobj(config["out_ibd_kc32_file"])
dim(ibd)

## if ibd has > 5000 rows, make png instead of pdf
if (nrow(ibd) > 5000) {
  plotfile <- function(filename, ...) {
    file <- sub("pdf$", "png", filename)
    png(file, width=650, height=650, ...)
    par(cex=1.5)
  }
} else {
  plotfile <- function(filename, ...) {
    pdf(filename, width=6, height=6)
  }
}

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)
samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_subjectCol"]))
names(samp) <- c("scanID", "Individ")
# keep only scans used in IBD
if (!is.na(config["scan_ibd_include_file"])) {
  scan.ids <- getobj(config["scan_ibd_include_file"])
  samp <- samp[samp$scanID %in% scan.ids,]
}
dim(samp)

# add subjectID so that we can bring in expected relation
ibd <- merge(ibd, samp, by.x="ID1", by.y="scanID")
ibd <- merge(ibd, samp, by.x="ID2", by.y="scanID", suffixes=c("1","2"))
ibd <- ibd[,c("ID1", "ID2", "Individ1", "Individ2", "k0", "k1", "kinship")]
ibd$ii <- paste(pmin(ibd$Individ1, ibd$Individ2),
                pmax(ibd$Individ1, ibd$Individ2))

if (!is.na(config["exp_rel_file"])) {
  relprs <- getobj(config["exp_rel_file"])
  print(table(relprs$relation))
  relprs$ii <- paste(pmin(relprs$Individ1, relprs$Individ2),
                     pmax(relprs$Individ1, relprs$Individ2))

  ibd <- merge(ibd, relprs[,c("ii","relation")], all.x=TRUE)
  names(ibd)[names(ibd) == "relation"] <- "exp.rel"
  print(table(ibd$exp.rel, useNA="ifany"))
} else {
  ibd$exp.rel <- NA
}

ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"
ibd$exp.rel[is.na(ibd$exp.rel)] <- "U"
table(ibd$exp.rel)

plotfile(config["out_ibd_exp_plot"])
ibdPlot(ibd$k0, ibd$k1, relation=ibd$exp.rel, main="IBD - expected")
dev.off()

# assign observed relationships
ibd$obs.rel <- ibdAssignRelatedness(ibd$k0, ibd$k1)
table(ibd$obs.rel)

plotfile(config["out_ibd_obs_plot"])
ibdPlot(ibd$k0, ibd$k1, relation=ibd$obs.rel, main="IBD - observed")
dev.off()

# plot of unexpected relationships (kinship > 0.1)
unexp <- ibd$exp.rel != ibd$obs.rel & ibd$kinship > 0.09833927
# check for Deg2 and Deg3
deg2.rel <-  c("HS", "HSr", "HSFC", "Av", "GpGc", "DFC")
deg2 <- ibd$exp.rel %in% deg2.rel & ibd$obs.rel == "Deg2"
deg3.rel <- c("FC", "HAv", "OAv", "OC")
deg3 <- ibd$exp.rel %in% deg3.rel & ibd$obs.rel == "Deg3"
unexp <- unexp & !deg2 & !deg3
table(ibd$obs.rel[unexp])

plotfile(config["out_ibd_unexp_plot"])
psym <- rep(1, nrow(ibd))
psym[unexp] <- 2
# make our own color code so we can modify the plot legend
prel <- ibd$exp.rel
prel[prel %in% deg2.rel] <- "Deg2"
pcol <- rep("black", nrow(ibd))
pcol[prel == "Dup"] <- "magenta"
pcol[prel == "PO"] <- "cyan"
pcol[prel == "FS"] <- "red"
pcol[prel == "Deg2"] <- "blue"
ibdPlot(ibd$k0, ibd$k1, color=pcol, pch=psym, rel.draw=c("FS", "Deg2"))
rel <- unique(prel)
col <- unique(pcol)
sym <- c(rep(1, length(rel)), 2)
ord <- order(rel)
rel[rel == "U"] <- "Unrel"
rel <- paste("Exp", rel)
rel <- c(rel[ord], "Unexp")
col <- c(col[ord], "black")
legend("topright", legend=rel, col=col, pch=sym)
dev.off()


# check for expected relationships not observed
if (!is.na(config["exp_rel_file"])) {
  relprs <- relprs[relprs$relation != "U",]
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
# check for expected duplicates not observed
dupsubj <- unique(samp$subjectID[duplicated(samp$subjectID)])
unobs.dup <- list()
for (d in dupsubj) {
  exp.scans <- samp$scanID[samp$subjectID == d]
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


# IBD connectivity
dupids <- ibd$ID2[ibd$obs.rel == "Dup"]
length(dupids)
ibd.nodup <- ibd[!(ibd$ID1 %in% dupids) & !(ibd$ID2 %in% dupids),]
(npr <- nrow(ibd.nodup))

allsamp <- c(ibd.nodup$ID1, ibd.nodup$ID2)
uniqsamp <- unique(allsamp)
(n <- length(uniqsamp))
con <- rep(NA, n)
for(i in 1:n) con[i] <- sum(allsamp %in% uniqsamp[i])

plotfile(config["out_ibd_con_plot"])
plot(sort(con), xlab="rank", ylab="sample connectivity", main=paste(npr, "pairs with KC > 1/32"))
dev.off()

con.df <- data.frame("scanID"=uniqsamp, "connectivity"=con)
save(con.df, file=config["out_ibd_con_file"])
