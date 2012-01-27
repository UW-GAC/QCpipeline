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
print(config)

ibd <- getobj(config["out_ibd_kc32_file"])
dim(ibd)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)
samp <- getVariable(scanAnnot, c("scanID", config["annot_scan_subjectCol"]))
names(samp) <- c("scanID", "subjectID")
# keep only scans used in IBD
if (!is.na(config["scan_ibd_include_file"])) {
  scan.ids <- getobj(config["scan_ibd_include_file"])
  samp <- samp[samp$scanID %in% scan.ids,]
}
dim(samp)

# add subjectID so that we can bring in expected relation
ibd <- merge(ibd, samp, by.x="sample1", by.y="scanID")
names(ibd)[names(ibd) == "subjectID"] <- "Individ1"
ibd <- merge(ibd, samp, by.x="sample2", by.y="scanID")
names(ibd)[names(ibd) == "subjectID"] <- "Individ2"
ibd$ii <- paste(ibd$Individ1, ibd$Individ2)

if (!is.na(config["exp_rel_file"])) {
  relprs <- getobj(config["exp_rel_file"])
  print(table(relprs$relation))
  relprs$i12 <- paste(relprs$Individ1, relprs$Individ2)
  relprs$i21 <- paste(relprs$Individ2, relprs$Individ1)

  ibd <- merge(ibd, relprs[,c("i12","relation")], by.x="ii", by.y="i12", all.x=TRUE)
  names(ibd)[names(ibd) == "relation"] <- "rel"
  ibd <- merge(ibd, relprs[,c("i21","relation")], by.x="ii", by.y="i21", all.x=TRUE)
  names(ibd)[names(ibd) == "relation"] <- "rel2"
  ibd$exp.rel <- ibd$rel
  ibd$exp.rel[is.na(ibd$exp.rel)] <- ibd$rel2[is.na(ibd$exp.rel)]
  print(table(ibd$exp.rel, exclude=NULL))
  ibd$rel <- NULL
  ibd$rel2 <- NULL
} else {
  ibd$exp.rel <- NA
}

ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"
ibd$exp.rel[is.na(ibd$exp.rel)] <- "U"
table(ibd$exp.rel)

pdf(config["out_ibd_exp_plot"], width=6, height=6)
ibdPlot(ibd$k0, ibd$k1, relation=ibd$exp.rel, main="IBD - expected")
dev.off()

# assign observed relationships
ibd$obs.rel <- ibdAssignRelatedness(ibd$k0, ibd$k1)
table(ibd$obs.rel)

pdf(config["out_ibd_obs_plot"], width=6, height=6)
ibdPlot(ibd$k0, ibd$k1, relation=ibd$obs.rel, main="IBD - observed")
dev.off()


# check for expected relationships not observed
if (!is.na(config["exp_rel_file"])) {
  relprs <- relprs[relprs$relation != "U",]
  unobs.sel <- !(relprs$i12 %in% ibd$ii | relprs$i21 %in% ibd$ii)
  if (sum(unobs.sel) > 0) {
    message(paste(sum(unobs.sel), "relative pairs not observed"))
    unobs.rel <- relprs[unobs.sel,]
    unobs.rel$i12 <- NULL
    unobs.rel$i21 <- NULL
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
  obs.set <- c(ibd$sample1[this.subj], ibd$sample2[this.subj])
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
dupids <- ibd$sample2[ibd$obs.rel == "Dup"]
length(dupids)
ibd.nodup <- ibd[!(ibd$sample1 %in% dupids) & !(ibd$sample2 %in% dupids),]
(npr <- nrow(ibd.nodup))

allsamp <- c(ibd.nodup$sample1, ibd.nodup$sample2)
uniqsamp <- unique(allsamp)
(n <- length(uniqsamp))
con <- rep(NA, n)
for(i in 1:n) con[i] <- sum(allsamp %in% uniqsamp[i])

pdf(config["out_ibd_con_plot"], width=6, height=6)
plot(sort(con), xlab="rank", ylab="sample connectivity", main=paste(npr, "pairs with KC > 1/32"))
dev.off()

con.df <- data.frame("scanID"=uniqsamp, "connectivity"=con)
save(con.df, file=config["out_ibd_con_file"])
