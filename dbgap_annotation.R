# make sample annotation files for dbGaP
# arguments:
#  annot - dataframe with sample annotation
#  vardata - dataframe with columns "varname", "newname", "annotation",
#     "analysis," "description", "type"
#  dir - directory for output (default current)

dbgap.sample.annotation <- function(annot, vardata, dir=".") {

  # check for essential columns
  stopifnot("consent" %in% names(annot),
            "subj.plink" %in% names(annot),
            "dup.post" %in% names(annot),
            "no.post" %in% names(annot))
  stopifnot("consent" %in% vardata$varname,
            "subj.plink" %in% vardata$varname,
            "dup.post" %in% vardata$varname,
            "no.post" %in% vardata$varname)
  
  # consent levels
  conslev <- unique(annot$consent)

  # annotation and analysis
  for (type in c("annotation", "analysis")) {
  
    # output files
    annotfile <- paste(dir, "/Sample_", type, ".csv", sep="")
    dupfile <- paste(dir, "/Sample_", type, "_duplicates.csv", sep="")
    ddfile <- paste(dir, "/Sample_", type, "_DD.txt", sep="")
    consent.annotfiles <- paste(dir, "/Sample_", type, "_consent_", conslev, ".csv", sep="")
    consent.dupfiles <- paste(dir, "/Sample_", type, "_duplicates_consent_", conslev, ".csv", sep="")

    # variables to include in sample annotation file
    sann <- vardata$varname[vardata[,type] == 1]
    nsann <- vardata$newname[vardata[,type] == 1]

    # remove any samples that will not be posted on dbGaP
    annot2 <- annot[!annot$no.post, sann]

    # change var names to a more generic form
    names(annot2) <- nsann

    # split the sample annotation into the main, unduplicated set
    # (one sample per subject) and the duplicated set
    subj <- annot2[annot2$subj.plink,]
    dups <- annot2[annot2$dup.post,]
    stopifnot(nrow(dups) + nrow(subj) == nrow(annot2))
  
    # remove the splitting variable(s)
    x <- which(is.element(names(subj),c("subj.plink","dup.post","no.post")))
    subj <- subj[,-x] 
    x <- which(is.element(names(dups),c("subj.plink","dup.post","no.post")))
    dups <- dups[,-x]

    # prepare data dictionary
    dd <- vardata[vardata[,type] == 1, c("newname", "type", "description")]
    dd <- dd[!is.element(dd$newname,c("subj.plink","dup.post","no.post")),]
    stopifnot(names(subj) == dd$newname)
    names(dd)[1] <- "variable"

    # write the files
    write.csv(subj, file=annotfile, quote=FALSE, row.names=FALSE)
    write.csv(dups, file=dupfile, quote=FALSE, row.names=FALSE)
    write.table(dd, file=ddfile, sep="\t", quote=FALSE, row.names=FALSE)

    # Divide samples by consent group
    if (type == "annotation") {
      for (i in 1:length(conslev)) {
        s <- subj[subj$consent == conslev[i],]
        write.csv(s, file=consent.annotfiles[i], quote=FALSE, row.names=FALSE)
      }

      for (i in 1:length(conslev)) {
        s <- dups[dups$consent == conslev[i],]
        write.csv(s, file=consent.dupfiles[i], quote=FALSE, row.names=FALSE)
      }
    }
  }
}


# make snp annotation files for dbGaP
# arguments:
#  annot - dataframe with snp annotation
#  vardata - dataframe with columns "varname", "newname", "annotation",
#     "analysis," "description", "type"
#  dir - directory for output (default current)

dbgap.snp.annotation <- function(annot, vardata, dir=".") {

  # annotation and analysis
  for (type in c("annotation", "analysis")) {
    
    # output files
    annotfile <- paste(dir, "/SNP_", type, ".csv", sep="")
    ddfile <- paste(dir, "/SNP_", type, "_DD.txt", sep="")

    # variables to include in sample annotation file
    sann <- vardata$varname[vardata[,type] == 1]
    nsann <- vardata$newname[vardata[,type] == 1]

    # keep selected variables
    annot2 <- annot[, sann]

    # change var names to a more generic form
    names(annot2) <- nsann

    # prepare data dictionary
    dd <- vardata[vardata[,type] == 1, c("newname", "type", "description")]
    stopifnot(names(annot2) == dd$newname)
    names(dd)[1] <- "variable"

    # write the files
    write.csv(annot2, file=annotfile, quote=FALSE, row.names=FALSE)
    write.table(dd, file=ddfile, sep="\t", quote=FALSE, row.names=FALSE)
  }
}
