# make sample annotation files for dbGaP
# arguments:
#  annot - ScanAnnotationDataFrame
#  dir - directory for output (default current)

dbgapScanAnnotation <- function(scanAnnot, dir=".",
   consentVar="consent", subjVar="subj.plink", dupVar="dup.post", omitVar="no.post",
   annotationCol="annotation", analysisCol="analysis") {

  stopifnot(is(scanAnnot, "ScanAnnotationDataFrame"))
  stopifnot(all(hasVariable(scanAnnot, c(consentVar, subjVar, dupVar, omitVar))))
  stopifnot(all(c(annotationCol, analysisCol) %in% names(varMetadata(scanAnnot))))
  
  # consent levels
  conslev <- unique(scanAnnot[[consentVar]])

  # annotation and analysis
  for (type in c("annotation", "analysis")) {
  
    # output files
    annotfile <- paste(dir, "/Sample_", type, ".csv", sep="")
    dupfile <- paste(dir, "/Sample_", type, "_duplicates.csv", sep="")
    ddfile <- paste(dir, "/Sample_", type, "_DD.txt", sep="")
    consent.annotfiles <- paste(dir, "/Sample_", type, "_consent_", conslev, ".csv", sep="")
    consent.dupfiles <- paste(dir, "/Sample_", type, "_duplicates_consent_", conslev, ".csv", sep="")

    # variables to include in sample annotation file
    if (type == "annotation") thistype <- annotationCol else thistype <- analysisCol
    sann <- varLabels(scanAnnot)[varMetadata(scanAnnot)[[thistype]]]

    # remove any samples that will not be posted on dbGaP
    annot2 <- scanAnnot[!scanAnnot[[omitVar]],]

    # split the sample annotation into the main, unduplicated set
    # (one sample per subject) and the duplicated set
    subj <- annot2[annot2[[subjVar]],]
    dups <- annot2[annot2[[dupVar]],]
    stopifnot(nrow(dups) + nrow(subj) == nrow(annot2))
  
    # prepare data dictionary
    meta <- varMetadata(subj)[sann,]
    meta$type <- unlist(lapply(pData(subj)[, sann], class))
    meta$variable <- row.names(meta)
    dd <- meta[,c("variable", "labelDescription", "type")]
    names(dd) <- c("variable", "description", "type")

    # write the files
    subj <- pData(subj)[, sann]
    write.csv(subj, file=annotfile, quote=FALSE, row.names=FALSE, na="")
    dups <- pData(dups)[, sann]
    write.csv(dups, file=dupfile, quote=FALSE, row.names=FALSE, na="")
    write.table(dd, file=ddfile, sep="\t", quote=FALSE, row.names=FALSE, na="")

    # Divide samples by consent group
    if (type == "annotation") {
      for (i in 1:length(conslev)) {
        s <- subj[subj[[consentVar]] == conslev[i],]
        write.csv(s, file=consent.annotfiles[i], quote=FALSE, row.names=FALSE, na="")
      }

      for (i in 1:length(conslev)) {
        s <- dups[dups[[consentVar]] == conslev[i],]
        write.csv(s, file=consent.dupfiles[i], quote=FALSE, row.names=FALSE, na="")
      }
    }
  }
}


# make snp annotation files for dbGaP
# arguments:
#  annot - SnpAnnotationDataFrame
#  dir - directory for output (default current)

dbgapSnpAnnotation <- function(snpAnnot, dir=".",
   annotationCol="annotation", analysisCol="analysis") {

  stopifnot(is(snpAnnot, "SnpAnnotationDataFrame"))
  stopifnot(all(c(annotationCol, analysisCol) %in% names(varMetadata(snpAnnot))))
  
  # annotation and analysis
  for (type in c("annotation", "analysis")) {
    
    # output files
    annotfile <- paste(dir, "/SNP_", type, ".csv", sep="")
    ddfile <- paste(dir, "/SNP_", type, "_DD.txt", sep="")
    
    # variables to include in sample annotation file
    if (type == "annotation") thistype <- annotationCol else thistype <- analysisCol
    sann <- varLabels(snpAnnot)[varMetadata(snpAnnot)[[thistype]]]

    # keep selected variables
    annot2 <- pData(snpAnnot)[, sann]

    # prepare data dictionary
    meta <- varMetadata(snpAnnot)[sann,]
    meta$type <- unlist(lapply(annot2, class))
    meta$variable <- row.names(meta)
    dd <- meta[,c("variable", "labelDescription", "type")]
    names(dd) <- c("variable", "description", "type")

    # write the files
    write.csv(annot2, file=annotfile, quote=FALSE, row.names=FALSE, na="")
    write.table(dd, file=ddfile, sep="\t", quote=FALSE, row.names=FALSE, na="")
  }
}
