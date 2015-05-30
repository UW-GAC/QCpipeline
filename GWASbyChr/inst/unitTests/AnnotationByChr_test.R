test_AnnotationByChr <- function(){
  
  directory <- file.path(tempdir(), paste(sample(c (letters, LETTERS), 10, replace=T), collapse=""))
  dir.create(directory)
  
  (prefix <- paste(sample(c(letters, LETTERS), 10, replace=T), collapse=""))
  
  simulateAnnotation(directory, prefix)
  
  bcAnnot <- AnnotationByChr(directory)
  checkEquals(bcAnnot@base, prefix)
  checkEquals(bcAnnot@chromSep, "_chr-")
  
  checkEquals(getValidChromosomes(bcAnnot), as.character(1:23))

  # check snp annotation methods
  chromosome <- sample(getValidChromosomes(bcAnnot), 1)
  snpAnnot <- getSnpAnnotation(bcAnnot, chromosome)
  snpAnnot.chk <- getobj(file.path(directory, paste(prefix, "_chr-", chromosome, ".RData", sep="")))
  checkEquals(snpAnnot, snpAnnot.chk)

  # read in multiple chromosomes
  chromosome <- sample(getValidChromosomes(bcAnnot), 3)
  chromosome <- as.numeric(chromosome)
  snpAnnot <- getSnpAnnotation(bcAnnot, chromosome)
  snpAnnot.1 <- getobj(file.path(directory, paste(prefix, "_chr-", sort(chromosome)[1], ".RData", sep="")))
  snpAnnot.2 <- getobj(file.path(directory, paste(prefix, "_chr-", sort(chromosome)[2], ".RData", sep="")))
  snpAnnot.3 <- getobj(file.path(directory, paste(prefix, "_chr-", sort(chromosome)[3], ".RData", sep="")))  
  snp.chk <- rbind(pData(snpAnnot.1), pData(snpAnnot.2), pData(snpAnnot.3))
  checkEquals(pData(snpAnnot), snp.chk)
  
  
  # look up SNPs
  annot.chk <- getSnpAnnotation(bcAnnot, chromosome=1)
  annot1 <- lookUpSnps(bcAnnot, chromosome=1)
  checkEquals(annot.chk[, varLabels(annot1)], annot1)
  
  annot.chk <- getSnpAnnotation(bcAnnot, chromosome=c(2,3))
  annot1 <- lookUpSnps(bcAnnot, chromosome=c(2,3))
  checkEquals(annot.chk[, varLabels(annot1)], annot1)
  
  snp.ids <- c(1:5, 101:105)
  annot1 <- lookUpSnps(bcAnnot, snps=snp.ids)
  checkEquals(annot1$snpID, snp.ids)
  annot2 <- lookUpSnps(bcAnnot, snps=paste0("rs", snp.ids), column="rsID")
  checkIdentical(annot1, annot2)
  annot1 <- lookUpSnps(bcAnnot, snps=snp.ids, chromosome=1)
  checkEquals(annot1$snpID, 1:5)
  
  snp.ids <- snpAnnot$snpID[1:5]
  annot1 <- lookUpSnps(bcAnnot, snps=snp.ids, extraCols="extra")
  head(pData(annot1))
  checkEquals(annot1$extra[annot1$snpID %in% snp.ids], snpAnnot$extra[snpAnnot$snpID %in% snp.ids])

  # wrong column
  annot1 <- lookUpSnps(bcAnnot, snps=snp.ids, column="rsID")
  checkTrue(nrow(annot1) == 0)
  
  # check accessor errors
  checkException(getSnpAnnotation(bcAnnot, -1))
  
  
  # check constructor errors
  checkException(AnnotationByChr(directory, base=substr(prefix, 1, 1)))
  checkException(AnnotationByChr(directory, chromSep=""))
  
  unlink(directory, recursive=T)
}
