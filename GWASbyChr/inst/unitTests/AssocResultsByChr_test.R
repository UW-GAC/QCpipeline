test_AssocResultsByChr <- function(){
  
  directory <- file.path(tempdir(), paste(sample(c (letters, LETTERS), 10, replace=T), collapse=""))
  dir.create(directory)
  
  prefix <- "assoc"
  
  simulateAssocResults(directory, prefix)
  
  bcAssoc <- AssocResultsByChr(directory)
  checkEquals(bcAssoc@base, prefix)
  checkEquals(bcAssoc@chromSep, "_chr")
  
  checkEquals(getValidChromosomes(bcAssoc), as.character(1:23))

  # check return of data frame
  chromosome <- sample(getValidChromosomes(bcAssoc), 1)
  assoc <- getAssocResults(bcAssoc, chromosome)
  assoc.chk <- getobj(file.path(directory, paste(prefix, "_chr", chromosome, ".RData", sep="")))
  checkEquals(assoc, assoc.chk)

  # check multiple chromosomes
  chromosome <- sample(getValidChromosomes(bcAssoc), 2)
  assoc <- getAssocResults(bcAssoc, chromosome)
  assoc.chk <- lapply(chromosome, function(x) getobj(file.path(directory, paste(prefix, "_chr", x, ".RData", sep=""))))
  checkEquals(assoc, do.call(rbind, assoc.chk))

  # check all chromosomes
  assoc <- getAssocResults(bcAssoc)
  checkEquals(length(unique(assoc$chromosome)), length(getValidChromosomes(bcAssoc)))

  # check X
  chromosome <- 23
  assoc <- getAssocResults(bcAssoc, chromosome)
  assoc.chk <- getobj(file.path(directory, paste(prefix, "_chr", chromosome, ".RData", sep="")))
  checkEquals(assoc, assoc.chk)
  
  snpID <- sample(assoc$snpID, 1)
  
  ## check looking up snps
  assoc <- lookUpSnps(bcAssoc, chromosome=c(1,2))
  assoc.chk <- getAssocResults(bcAssoc, chromosome=c(1,2))
  checkEquals(assoc, assoc.chk)
  
  snp.ids <- c(1:5, 101:105, snpID)
  assoc <- lookUpSnps(bcAssoc, snps=snp.ids)
  checkEquals(assoc$snpID, snp.ids)
  assoc <- lookUpSnps(bcAssoc, snps=snp.ids, chromosome=1)
  checkEquals(assoc$snpID, 1:5)

  ## add annotation
  annot.dir <- file.path(tempdir(), paste(sample(c (letters, LETTERS), 10, replace=T), collapse=""))
  dir.create(annot.dir)
  prefix <- paste(sample(c(letters, LETTERS), 10, replace=T), collapse="") 
  simulateAnnotation(annot.dir, prefix) 
  bcAnnot <- AnnotationByChr(annot.dir)
  assoc <- lookUpSnps(bcAnnot, bcAssoc, snps=paste0("rs", snp.ids), column="rsID")
  checkEquals(assoc$snpID, snp.ids)
  # add an extra column
  assoc <- lookUpSnps(bcAnnot, bcAssoc, snps=paste0("rs", snp.ids), column="rsID", extraCols="extra")
  checkTrue("extra" %in% names(assoc))
  assoc <- lookUpSnps(bcAnnot, bcAssoc, snps=paste0("rs", snp.ids), column="rsID", chromosome=1)
  checkEquals(assoc$snpID, 1:5)
  
  # check accessor errors
  checkException(getAssocResults(bcAssoc, -1))
  
  
  # check constructor errors
  checkException(AssocResultsByChr(directory, base=substring(prefix, 1, 1)))
  checkException(AssocResultsByChr(directory, chromSep=""))
  
  unlink(c(directory), recursive=T)
}


test_AssocResultsByChr_dupSnpID <- function(){
  
  directory <- file.path(tempdir(), paste(sample(c (letters, LETTERS), 10, replace=T), collapse=""))
  dir.create(directory)
  
  prefix <- "assoc"
  
  simulateAssocResults(directory, prefix)
  
  bcAssoc <- AssocResultsByChr(directory)
  
  # add a duplicated snpID
  chromosome <- 1
  assoc.orig <- getAssocResults(bcAssoc, chromosome=chromosome)
  assoc.dup <- rbind(assoc.orig, assoc.orig[1, ])
  save(assoc.dup, file=file.path(directory, paste0(bcAssoc@base, bcAssoc@chromSep, 1, ".RData")))
  
  assoc.check <- getAssocResults(bcAssoc, chromosome=chromosome)
  checkEquals(assoc.orig, assoc.check)
  
  assoc <- getAssocResults(bcAssoc)
  checkTrue(all(!duplicated(assoc$snpID)))
  
  snpID <- c(assoc.orig$snpID[1], sample(assoc$snpID[!(assoc$snpID %in% assoc.orig$snpID)], 1))
  
  ## check looking up snps
  snp.ids <- c(1:5, 101:105, snpID)
  assoc <- lookUpSnps(bcAssoc, snps=snp.ids)
  checkEquals(assoc$snpID, sort(unique(snp.ids)))
  assoc <- lookUpSnps(bcAssoc, snps=snp.ids, chromosome=1)
  checkEquals(assoc$snpID, 1:5)
  
  unlink(directory, recursive=T)
}
