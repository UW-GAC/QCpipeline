test_GenotypeDataByChr <- function(){
  
  directory <- file.path(tempdir(), paste(sample(c (letters, LETTERS), 10, replace=T), collapse=""))
  dir.create(directory)
  
  prefix <- paste(sample(c(letters, LETTERS), 10, replace=T), collapse="")
  
  simulateImputedGenotypeData(directory, prefix)
  
  bcData <- GenotypeDataByChr(directory)
  checkEquals(bcData@base, prefix)
  checkEquals(bcData@chromSep, "_chr-")
  
  checkEquals(getValidChromosomes(bcData), as.character(1:23))

  # check snp annotation methods
  chromosome <- sample(getValidChromosomes(bcData), 1)
  snpAnnot <- getSnpAnnotation(bcData, chromosome)
  snpAnnot.chk <- getobj(file.path(directory, paste(prefix, "_chr-", chromosome, "_snpAnnot.RData", sep="")))
  checkEquals(pData(snpAnnot), pData(snpAnnot.chk))

  # check scanID
  gds <- getGenoData(bcData, 1)
  scanID <- getScanID(gds)
  close(gds)
  checkEquals(scanID, getScanID(bcData))

  # check genotypes from snpID
  gds <- getGenoData(bcData, chromosome)
  snpID <- getSnpID(gds)
  snp <- sample(snpID, 1)
  checkTrue(hasSnpID(gds, snp))
  geno <- getGenotype(gds, snp=c(which(snpID == snp), 1), scan=c(1,-1), use.names=TRUE)
  checkEquals(geno, getGenotypeFromSnpID(gds, snp))
  close(gds)
  checkTrue(hasSnpID(bcData, snp))
  checkEquals(geno, getGenotypeFromSnpID(bcData, snp)[1,])

  # multiple snps
  gds <- getGenoData(bcData, chromosome=1)
  snpID <- getSnpID(gds)
  snp1 <- sample(snpID, 2)
  geno <- list()
  for (i in snp1) {
      geno[[as.character(i)]] <- getGenotype(gds, snp=c(which(snpID == i), 1), scan=c(1,-1))
  }
  close(gds)
  
  gds <- getGenoData(bcData, chromosome=2)
  snpID <- getSnpID(gds)
  snp2 <- sample(snpID, 1)
  geno[[as.character(snp2)]] <- getGenotype(gds, snp=c(which(snpID == snp2), 1), scan=c(1,-1))

  geno <- matrix(unlist(geno), nrow=length(geno), byrow=TRUE, dimnames=list(names(geno), getScanID(gds)))
  close(gds)

  snp <- c(snp1, snp2)
  checkTrue(all(hasSnpID(bcData, snp)))
  checkEquals(geno, getGenotypeFromSnpID(bcData, snp))
  
  checkTrue(!hasSnpID(bcData, -1))
  
  # check accessor errors
  checkException(getSnpAnnotation(bcData, -1))
  
  
  # check constructor errors
  checkException(GenotypeDataByChr(directory, base=substring(prefix, 1, 1)))
  checkException(GenotypeDataByChr(directory, chromSep=""))
  
  unlink(directory, recursive=T)
}
