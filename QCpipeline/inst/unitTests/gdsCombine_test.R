library(gdsfmt)

.makeTestGds <- function(filename, nsamp=50, nsnp=100, scanBySnp=T, snpIDtype="diff1"){
  
  nchr <- 23
  snp.df <- data.frame(chromosome=rep(1:nchr, each=nsnp), stringsAsFactors=F)
  
  snp.df$position <- sample(1:(nsnp*10), nsnp*nchr, replace=T)
  snp.df <- snp.df[order(snp.df$chromosome, snp.df$position), ]
  snp.df$alleleA <- sample(c("A", "C", "G", "T"), nrow(snp.df), replace=T)
  snp.df$alleleB <- sample(c("A", "C", "G", "T"), nrow(snp.df), replace=T)
  if (snpIDtype %in% "diff1"){
    snp.df$snpID <- 1:nrow(snp.df)
  } else if (snpIDtype %in% "any"){
    snp.df$snpID <- sort(sample(1:(nrow(snp.df)*10), nrow(snp.df)))
  }
  
  
  geno <- matrix(sample(0:3, nsnp*nchr*nsamp, replace=T), nrow=nsamp)
  
  if (file.exists(filename)) stop(paste("filename", filename, "already exists!"))
  gfile <- createfn.gds(filename)
  add.gdsn(gfile, "sample.id", 1:nsamp)
  add.gdsn(gfile, "snp.id", snp.df$snpID)
  add.gdsn(gfile, "snp.chromosome", snp.df$chromosome)
  add.gdsn(gfile, "snp.position", snp.df$position)
  add.gdsn(gfile, "snp.allele", paste(snp.df$alleleA, snp.df$alleleB, sep="/"))
  add.gdsn(gfile, "genotype", geno, storage="bit2")
  
  closefn.gds(gfile)
}

test_gdsCombine <- function(){
  
  for (snptype in c("diff1", "any")){
    
    message(snptype)
    
    fileA <- tempfile()
    .makeTestGds(fileA, snpIDtype=snptype)
    gdsA <- GdsGenotypeReader(fileA)
    
    fileB <- tempfile()
    .makeTestGds(fileB, snpIDtype=snptype)
    gdsB <- GdsGenotypeReader(fileB)
    
    fileC <- tempfile()
    .makeTestGds(fileC, snpIDtype=snptype)
    gdsC <- GdsGenotypeReader(fileC)
    
    gdsList <- list(A=gdsA, B=gdsB, C=gdsC)
    
    filename <- tempfile()
    
    for (b in c(5000, (nsnp(gdsA) + nsnp(gdsB) + nsnp(gdsC) - 1))){
      message(b)
      snp <- gdsCombine(gdsList, filename, blockSize=b, dryRun=TRUE)  
      checkTrue(!file.exists(filename))
      snp <- gdsCombine(gdsList, filename, blockSize=b)  
      
      gds <- GdsGenotypeReader(filename)
      
      geno <- getGenotype(gds)
      
      # check individual datasets A
      for (n in names(gdsList)){
        gds.set <- gdsList[[n]]
        sel <- snp$dataset %in% n
        geno.set <- getGenotype(gds.set)
        checkEquals(geno.set, geno[sel, ])
        checkEquals(getAlleleA(gds)[sel], getAlleleA(gds.set))  
        checkEquals(getAlleleB(gds)[sel], getAlleleB(gds.set))  
        checkEquals(getChromosome(gds)[sel], getChromosome(gds.set))  
        checkEquals(getPosition(gds)[sel], getPosition(gds.set))  
        # check annotation
        checkEquals(snp$alleleA[sel], getAlleleA(gds.set))  
        checkEquals(snp$alleleB[sel], getAlleleB(gds.set))  
        checkEquals(snp$chromosome[sel], getChromosome(gds.set))  
        checkEquals(snp$position[sel], getPosition(gds.set))  
        
      }
      close(gds)
      unlink(filename)
    }
    
    lapply(gdsList, close)
    
    unlink(fileA)
    unlink(fileB)
    unlink(fileC)
    
  }
  
}


test_checkGdsCombine <- function() {
  
  fileA <- tempfile()
  .makeTestGds(fileA)
  gdsA <- GdsGenotypeReader(fileA)
  
  fileB <- tempfile()
  .makeTestGds(fileB)
  gdsB <- GdsGenotypeReader(fileB)
  
  fileC <- tempfile()
  .makeTestGds(fileC)
  gdsC <- GdsGenotypeReader(fileC)
  
  gdsList <- list(A=gdsA, B=gdsB, C=gdsC)
  
  filename <- tempfile()
  filename.fail <- tempfile()
  
  # combine them
  snp <- gdsCombine(gdsList, filename)  
  
  for (blockSize in c(5000, nrow(snp)-1)){
    for (bySnp in c(FALSE, TRUE)){
      # make sure it works
      gds <- GdsGenotypeReader(filename)
      checkGdsCombine(gds, gdsList, snp, bySnp = bySnp, blockSize=blockSize)
      close(gds)
      
      # change a snp chromosome
      file.copy(filename, filename.fail)
      gds <- openfn.gds(filename.fail, readonly=F)
      write.gdsn(index.gdsn(gds, "snp.chromosome"), val=30, start=1, count=1)
      closefn.gds(gds)
      # make sure it fails
      genoData <- GdsGenotypeReader(filename.fail)
      checkException(checkGdsCombine(genoData, gdsList, snp, bySnp = bySnp, blockSize=blockSize))
      close(genoData)
      unlink(filename.fail)
      
      # change a genotype
      file.copy(filename, filename.fail)
      gds <- openfn.gds(filename.fail, readonly=F)
      gGeno <- index.gdsn(gds, "genotype")
      geno <- read.gdsn(gGeno)
      i.na <- which(geno == get.attr.gdsn(gGeno)$missing.value, arr.ind=T)
      # change the first one to non-zero
      write.gdsn(gGeno, val=1, start=(i.na[1, ]), count=c(1,1))
      closefn.gds(gds)
      # make sure it fails
      genoData <- GdsGenotypeReader(filename.fail)
      checkException(checkGdsCombine(genoData, gdsList, snp, bySnp = bySnp, blockSize=blockSize))
      close(genoData)
      unlink(filename.fail)
      
      # change a non-missing value to another non-missing value
      geno.map <- c("0"=1, "1"=2, "2"=0)
      file.copy(filename, filename.fail)
      gds <- openfn.gds(filename.fail, readonly=F)
      gGeno <- index.gdsn(gds, "genotype")
      geno <- read.gdsn(gGeno)
      i.na <- which(geno != get.attr.gdsn(gGeno)$missing.value, arr.ind=T)
      # change the first one to non-zero
      oldval <- geno[i.na[2, "row"], i.na[2, "col"]]
      newval <- geno.map[as.character(oldval)]
      write.gdsn(gGeno, val=newval, start=(i.na[2, ]), count=c(1,1))
      closefn.gds(gds)
      # make sure it fails
      genoData <- GdsGenotypeReader(filename.fail)
      checkException(checkGdsCombine(genoData, gdsList, snp, bySnp = bySnp, blockSize=blockSize))
      close(genoData)
      unlink(filename.fail)
    } 
  }
  
  lapply(gdsList, close)
  
  unlink(fileA)
  unlink(fileB)
  unlink(fileC)
  unlink(filename)
  
}



test_gdsCombine_snpExclude <- function(){
  
  snptype <- 'any'
  
  fileA <- tempfile()
  .makeTestGds(fileA, snpIDtype=snptype)
  gdsA <- GdsGenotypeReader(fileA)
  
  fileB <- tempfile()
  .makeTestGds(fileB, snpIDtype=snptype)
  gdsB <- GdsGenotypeReader(fileB)
  
  fileC <- tempfile()
  .makeTestGds(fileC, snpIDtype=snptype)
  gdsC <- GdsGenotypeReader(fileC)
  
  gdsList <- list(A=gdsA, B=gdsB, C=gdsC)
  
  # snp.exclude
  snp.exclude <- list()
  # do not exclude any from the first file
  snp.exclude[["B"]] <- getSnpID(gdsB)[2]
  snp.exclude[["C"]] <- sample(getSnpID(gdsC), 5)
  
  filename <- tempfile()
  
  snp <- gdsCombine(gdsList, filename, snpExcludeList=snp.exclude)
  
  gds <- GdsGenotypeReader(filename)
  
  geno <- getGenotype(gds)
  
  # check individual datasets A
  for (n in names(gdsList)){
    gds.set <- gdsList[[n]]
    sel <- snp$dataset %in% n
    sel.set <- getSnpID(gds.set) %in% snp$snpID.original[sel]
    geno.set <- getGenotype(gds.set)
    checkEquals(geno[sel, ], geno.set[sel.set, ])
    checkEquals(getAlleleA(gds)[sel], getAlleleA(gds.set)[sel.set])  
    checkEquals(getAlleleB(gds)[sel], getAlleleB(gds.set)[sel.set])  
    checkEquals(getChromosome(gds)[sel], getChromosome(gds.set)[sel.set])  
    checkEquals(getPosition(gds)[sel], getPosition(gds.set)[sel.set])  
    # check annotation
    checkEquals(snp$alleleA[sel], getAlleleA(gds.set)[sel.set])  
    checkEquals(snp$alleleB[sel], getAlleleB(gds.set)[sel.set])  
    checkEquals(snp$chromosome[sel], getChromosome(gds.set)[sel.set])  
    checkEquals(snp$position[sel], getPosition(gds.set)[sel.set])  
    
  }
  close(gds)
  unlink(filename)
  
  lapply(gdsList, close)
  
  unlink(fileA)
  unlink(fileB)
  unlink(fileC)
  
  
}



test_checkGdsCombine_snpExclude <- function() {
  
  fileA <- tempfile()
  .makeTestGds(fileA)
  gdsA <- GdsGenotypeReader(fileA)
  
  fileB <- tempfile()
  .makeTestGds(fileB)
  gdsB <- GdsGenotypeReader(fileB)
  
  fileC <- tempfile()
  .makeTestGds(fileC)
  gdsC <- GdsGenotypeReader(fileC)
  
  gdsList <- list(A=gdsA, B=gdsB, C=gdsC)
  
  filename <- tempfile()
  filename.fail <- tempfile()
  
  # snp exclude
  snp.exclude <- list()
  # do not exclude any from the first file
  snp.exclude[["B"]] <- getSnpID(gdsB)[2]
  snp.exclude[["C"]] <- sample(getSnpID(gdsC), 5)
  
  # combine them
  snp <- gdsCombine(gdsList, filename, snpExcludeList=snp.exclude)  
  
  
  for (bySnp in c(FALSE, TRUE)){
    for (blockSize in c(5000, nrow(snp)-1)){
      
      # make sure it works
      gds <- GdsGenotypeReader(filename)
      checkGdsCombine(gds, gdsList, snp, snpExcludeList=snp.exclude, bySnp=bySnp, blockSize=blockSize)

      # different order for snp exclude list
      checkGdsCombine(gds, gdsList, snp, snpExcludeList=snp.exclude[c("C", "B", "A")], bySnp=bySnp, blockSize=blockSize)
      
      # try running without the snpExcludeList
      checkException(checkGdsCombine(gds, gdsList, snp, bySnp=bySnp, blockSize=blockSize))
      
      # try permuting the snpExcludeList
      snp.exclude.fail <- snp.exclude
      names(snp.exclude.fail) <-c("C", "A")
      checkException(checkGdsCombine(gds, gdsList, snp, snpExcludeList = snp.exclude.fail, bySnp=bySnp, blockSize=blockSize))
      close(gds)
      
      # change a genotype
      file.copy(filename, filename.fail)
      gds <- openfn.gds(filename.fail, readonly=F)
      gGeno <- index.gdsn(gds, "genotype")
      geno <- read.gdsn(gGeno)
      i.na <- which(geno == get.attr.gdsn(gGeno)$missing.value, arr.ind=T)
      # change the first one to non-zero
      write.gdsn(gGeno, val=1, start=(i.na[1, ]), count=c(1,1))
      closefn.gds(gds)
      # make sure it fails
      genoData <- GdsGenotypeReader(filename.fail)
      checkException(checkGdsCombine(genoData, gdsList, snp, snpExcludeList = snp.exclude, bySnp=bySnp, blockSize=blockSize))
      close(genoData)
      unlink(filename.fail)
      
      # change a non-missing value to another non-missing value
      geno.map <- c("0"=1, "1"=2, "2"=0)
      file.copy(filename, filename.fail)
      gds <- openfn.gds(filename.fail, readonly=F)
      gGeno <- index.gdsn(gds, "genotype")
      geno <- read.gdsn(gGeno)
      i.na <- which(geno != get.attr.gdsn(gGeno)$missing.value, arr.ind=T)
      # change the first one to non-zero
      oldval <- geno[i.na[2, "row"], i.na[2, "col"]]
      newval <- geno.map[as.character(oldval)]
      write.gdsn(gGeno, val=newval, start=(i.na[2, ]), count=c(1,1))
      closefn.gds(gds)
      # make sure it fails
      genoData <- GdsGenotypeReader(filename.fail)
      checkException(checkGdsCombine(genoData, gdsList, snp, snpExcludeList = snp.exclude, bySnp=bySnp, blockSize=blockSize))
      close(genoData)
      unlink(filename.fail)
    }
  }
  
  lapply(gdsList, close)
  
  unlink(fileA)
  unlink(fileB)
  unlink(fileC)
  unlink(filename)
  
}
