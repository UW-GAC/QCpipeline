library(gdsfmt)
library(S4Vectors)

.writeTestLaiFile <- function(outfile) {
  
  # make the data
  dat <- data.frame(CHROM=rep(c(1,2), each=10), POS=rep(1:10, times=2), ID=paste0("rs", 1:20), REF="-", ALT="-", QUAL="-", FILTER="-", INFO="-", FORMAT="-")
  dat$zz100 <- c(rep("2/0", 5), rep("0/2", 11), rep("1/1", 4))
  dat$zz101 <- c(rep("2/0", 19), "1/1")
  dat$zz102 <- "1/0"
  names(dat) <- gsub("^zz", "", names(dat))
  
  
  gz1 <- gzfile(outfile, "w")
  write.table(dat, file=gz1, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, na="")
  close(gz1)
  
  # read it back in, add header rows, and write back out
  gz <- file(outfile, "r")
  lines <- readLines(gz)
  lines <- c("##file format: similar to VCFv4.1",
             "##fileDate=2014-11-28",
             "##description: 0=European ancestry, 1=African ancestry, 2=Native American ancestry",
             "#CHROM  POS  ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	100	101	102", lines)
  close(gz)
  
  gz1 <- gzfile(outfile, "w")
  writeLines(lines, gz1)
  close(gz1)
  
}

test_parseLinesRLE <- function(){
  
  pops <- c("eur"=0, "afr"=1, "amer"=2)
  
  lafile <- tempfile()
  .writeTestLaiFile(lafile)
  
  lines <- readLines(lafile)
  lines <- lines[5:length(lines)]
  
  previous <- paste(sample(letters, 10), collapse="")
  
  dat <- read.table(lafile, colClasses="character", as.is=TRUE)
  dat <- dat[, c(1, 10:12)]
  m <- expand.grid(a1=pops, a2=pops)
  m$pattern <- paste(m$a1, m$a2, sep="/")
  m$new.pattern <- pasteSorted(m$a1, m$a2, sep="/")
  for (i in 1:nrow(m)) {
    dat[, 2] <- gsub(m$pattern[i], m$new.pattern[i], dat[, 2])
    dat[, 3] <- gsub(m$pattern[i], m$new.pattern[i], dat[, 3])
    dat[, 4] <- gsub(m$pattern[i], m$new.pattern[i], dat[, 4])
  }
  chk <- Rle(c(previous, apply(dat, 1, paste, collapse=" ")))
  
  rle <- QCpipeline:::.parseLinesRLE(lines, previous, pops)
  checkIdentical(chk, rle)

  # check previous vector repeated
  previous <- as.vector(rle)[2]
  chk <- Rle(c(previous, apply(dat, 1, paste, collapse=" ")))
  rle <- QCpipeline:::.parseLinesRLE(lines, previous, pops)
  checkIdentical(chk, rle)
  
}


test_ancestryDosage <- function(){
  
  pops <- c("eur"=0, "afr"=1, "amer"=2)
  
  # a known sample
  geno <- matrix("0/0")
  chk <- QCpipeline:::.ancestryDosage(geno, pops)
  checkEquals(chk$eur, matrix(2))
  checkEquals(chk$afr, matrix(0))
  checkEquals(chk$amer, matrix(0))
  
  chk <- QCpipeline:::.ancestryDosage(geno, c("afr"=0, "eur"=1, "amer"=2))
  checkEquals(chk$eur, matrix(0))
  checkEquals(chk$afr, matrix(2))
  checkEquals(chk$amer, matrix(0))
  
  
  geno <- matrix(c("0/0", "0/1", "0/2", "1/0", "1/1", "1/2", "2/0", "2/1", "2/2"), ncol=3)
  chk <- QCpipeline:::.ancestryDosage(geno, pops)
  checkEquals(names(chk), names(pops))
  checkEquals(chk$eur, matrix(c(2,1,1,1,0,0,1,0,0), ncol=3))
  checkEquals(chk$afr, matrix(c(0,1,0,1,2,1,0,1,0), ncol=3))
  checkEquals(chk$amer, matrix(c(0,0,1,0,0,1,1,1,2), ncol=3))
  
  geno <- matrix(c("NA/NA", "0/1", "0/2"))
  checkException(QCpipeline:::.ancestryDosage(geno, pops))
  
  # check other exceptions
  checkException(QCpipeline:::.ancestryDosage(matrix(c("11")), pops))
  checkException(QCpipeline:::.ancestryDosage(matrix(c("3/3")), pops))
  
  # a random sample
  nsamp <- 10
  ndos <- 20
  
  allele1 <- matrix(sample(names(pops), nsamp*ndos, replace=TRUE), ncol=nsamp, nrow=ndos)
  allele2 <- matrix(sample(names(pops), nsamp*ndos, replace=TRUE), ncol=nsamp, nrow=ndos)
  geno <- matrix(paste(pops[allele1], pops[allele2], sep="/"), ncol=nsamp, nrow=ndos)
  chk <- QCpipeline:::.ancestryDosage(geno, pops)
  checkEquals(chk$eur, (allele1 == "eur") + (allele2 == "eur"))
  checkEquals(chk$afr, (allele1 == "afr") + (allele2 == "afr"))
  checkEquals(chk$amer, (allele1 == "amer") + (allele2 == "amer"))
  
}


test_laiDosageFile <- function(){
  
  lafile <- tempfile()
  gdsfile <- tempfile()
  
  .writeTestLaiFile(lafile)
  
  # read in test file
  dat <- read.table(lafile, stringsAsFactors=FALSE)
  names(dat)[1] <- "chromosome"
  names(dat)[2] <- "position"
  names(dat)[3] <- "rsID"
  ancestries <- as.matrix(dat[, 10:12])
  
  pops <- c("eur"=0, "afr"=1, "amer"=2)
  
  for (b in c(1000, 4, 3, 1)) {
    # make the file
    snps <- laiDosageFile(lafile, gdsfile, pops, onlyUnique=FALSE, blockSize=b)
    
    checkEquals(snps$snpID, 1:nrow(dat))
    checkEquals(snps$rsID, dat$rsID)
    checkEquals(snps$position, dat$position)
    
    gds <- openfn.gds(gdsfile)
    gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
    gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
    gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")
    genoData.eur <- GenotypeData(gds.eur, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.afr <- GenotypeData(gds.afr, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.amer <- GenotypeData(gds.amer, snpAnnot=SnpAnnotationDataFrame(snps))
    
    # check dosages
    dat.chk <- QCpipeline:::.ancestryDosage(ancestries, pops)
    checkEquals(dat.chk[["eur"]], getGenotype(genoData.eur))
    checkEquals(dat.chk[["afr"]], getGenotype(genoData.afr))
    checkEquals(dat.chk[["amer"]], getGenotype(genoData.amer))
    
    # check scanID
    checkEquals(getScanID(genoData.eur), c("100", "101", "102"))
    checkEquals(getScanID(genoData.afr), c("100", "101", "102"))
    checkEquals(getScanID(genoData.amer), c("100", "101", "102"))
    
    genoDataList <- list("eur"=genoData.eur, "afr"=genoData.afr, "amer"=genoData.amer)
    checkLaiDosageFile(lafile, genoDataList, pops, snps=snps, onlyUnique=FALSE, blockSize=b)
    
    close(gds.eur)
    unlink(gdsfile)
  }
  
}


test_laiDosageFile_onlyUnique <- function(){
  
  lafile <- tempfile()
  gdsfile <- tempfile()
  
  .writeTestLaiFile(lafile)
  
  # read in test file
  dat <- read.table(lafile, stringsAsFactors=FALSE)
  names(dat)[1] <- "chromosome"
  names(dat)[2] <- "position"
  names(dat)[3] <- "rsID"
  ancestries <- as.matrix(dat[, 10:12])
  
  keep.row <- c(1, 11, 17, 20)
  pops <- c("eur"=0, "afr"=1, "amer"=2)
  
  for (b in c(1000, 4, 3, 1)) {
    # make the file
    snps <- laiDosageFile(lafile, gdsfile, pops, onlyUnique=TRUE, blockSize=b)
    
    checkEquals(snps$snpID, 1:4)
    checkEquals(snps$rsID, dat$rsID[keep.row])
    checkEquals(snps$pos.start, dat$position[keep.row])
    checkEquals(snps$pos.end, dat$position[c((keep.row-1)[2:length(keep.row)], keep.row[length(keep.row)])])
    checkEquals(snps$position, round((snps$pos.start + snps$pos.end)/2))
    
    gds <- openfn.gds(gdsfile)
    gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
    gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
    gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")
    genoData.eur <- GenotypeData(gds.eur, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.afr <- GenotypeData(gds.afr, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.amer <- GenotypeData(gds.amer, snpAnnot=SnpAnnotationDataFrame(snps))
    
    # check dosages
    dat.chk <- QCpipeline:::.ancestryDosage(ancestries, pops)
    dat.chk <- lapply(dat.chk, '[', keep.row, )
    checkEquals(dat.chk[["eur"]], getGenotype(genoData.eur))
    checkEquals(dat.chk[["afr"]], getGenotype(genoData.afr))
    checkEquals(dat.chk[["amer"]], getGenotype(genoData.amer))
    
    # check scanID
    checkEquals(getScanID(genoData.eur), c("100", "101", "102"))
    checkEquals(getScanID(genoData.afr), c("100", "101", "102"))
    checkEquals(getScanID(genoData.amer), c("100", "101", "102"))
    
    genoDataList <- list("eur"=genoData.eur, "afr"=genoData.afr, "amer"=genoData.amer)
    checkLaiDosageFile(lafile, genoDataList, pops, snps=snps, onlyUnique=TRUE, blockSize=b)
    
    close(gds.eur)
    unlink(gdsfile)
  }
  
}


test_checkLaiDosageFile <- function() {
  
  lafile <- tempfile()
  gdsfile <- tempfile()
  
  .writeTestLaiFile(lafile)
  pops <- c("eur"=0, "afr"=1, "amer"=2)
  
  for (b in c(1000, 4, 2, 1)){
    snps <- laiDosageFile(lafile, gdsfile, pops, onlyUnique=FALSE, blockSize=b)
    
    gds <- openfn.gds(gdsfile)
    gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
    gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
    gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")
    genoData.eur <- GenotypeData(gds.eur, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.afr <- GenotypeData(gds.afr, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.amer <- GenotypeData(gds.amer, snpAnnot=SnpAnnotationDataFrame(snps))
    
    genoDataList <- list("eur"=genoData.eur, "afr"=genoData.afr, "amer"=genoData.amer)
    # make sure it works
    checkLaiDosageFile(lafile, genoDataList, pops, snps, blockSize=b)  
    
    # edit the file
    close(genoDataList[[1]])
    gds <- openfn.gds(gdsfile, readonly=FALSE)
    node <- index.gdsn(gds, "dosage_eur")
    val <- read.gdsn(node, start=c(1,1), count=c(1,1))
    write.gdsn(node, abs(val-1), start=c(1,1), count=c(1,1))
    closefn.gds(gds)
    # now it should crash
    gds <- openfn.gds(gdsfile)
    gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
    gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
    gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")
    genoData.eur <- GenotypeData(gds.eur, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.afr <- GenotypeData(gds.afr, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.amer <- GenotypeData(gds.amer, snpAnnot=SnpAnnotationDataFrame(snps))
    
    genoDataList <- list("eur"=genoData.eur, "afr"=genoData.afr, "amer"=genoData.amer)
    checkException(checkLaiDosageFile(lafile, genoDataList, pops, snps, blockSize=b))  
    
    close(genoDataList[[1]])  
    
    unlink(gdsfile)
  }
  unlink(lafile)
  
}



test_checkLaiDosageFile_onlyUnique <- function() {
  
  lafile <- tempfile()
  gdsfile <- tempfile()
  
  .writeTestLaiFile(lafile)
  pops <- c("eur"=0, "afr"=1, "amer"=2)
  
  for (b in c(1000, 4, 2, 1)){
    snps <- laiDosageFile(lafile, gdsfile, pops, onlyUnique=TRUE, blockSize=b)
    
    gds <- openfn.gds(gdsfile)
    gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
    gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
    gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")
    genoData.eur <- GenotypeData(gds.eur, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.afr <- GenotypeData(gds.afr, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.amer <- GenotypeData(gds.amer, snpAnnot=SnpAnnotationDataFrame(snps))
    
    genoDataList <- list("eur"=genoData.eur, "afr"=genoData.afr, "amer"=genoData.amer)
    # make sure it works
    checkLaiDosageFile(lafile, genoDataList, pops, snps, onlyUnique=TRUE, blockSize=b)
    
    # make it crash with a correct file but incorrect annotation
    snps.fail <- snps
    snps.fail$pos.end[2] <- snps.fail$pos.start[3]
    snps.fail$pos.start[3] <- snps$pos.start[3]+1
    checkException(checkLaiDosageFile(lafile, genoDataList, pops, snps.fail, onlyUnique=TRUE, blockSize=b))
    
    # edit the file
    close(genoDataList[[1]])
    gds <- openfn.gds(gdsfile, readonly=FALSE)
    node <- index.gdsn(gds, "dosage_eur")
    val <- read.gdsn(node, start=c(1,1), count=c(1,1))
    write.gdsn(node, abs(val-1), start=c(1,1), count=c(1,1))
    closefn.gds(gds)
    # now it should crash
    gds <- openfn.gds(gdsfile)
    gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
    gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
    gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")
    genoData.eur <- GenotypeData(gds.eur, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.afr <- GenotypeData(gds.afr, snpAnnot=SnpAnnotationDataFrame(snps))
    genoData.amer <- GenotypeData(gds.amer, snpAnnot=SnpAnnotationDataFrame(snps))
    
    genoDataList <- list("eur"=genoData.eur, "afr"=genoData.afr, "amer"=genoData.amer)
    
    checkException(checkLaiDosageFile(lafile, genoDataList, pops, snps, onlyUnique=TRUE, blockSize=b))
    close(genoDataList[[1]])  
    
    unlink(gdsfile)
  }
  unlink(lafile)
  
}