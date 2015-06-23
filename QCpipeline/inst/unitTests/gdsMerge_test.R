library(gdsfmt)

.testgds <- function(filename, snp.df, scanID, dosage=FALSE) {
    nsnp <- nrow(snp.df)
    nsamp <- length(scanID)
    if (dosage) {
        gdist <- runif(nsnp*nsamp, 0, 2)
        gdist[sample(1:length(gdist), round(0.1*nsnp*nsamp))] <- -1
    } else {
        gdist <- sample(0:3, nsnp*nsamp, replace=TRUE)
    }
    geno <- matrix(gdist, nrow=nsnp, ncol=nsamp)
  
    if (file.exists(filename)) stop(paste("filename", filename, "already exists!"))
    gfile <- createfn.gds(filename)
    add.gdsn(gfile, "sample.id", scanID)
    add.gdsn(gfile, "snp.id", snp.df$snpID)
    add.gdsn(gfile, "snp.rs.id", snp.df$rsID)
    add.gdsn(gfile, "snp.chromosome", snp.df$chromosome)
    add.gdsn(gfile, "snp.position", snp.df$position)
    add.gdsn(gfile, "snp.allele", paste(snp.df$alleleA, snp.df$alleleB, sep="/"))
    geno.node <- add.gdsn(gfile, "genotype", geno, storage=ifelse(dosage, "float32", "bit2"))
    if (dosage) put.attr.gdsn(geno.node, "missing.value", -1)
    closefn.gds(gfile)
}

.testGenoData <- function(filename, snp.df, scanID, dosage=FALSE) {
    .testgds(filename, snp.df, scanID, dosage)
    GenotypeData(GdsGenotypeReader(filename),
                 snpAnnot=SnpAnnotationDataFrame(snp.df))
}

.snp1 <- function() {
    data.frame(snpID=1:20,
               chromosome=1L,
               position=1:20,
               rsID=paste0("rs", 1:20),
               alleleA="A",
               alleleB="G",
               stringsAsFactors=FALSE)
}

## first 5 match on everything
## second 5 match on position but not name
## third 5 match on name but not position
## last 5 match on position and name but not alleles
.snp2 <- function() {
    snp2 <- .snp1()
    snp2$rsID[6:10] <- letters[1:5]
    snp2$position[11:15] <- 111:115
    snp2$alleleA[16:20] <- "C"
    snp2
}

## allele swap for 11-20
.snp3 <- function() {
    snp3 <- .snp1()
    a <- snp3$alleleA[11:20]
    b <- snp3$alleleB[11:20]
    snp3$alleleA[11:20] <- b
    snp3$alleleB[11:20] <- a
    snp3
}

test_gdsMerge <- function() {
    newfile <- tempfile()
    file1 <- tempfile()
    file2 <- tempfile()
    genoDataList <- list(f1=.testGenoData(file1, .snp1(), 1:50),
                         f2=.testGenoData(file2, .snp2(), 51:100))
    
    gdsMerge(genoDataList, outPrefix=newfile, match.snps.on="position", newSnpID=FALSE)
    gdsMergeCheck(genoDataList, outPrefix=newfile)

    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(1:10, getSnpID(new))
    close(new)

    nameList <- list(f1="rsID", f2="rsID")
    gdsMerge(genoDataList, outPrefix=newfile, match.snps.on="name",
             snpNameList=nameList, newSnpID=FALSE)
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(c(1:5,11:15), getSnpID(new))
    close(new)

    gdsMerge(genoDataList, outPrefix=newfile, match.snps.on=c("position", "name"),
             snpNameList=nameList, newSnpID=FALSE)
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(1:5, getSnpID(new))
    close(new)

    lapply(genoDataList, close)
    unlink(c(file1, file2, paste0(newfile, "*")))
}

test_gdsMerge_swap <- function() {
    newfile <- tempfile()
    file1 <- tempfile()
    file2 <- tempfile()
    genoDataList <- list(f1=.testGenoData(file1, .snp1(), 1:50),
                         f2=.testGenoData(file2, .snp3(), 51:100))
    
    gdsMerge(genoDataList, outPrefix=newfile, newSnpID=FALSE)
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(1:20, getSnpID(new))
    close(new)

    lapply(genoDataList, close)
    unlink(c(file1, file2, paste0(newfile, "*")))
}

test_gdsMerge_snpList <- function() {
    newfile <- tempfile()
    file1 <- tempfile()
    file2 <- tempfile()
    genoDataList <- list(f1=.testGenoData(file1, .snp1(), 1:50),
                         f2=.testGenoData(file2, .snp2(), 51:100))

    snpList <- list(f1=1:10, f2=3:20)
    gdsMerge(genoDataList, outPrefix=newfile, snpList=snpList, newSnpID=FALSE)   
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(3:10, getSnpID(new))
    close(new)

    gdsMerge(genoDataList, outPrefix=newfile, snpList=snpList, newSnpID=TRUE)
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(1:8, getSnpID(new))
    close(new)
    
    lapply(genoDataList, close)
    unlink(c(file1, file2, paste0(newfile, "*")))
}

test_gdsMerge_sampleList <- function() {
    newfile <- tempfile()
    file1 <- tempfile()
    file2 <- tempfile()
    genoDataList <- list(f1=.testGenoData(file1, .snp1(), 51:100),
                         f2=.testGenoData(file2, .snp2(), 1:50))

    sampleList <- list(f1=91:100, f2=1:10)
    gdsMerge(genoDataList, outPrefix=newfile, sampleList=sampleList, sortByScanID=FALSE)
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(unlist(sampleList, use.names=FALSE), getScanID(new))
    close(new)

    gdsMerge(genoDataList, outPrefix=newfile, sampleList=sampleList, sortByScanID=TRUE)
    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(sort(unlist(sampleList, use.names=FALSE)), getScanID(new))
    close(new)
    
    lapply(genoDataList, close)
    unlink(c(file1, file2, paste0(newfile, "*")))
}

test_gdsMerge_dosage <- function() {
    newfile <- tempfile()
    file1 <- tempfile()
    file2 <- tempfile()
    genoDataList <- list(f1=.testGenoData(file1, .snp1(), 1:50, dosage=TRUE),
                         f2=.testGenoData(file2, .snp3(), 51:100, dosage=FALSE))
    
    gdsMerge(genoDataList, outPrefix=newfile, newSnpID=FALSE, dosage=TRUE)
    gdsMergeCheck(genoDataList, outPrefix=newfile)

    new <- GdsGenotypeReader(paste0(newfile, ".gds"))
    checkEquals(1:20, getSnpID(new))
    checkEquals("dFloat32", getNodeDescription(new, "genotype")$storage)
    close(new)

    lapply(genoDataList, close)
    unlink(c(file1, file2, paste0(newfile, "*")))
}
