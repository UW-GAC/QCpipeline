simulateImputedGenotypeData <- function(path, prefix, nsamp=50, nsnp=100, nchromosome=23, nsegments=3, missing=NULL, filterColumn="oevar") {
  
  # check that path is a directory
  stopifnot(file_test("-d", path))
  
  # all snps will have sequential snpIDs starting at 1
  snpID.start <- 0
  
  segList = list()
  
  # loop over chromosomes and make one file per chromosome
  for (chromosome in 1:nchromosome){
    
    # make a snp annotation
    snp.df <- data.frame("snpID"=as.integer(1:nsnp + snpID.start), "position"=as.integer(1:nsnp), "chromosome"=as.integer(rep(chromosome, nsnp)))
    snp.df$segment <- ceiling(1:nsnp / nsnp * nsegments)
    snp.df[[filterColumn]] <- sqrt(runif(nsnp))
    
    filename = file.path(path, paste(prefix, "_chr-", chromosome, ".gds", sep=""))
    # simulate genotypes
    simulateGenotypeFile(filename, snp.df, nsamp, missing=missing)
    
    # save the snp annotation
    snpAnnot <- SnpAnnotationDataFrame(snp.df)
    save(snpAnnot, file=file.path(path, paste(prefix, "_chr-", chromosome, "_snpAnnot.RData", sep="")))
    
    sel <- !duplicated(paste(snpAnnot$chromosome, snpAnnot$segment))
    segList[[chromosome]] <- data.frame(chrom=snpAnnot$chromosome[sel], segment=snpAnnot$segment[sel])
    
    snpID.start <- snpID.start + nsnp
  }
  
  segments <- do.call(rbind, segList)
  
  segment.filename <- file.path(path, paste(prefix, "_snp_segment_map.csv", sep=""))
  write.csv(segments, file=segment.filename, quote=FALSE, row.names=FALSE)
  
}


simulateGenotypeFile <- function(filename, snp, nsamp, missing=NULL){

  stopifnot(c("snpID", "chromosome", "position") %in% names(snp))
  
  nsnp <- nrow(snp)
  
  geno <- matrix(NA, nsnp, nsamp)
  
  # get allelic freq by random sampling of a uniform dist on (0,1)
  afreq <- runif(nsnp)
  # simulate genotypes by random sampling 2 gametes for each sample from binomial dist
  for (i in 1:nsnp) geno[i,] <- rbinom(nsamp, 2, afreq[i])
  
  # set missing values
  if (!is.null(missing)) {
    i_miss <- sample(1:length(geno), round(length(geno)*missing))
    geno[i_miss] <- 3
  }
  
  gfile <- createfn.gds(filename)
  add.gdsn(gfile, "sample.id", 1:nsamp)
  add.gdsn(gfile, "snp.id", snp$snpID)
  add.gdsn(gfile, "snp.position", snp$position)
  add.gdsn(gfile, "snp.chromosome", snp$chromosome)
  add.gdsn(gfile, "genotype", geno, storage="bit2")
  closefn.gds(gfile)
  
}


simulateAnnotation <- function(path, prefix, nsnp=100, nchromosome=23) {
  
  # check that path is a directory
  stopifnot(file_test("-d", path))
  
  # all snps will have sequential snpIDs starting at 1
  snpID.start <- 0
  
  # loop over chromosomes and make one file per chromosome
  for (chromosome in 1:nchromosome){
    
    # make a snp annotation
    id <- as.integer(1:nsnp + snpID.start)
    snp.df <- data.frame("snpID"=id,
                         "position"=as.integer(1:nsnp),
                         "chromosome"=as.integer(rep(chromosome, nsnp)),
                         "snpName"=paste0("name", id),
                         "rsID"=paste0("rs",id),
                         "alleleA"=rep("A", nsnp),
                         "alleleB"=rep("G", nsnp),
                         "type"=rep(2, nsnp),
                         "oevar"=runif(nsnp, 0, 1),
                         "extra"=sample(letters, nsnp, replace=TRUE), stringsAsFactors=F)
    
    # save the snp annotation
    snpAnnot <- SnpAnnotationDataFrame(snp.df)
    save(snpAnnot, file=file.path(path, paste(prefix, "_chr-", chromosome, ".RData", sep="")))
    
    snpID.start <- snpID.start + nsnp
  }
}


simulateAssocResults <- function(path, prefix, nsnp=100, nchromosome=23) {
  
  # check that path is a directory
  stopifnot(file_test("-d", path))
  
  # all snps will have sequential snpIDs starting at 1
  snpID.start <- 0
  
  # loop over chromosomes and make one file per chromosome
  for (chromosome in 1:nchromosome){
    
    # make a data frame
    id <- as.integer(1:nsnp + snpID.start)
    assoc <- data.frame("snpID"=id,
                        "position"=as.integer(1:nsnp),
                        "chromosome"=rep(chromToChar(chromosome), nsnp),
                        "alleleA"=rep("A", nsnp),
                        "alleleB"=rep("G", nsnp),
                        "type"=rep(2, nsnp),
                        "pval"=runif(nsnp,0,1),
                        "oevar"=runif(nsnp, 0, 1),
                        "effN"=runif(nsnp, 0, 100),
                        stringsAsFactors=F)
    
    # save 
    save(assoc, file=file.path(path, paste(prefix, "_chr", chromosome, ".RData", sep="")))
    
    snpID.start <- snpID.start + nsnp
  }
}
