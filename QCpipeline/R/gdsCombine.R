

gdsCombine <- function(gdsGenoList,filename,
                       blockSize=5000,
                       snpExcludeList=NULL,
                       genotypeDim="scan,snp",
                       verbose=TRUE,
                       dryRun=FALSE){
  
  if (is.null(names(gdsGenoList))) stop("Please supply names for gdsGenoList")
  
  # check that scanIDs are the same -- this is a requirement
  scanID <- getScanID(gdsGenoList[[1]])
  storage <- getNodeDescription(gdsGenoList[[1]], "genotype")$storage
  miss.val <- getAttribute(gdsGenoList[[1]], "missing.value", "genotype")
  for (x in gdsGenoList){
    if (!allequal(scanID, getScanID(x))) stop("scanIDs are not all equal!")
    if (!allequal(storage, getNodeDescription(x, "genotype")$storage)) stop("storage modes ar enot all equal!")
    if (!allequal(order(getSnpID(x)), 1:nsnp(x))) stop("snpIDs are not sorted!")
    if (!allequal(order(getChromosome(x), getPosition(x)), 1:nsnp(x))) stop("snps are not sorted!")
    if (!all.equal(miss.val, getAttribute(x, "missing.value", "genotype"))) stop("missing.values are different!")
  }
  if (is.null(miss.val)) miss.val <- 3 # standard for bit2 files
  
  # check snpExcludeList
  if (!is.null(snpExcludeList)){
    if (!all(names(snpExcludeList) %in% names(gdsGenoList))) stop("names(snpExcludeList) must be in names(gdsGenoList)")
  }
  
  nsamp <- length(scanID)
  
  snpInfoList <- list()
  # extract snp information
  for (n in names(gdsGenoList)) {
    x <- gdsGenoList[[n]]
    snp.df <- data.frame(snpID.original=getSnpID(x),
                         chromosome=getChromosome(x),
                         position=getPosition(x),
                         alleleA=getAlleleA(x),
                         alleleB=getAlleleB(x), stringsAsFactors=F)
    snp.df$dataset <- n
    if (!all(snpExcludeList[[n]] %in% snp.df$snpID.original)){
      stop(paste("snpExcludeList for dataset", n, "are not all in the gds file"))
    }
    snp.df <- snp.df[!(snp.df$snpID.original %in% snpExcludeList[[n]]),]
    snpInfoList[[n]] <- snp.df
    rm(snp.df)
  }
  
  # order snps
  snp <- do.call(rbind, snpInfoList)
  snp <- snp[order(snp$chromosome, snp$position), ]
  snp$snpID <- 1:nrow(snp)
  snp$id <- paste(snp$dataset, snp$snpID.original) # unique id for merging later.. do we still need this? probably
  
  if (!dryRun){
    # create the new file
    # create the new file
    if (file.exists(filename)) stop(paste("filename", filename, "already exists!"))
    gfile <- createfn.gds(filename)
    add.gdsn(gfile, "sample.id", scanID)
    add.gdsn(gfile, "snp.id", snp$snpID)
    add.gdsn(gfile, "snp.chromosome", snp$chromosome)
    add.gdsn(gfile, "snp.position", snp$position)
    add.gdsn(gfile, "snp.allele", paste(snp$alleleA, snp$alleleB, sep="/"))

    if (genotypeDim == "scan,snp") {
        gGeno <- add.gdsn(gfile, "genotype", valdim=c(length(scanID), nrow(snp)), storage=storage)
        # add order for gds file
        put.attr.gdsn(gGeno, "sample.order")
    } else {
        gGeno <- add.gdsn(gfile, "genotype", valdim=c(nrow(snp), length(scanID)), storage=storage)
        put.attr.gdsn(gGeno, "snp.order")
    }
        
    # add missing value
    put.attr.gdsn(gGeno, "missing.value", miss.val)
    
  }
  
  nblocks <- ceiling(nrow(snp)/blockSize)
  
  counter <- 1
  
  if (verbose){
    if (dryRun) message("Looping over SNPs...") else message("Writing genotypes...")
  }
  
  for (block in 1:nblocks){
    snp.start <- counter
    snp.end <- counter + blockSize - 1
    if (snp.end > nrow(snp)) snp.end <- nrow(snp)
    snp.block <- snp[snp.start:snp.end, ]
    nsnp.block <- snp.end - snp.start + 1
    
    if (verbose) message(paste("Block", block, "of", nblocks))
    
    if (genotypeDim == "scan,snp") {
        geno <- matrix(NA, ncol=nsnp.block, nrow=nsamp)
    } else {
        geno <- matrix(NA, ncol=nsamp, nrow=nsnp.block)
    }
    
    # create the combined genotype array for this block
    for (n in names(gdsGenoList)){
      genoData <- gdsGenoList[[n]]
      
      # which SNPs are in this dataset?
      sel.comb <- snp.block$dataset %in% n
      
      # select those SNPs in the gds file
      sel.set <- getSnpID(genoData) %in% snp.block$snpID.original[sel.comb]

      if (genotypeDim == "scan,snp") {
          # read them
          ggeno <- getGenotypeSelection(genoData, snp=sel.set, transpose=TRUE)
          # put them in the combined matrix
          geno[, sel.comb] <- ggeno
      } else {
          # read them
          ggeno <- getGenotypeSelection(genoData, snp=sel.set, transpose=FALSE)
          # put them in the combined matrix
          geno[sel.comb, ] <- ggeno
      }
      
    }
    
    # set NAs to missing value for gds
    geno[is.na(geno)] <- miss.val
    
    # write the genotype to the file
    if (!dryRun){
      if (genotypeDim == "scan,snp") {
          start = c(1, snp.start)
          count <- c(-1, snp.end-snp.start+1)
      } else {
          start = c(snp.start, 1)
          count <- c(snp.end-snp.start+1, -1)
      }
      write.gdsn(gGeno, geno, start=start, count=count)
    }
    
    counter <- snp.end + 1
  }
  
  if (!dryRun){
    closefn.gds(gfile)
  }
  
  return(snp)
  
}



# this checks sample by sample
checkGdsCombine <- function(genoData, gdsGenoList, snp, snpExcludeList=NULL, verbose=TRUE, bySnp=TRUE, blockSize=5000){
  
  if (is.null(names(gdsGenoList))) stop("Please supply names for gdsGenoList")
  
  stopifnot(all(c("snpID", "dataset", "snpID.original", "id") %in% names(snp)))
  
  if (!setequal(names(gdsGenoList), unique(snp$dataset))) stop("names of gdsGenoList do not match snp$dataset!")
  if (!all(names(snpExcludeList %in% names(gdsGenoList)))) stop("names of snpExcludeList are not all in gdsGenoList")
  
  chromosome <- getChromosome(genoData)
  position <- getPosition(genoData)
  alleleA <- getAlleleA(genoData)
  alleleB <- getAlleleB(genoData)
  
  scanID <- getScanID(genoData)
  
  # check annotation data
  for (n in names(gdsGenoList)){
    message(paste("Checking snp and sample info for dataset", n))
    genoData.set <- gdsGenoList[[n]]
    
    sel <- snp$dataset %in% n
    sel.set <- getSnpID(genoData.set) %in% snp$snpID.original[sel]
    
    # snp data
    if (!allequal(getChromosome(genoData.set)[sel.set], chromosome[sel])) stop(paste("chromosome does not match for dataset", n))
    if (!allequal(getPosition(genoData.set)[sel.set], position[sel])) stop(paste("position does not match for dataset", n))
    if (!allequal(getAlleleA(genoData.set)[sel.set], alleleA[sel])) stop(paste("alleleA does not match for dataset", n))
    if (!allequal(getAlleleB(genoData.set)[sel.set], alleleB[sel])) stop(paste("alleleB does not match for dataset", n))
    
    # sample data
    if (!allequal(getScanID(genoData.set), scanID)) stop(paste("scanID does not match for dataset", n))
  }
  message("ok!")
  
  nsamp <- nscan(genoData)
  
  # check excluded snps
  for (n in names(gdsGenoList)){
    snpID.set <- getSnpID(gdsGenoList[[n]])
    snp.set <- snp[snp$dataset %in% n, ]
    if (!setequal(snp.set$snpID.original, setdiff(snpID.set, snpExcludeList[[n]]))) stop(paste("included snpIDs are different for dataset", n))
  }
  
  message("Checking genotype data...")
  if (bySnp){
    
    
    
    nblocks <- ceiling(nrow(snp)/blockSize)
    
    counter <- 1
    
    for (block in 1:nblocks){
      
      snp.start <- counter
      snp.end <- counter + blockSize-1
      if (snp.end > nrow(snp)) snp.end <- nrow(snp)
      snp.block <- snp[snp.start:snp.end, ]
      nsnp.block <- snp.end - snp.start + 1
      
      if (verbose) message(paste("Block", block, "of", nblocks))
      
      geno <- getGenotype(genoData, snp=c(snp.start, nsnp.block), scan=c(1,-1))
      if (class(geno) != "matrix") geno <- t(as.matrix(geno))
      
      for (n in names(gdsGenoList)){
        
        genoDataSet <- gdsGenoList[[n]] # because we will use it later...
        
        # snp selection
        sel.block <- snp.block$dataset %in% n
        
        snp.sel <- getSnpID(genoDataSet) %in% snp.block$snpID.original[sel.block]
        geno.set <- getGenotypeSelection(genoDataSet, snp=snp.sel)
        
        chk <- allequal(geno.set, geno[sel.block, ])
        if (!chk) stop(paste("genotypes not equal for dataset", n, "for SNPs", snp.start, "-", snp.end))
      } 
      counter <- snp.end + 1
      
    }
  } else {
    
    # check snpExclude
    # set up snp-mapping for checking by sample
    # also check snp exclude?
    snp.map <- list()
    snp.set.map <- list()
    for (n in names(gdsGenoList)){
      snpID.set <- getSnpID(gdsGenoList[[n]])
      snpID.set.incl <- setdiff(snpID.set, snpExcludeList[[n]])
      j <- match(paste(n, snpID.set.incl), snp$id)
      snp.map[[n]] <- j
      k <- match(snpID.set.incl, snpID.set)
      snp.set.map[[n]] <- k
    }
    
    # loop over samples
    for (i in 1:length(scanID)){
      
      if (verbose & (i %% 10 == 0)) message(paste("  Sample", i, "of", length(scanID)))
      sid <- scanID[i]
      
      geno <- getGenotype(genoData, scan=c(i,1), snp=c(1,-1))
      
      # check each individual dataset
      for (n in names(gdsGenoList)){
        geno.set <- getGenotype(gdsGenoList[[n]], scan=c(i,1), snp=c(1,-1))
        idx <- snp.map[[n]]
        idx.set <- snp.set.map[[n]]
        if (!allequal(geno.set[idx.set], geno[idx])) stop(paste("genotypes not equal for scanID", sid, "for dataset", n))
      }
      
    }
  }
}
