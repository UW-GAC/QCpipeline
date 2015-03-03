# add data to the file
.addDataLai <- function(gds, snp, geno.list) {
  append.gdsn(index.gdsn(gds, "snp.id"), val=snp$snpID)
  append.gdsn(index.gdsn(gds, "snp.chromosome"), val=snp$chrom)
  #append.gdsn(index.gdsn(gds, "snp.position"), val=snp$pos)
  append.gdsn(index.gdsn(gds, "snp.rs.id"), val=snp$rsID)
  for (a in names(geno.list)){
    vals <- t(geno.list[[a]])
    vals[is.na(vals)] <- 3
    append.gdsn(index.gdsn(gds, paste0("dosage_", a)), val=vals)
  }
}


# create gds file
.createGdsLai <- function(gds.filename, samples, pops, zip="ZIP.max:"){
  
  gds <- createfn.gds(gds.filename)
  add.gdsn(gds, "sample.id", samples, compress=zip, closezip=TRUE)
  add.gdsn(gds, "snp.id", storage="integer", valdim=0, compress=zip)
  add.gdsn(gds, "snp.chromosome", storage="uint8", valdim=0, compress=zip)
  #add.gdsn(gds, "snp.position", storage="integer", valdim=0, compress=zip)
  add.gdsn(gds, "snp.rs.id", storage="string", valdim=0, compress=zip)
  
  for (pop in names(pops)){
    node <- add.gdsn(gds, paste0("dosage_", pop), storage="bit2", valdim=c(length(samples), 0))
    put.attr.gdsn(node, "sample.order")
  }
  
  sync.gds(gds)
  gds
  
}

# return list of dosages of each ancestry
.ancestryDosage <- function(ancestries, pops) {
  
  map <- expand.grid(a1=names(pops), a2=names(pops))
  map$pattern <- paste(pops[map$a1], pops[map$a2], sep="/")
  map <- rbind(map, data.frame(a1="missing", a2="missing", pattern="./."))
  
  dat.list <- list()
  for (n in c(names(pops), "missing")) {
    map$dos <- rowSums(map[, c("a1", "a2")] == n)
    geno.map <- setNames(map$dos, map$pattern)
    
    val <- pops[n]
    dat.list[[n]] <- matrix(geno.map[ancestries], ncol=ncol(ancestries))
    if (!all(dat.list[[n]] %in% 0:2)) stop("ancestries do not equal 0, 1, or 2 for population ", n)
  }
  
  # make sure it all sums to 2
  
  if(!all(Reduce("+", dat.list) == 2)) stop("dosages for all ancestries sum do not sum to 2!")
  
  # set missing values
  missing <- dat.list[["missing"]]
  dat.list <- dat.list[1:length(pops)]
  dat.list <- lapply(dat.list, "[<-", missing == 2, NA)
  
  dat.list
}


.parseLinesRLE <- function(lines, previous, pops) {
  
  map <- expand.grid(a1=names(pops), a2=names(pops))
  map$pattern <- paste(pops[map$a1], pops[map$a2], sep="/")
  map$new.pattern <- pasteSorted(pops[map$a1], pops[map$a2])
  map <- map[map$new.pattern != map$pattern, ]
  
  loc <- str_locate_all(lines, "\t")
  
  # replace 2nd through 8th with ""
  tmp <- lines
  
  loc.start <- sapply(loc, "[", 1)
  loc.end <- sapply(loc, "[", 9)
  
  str_sub(tmp, start=loc.start, end=loc.end) <- " "
  tmp <- str_replace_all(tmp, "\t", " ")
  
  # gsub the patterns -- this is because 0/2 is the same dosage as 2/0
  for (i in 1:nrow(map)){
    tmp <- gsub(map$pattern[i], map$new.pattern[i], tmp)
  }
  
  # need to include the previous row in the rle so we don't run into problems on block edges where duplicated ancestries are included twice
  rle <- Rle(c(previous, tmp))
  
  rle
  
}


laiDosageFile <- function(lafile,
                         outfile,
                         pops,
                         blockSize=5000,
                         onlyUnique=FALSE,
                         dryRun=FALSE){
  
  if (is.null(names(pops))) stop("pops must be a named vector")
  
  # pre-process input file
  la <- file(lafile, "r")

  while (length(s <- readLines(la, n=1)) > 0) {
    if (substr(s, 1, 6) == "#CHROM") {
      header <- scan(text=s, what=character(), quiet=TRUE)
      break
    }
  }
  ncol <- length(header)
  samples <- header[10:ncol]
  
  # create the gds file
  if (!dryRun){
    genoFile <- .createGdsLai(outfile, samples, pops)
  }
  
  # setup for the loop
  last.snpID <- 0
  snp.list <- list()
  counter <- 1
  nlines <- 0
  
  # need a "previous" value before starting the loop, so make it just a random string
  previous <- paste(sample(letters, 10), collapse="")

  while(length(lines <- readLines(la, n = blockSize)) > 0){
    nlines <- nlines + length(lines)

    x <- scan(text=lines, what=character(), nlines=blockSize, quiet=TRUE)
    
    message("Read ", nlines, " lines")
    
    geno <- matrix(x, ncol=ncol, byrow=TRUE)
    chrom <- as.integer(geno[,1, drop=FALSE])
    pos <- as.integer(geno[,2, drop=FALSE])
    rsID <- geno[,3, drop=FALSE]
    geno <- geno[,10:ncol, drop=FALSE]
    
    if (onlyUnique){
      # we want to exclude the rows that are identical to the previous rows
      # to do that, need to paste the local ancestries together to get a string for an RLE
      #tmp <- apply(cbind(chrom, geno), 1, paste, collapse=" ") # this is slow but I don't know if there's a faster way...
      
      # instead of this, do some substringing
      # find the whitespace
      rle <- .parseLinesRLE(lines, previous, pops)
      
      # skip adding anything if all snps are repeat of the previous block's last snp
      # need to update snp annotation though
      if (length(start(rle)) == 1) {
        # THE SAME as previous - need to edit the snp annotation
        n.prev <- nrow(snp.list[[counter-1]])
        snp.list[[counter-1]]$pos.end[n.prev] <- pos[end(rle)[1]-1]
        snp.list[[counter-1]]$n.markers[n.prev] <- snp.list[[counter-1]]$n.markers[n.prev] + runLength(rle)[1]-1
        next
      }
      
      # get start and end values - need to remove the first because that's the previous
      start <- start(rle)[-1] - 1
      end <- end(rle)[-1] - 1
      lengths <- runLength(rle)[-1]
      
      # check if the first row of this block is the same as the last row of the last block
      if (start[1] != 1){
        # THE SAME as previous - need to edit the snp annotation
        n.prev <- nrow(snp.list[[counter-1]])
        snp.list[[counter-1]]$pos.end[n.prev] <- pos[end(rle)[1]-1]
        snp.list[[counter-1]]$n.markers[n.prev] <- snp.list[[counter-1]]$n.markers[n.prev] + runLength(rle)[1]-1
      }
      
      
      # make the snp annotation that will be added
      snp.start <- 1 + last.snpID
      snp.end <- last.snpID + length(start)
      snpID <- snp.start:snp.end
      last.snpID <- snp.end
      
      snp.block <- data.frame(snpID, chromosome=chrom[start],
                              pos.start=pos[start],
                              pos.end=pos[end],
                              rsID=rsID[start],
                              n.markers=lengths,
                              stringsAsFactors=FALSE)
      
      
      # only write the first non-duplicated ancestry
      geno <- geno[start, , drop=FALSE]

      # record last ancestries
      previous <- as.vector(rle[length(rle)])

    } else {
      snp.start <- 1 + last.snpID
      snp.end <- last.snpID + nrow(geno)
      snpID <- snp.start:snp.end
      last.snpID <- snp.end
      
      snp.block <- data.frame(snpID, chromosome=chrom, position=pos, rsID, stringsAsFactors=FALSE)
    }
    
    snp.list[[counter]] <- snp.block
    
    dat.list <- .ancestryDosage(geno, pops)
    
    if (!dryRun){
      .addDataLai(genoFile, snp.block, dat.list)
    }
    
    # end of block updates
    counter <- counter + 1
  }
  ## end of loop
  
  # concatenate snp annotation data frame
  snp.df <- do.call(rbind, snp.list)
  
  # calculate average positions if only unique ariants are being written
  if (onlyUnique){
    snp.df$position <- as.integer(round((snp.df$pos.start + snp.df$pos.end)/2))
  }
  
  
  # clean up gds file -- need to add position info
  if (!dryRun){
    
    # compress nodes that were added as we go
    vars <- ls.gdsn(genoFile)
    vars <- vars[!grepl("^sample", vars)]
    for (v in vars) readmode.gdsn(index.gdsn(genoFile, v))
    
    add.gdsn(genoFile, "snp.position", snp.df$position, compress="ZIP.max", closezip=TRUE)
    
    sync.gds(genoFile)
    
    ## close and cleanup
    filename <- genoFile$filename
    closefn.gds(genoFile)
    cleanup.gds(filename)
    
  }
  
  close(la)
  
  
  # add position here.
  
  snp.df
}



checkLaiDosageFile_working <- function(lafile, genoDataList, pops, snps, blockSize=5000, onlyUnique=FALSE){
  
  if (!setequal(names(genoDataList), names(pops))) stop("names(pops) must equal names(genoDataList)")
  
  # snp checks
  for (g in genoDataList){
    if (!allequal(snps$snpID, getSnpID(g))) stop("snpIDs are not equal")
    if (!allequal(snps$position, getPosition(g))) stop("snpIDs are not equal")
    if (!allequal(snps$chromosome, getChromosome(g))) stop("snpIDs are not equal")
  }
  
  # loop through file
  # pre-process input file
  la <- file(lafile, "r")
  
  while (length(s <- readLines(la, n=1)) > 0) {
    if (substr(s, 1, 6) == "#CHROM") {
      header <- scan(text=s, what=character(), quiet=TRUE)
      break
    }
  }
  ncol <- length(header)
  samples <- header[10:ncol]
  
  # check sampleIDs
  if (!allequal(samples, getScanID(g))) stop("scanID does not match samples in file")
  
  
  # need a "previous" value before starting the loop, so make it just a random string
  previous <- paste(sample(letters, 10), collapse="")
  
  nlines <- 0
  last.snpID <- 0
  # start looping.
  while(length(lines <- readLines(la, n = blockSize)) > 0){
    nlines <- nlines + length(lines)
    
    x <- scan(text=lines, what=character(), nlines=blockSize, quiet=TRUE)
    
    message("Read ", nlines, " lines")
    
    geno <- matrix(x, ncol=ncol, byrow=TRUE)
    chrom <- as.integer(geno[,1, drop=FALSE])
    pos <- as.integer(geno[,2, drop=FALSE])
    rsID <- geno[,3, drop=FALSE]
    geno <- geno[,10:ncol, drop=FALSE]
    
    if (onlyUnique){
      # get rle for these variants
      rle <- .parseLinesRLE(lines, previous, pops)
      
      
      # skip adding anything if all snps are repeat of the previous block's last snp
      # need to update snp annotation though
      if (length(start(rle)) == 1) {
        next
      }
      
      # get start and end values
      start <- start(rle)[2:length(start(rle))] - 1
      end <- end(rle)[2:length(end(rle))] - 1
      
      # only write the first non-duplicated ancestry
      geno <- geno[start, , drop=FALSE]
      
      # save previous ancestries here
      previous <- as.vector(rle[length(rle)])
      
    }
    
    snp.start <- 1 + last.snpID
    snp.end <- last.snpID + nrow(geno)
    snpID <- snp.start:snp.end
    last.snpID <- snp.end
    
    dat <- .ancestryDosage(geno, pops)
    for (n in names(genoDataList)){
      g.pop <- genoDataList[[n]]
      
      geno.chk <- getGenotype(g.pop, snp=c(snp.start, nrow(geno)), scan=c(1,-1))
      if (!allequal(geno.chk, dat[[n]])) stop("genotypes are not equal for population ", n)
    }
    
    
    
  }
  
  close(la)
  message("ok!")
}


# re-written to use GRanges!
checkLaiDosageFile <- function(lafile, genoDataList, pops, snps, blockSize=5000, onlyUnique=FALSE){
 
  if (!setequal(names(genoDataList), names(pops))) stop("names(pops) must equal names(genoDataList)")
  
  # snp checks
  for (g in genoDataList){
    if (!allequal(snps$snpID, getSnpID(g))) stop("snpIDs are not equal")
    if (!allequal(snps$position, getPosition(g))) stop("snpIDs are not equal")
    if (!allequal(snps$chromosome, getChromosome(g))) stop("snpIDs are not equal")
  }
  
  # loop through file
  # pre-process input file
  la <- file(lafile, "r")
  
  while (length(s <- readLines(la, n=1)) > 0) {
    if (substr(s, 1, 6) == "#CHROM") {
      header <- scan(text=s, what=character(), quiet=TRUE)
      break
    }
  }
  ncol <- length(header)
  samples <- header[10:ncol]
  
  # check sampleIDs
  if (!allequal(samples, getScanID(g))) stop("scanID does not match samples in file")

  # make granges object for the snp annotation
  if (onlyUnique){
    gr.snps <- GRanges(seqnames=paste0("chr", snps$chromosome), ranges=IRanges(snps$pos.start, end=snps$pos.end))
  } else {
    gr.snps <- GRanges(seqnames=paste0("chr", snps$chromosome), ranges=IRanges(snps$position, end=snps$position))
  }
  
  # need a "previous" value before starting the loop, so make it just a random string
  previous <- paste(sample(letters, 10), collapse="")
  
  nlines <- 0
  while(length(lines <- readLines(la, n = blockSize)) > 0){
    nlines <- nlines + length(lines)
    
    x <- scan(text=lines, what=character(), nlines=blockSize, quiet=TRUE)
    
    message("Read ", nlines, " lines")
    
    geno <- matrix(x, ncol=ncol, byrow=TRUE)
    chrom <- as.integer(geno[,1, drop=FALSE])
    pos <- as.integer(geno[,2, drop=FALSE])
    rsID <- geno[,3, drop=FALSE]
    geno <- geno[,10:ncol, drop=FALSE]
    
    # iranges for annotation
    gr.block <- GRanges(seqnames=paste0("chr", chrom), ranges=IRanges(start=pos, end=pos))
    
    # find overlaps with annotation
    overlaps <- as.data.frame(findOverlaps(gr.block, gr.snps))

    # get overlaps between annotation and this block
    gr.start <- overlaps$query[!duplicated(overlaps$subjectHits)]
    
    dat <- .ancestryDosage(geno, pops)
    for (n in names(genoDataList)){
      g.pop <- genoDataList[[n]]
      
      minr <- min(snps$snpID[overlaps$subjectHits])
      maxr <- max(snps$snpID[overlaps$subjectHits])
      snp.sel <- c(which(getSnpID(g.pop) %in% minr), maxr-minr+1)
      geno.chk <- getGenotype(g.pop, snp=snp.sel, scan=c(1,-1), drop=FALSE)
      if (!isTRUE(all.equal(geno.chk[overlaps$subjectHits-minr+1, , drop=FALSE], dat[[n]]))) stop("genotypes are not equal for population ", n)
    }
 
  }

  if (onlyUnique){
    stopifnot(nlines == sum(snps$n.markers))
  }
  
  close(la)
  message("ok!")
}
