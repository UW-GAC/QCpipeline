# sets up anomalies granges object
.getAnomalyGRanges <- function(anomalies){
  
  # names
  if (!all(c("subjectID", "chromosome", "left.base", "right.base") %in% names(anomalies))) stop("anomalies data frame must have subjectID, left.base, right.base")
  
  grAnom <- GRanges(seqnames=paste0("chr", anomalies$chromosome),
                    ranges=IRanges(start=anomalies$left.base, end=anomalies$right.base))
  
  grAnom
  
}

# sets up segment granges object
.getSegmentGRanges <- function(segments){
  
  if (!(all(c("segment", "chromosome", "bp.start", "bp.end") %in% names(segments)))) stop("segments must have segment, chromosome", "bp.start, bp.end")
  
  # convert to GRanges
  grSeg <- GRanges(seqnames=paste0("chr", segments$chromosome),
                   ranges=IRanges(start=segments$bp.start, end=segments$bp.end))
  
  grSeg
  
}

# finds oerlap between anomalies and segments
getAnomSegOverlap <- function(anomalies, segments){
  
  grAnom <- .getAnomalyGRanges(anomalies)
  grSeg <- .getSegmentGRanges(segments)
  
  anomalies$id <- 1:nrow(anomalies)
  # map of segments to filter
  overlaps <- findOverlaps(grAnom, grSeg)
  anom.seg <- anomalies[match(queryHits(overlaps), anomalies$id), ]
  anom.seg$segment <- subjectHits(overlaps)
  anomalies$id <- NULL
  
  anom.seg
}


filterImputationSegments <- function(infile,
                           outfile,
                           anomalies, # must have subjectID
                           sample.annot, # must have subjectID, corresponds to gprobs sample file
                           segments, # segments on this chromosome
                           chromosome, # the chromosome to operate on
                           filteredValue=-1,
                           block.size=5000,
                           overwrite=FALSE,
                           verbose=TRUE,
                           checkLineCount=FALSE
                           ){
  
  # write to a tempfile then copy.
  #tmpfile <- tempfile()
  
  # specifies number of header columns in imputation output
  # hard-coded for now
  n.col.skip <- 5
  
  # check that outfile does not exist
  if (file.exists(outfile)) {
    
    if (overwrite){
      if (verbose) {
        message(paste("removing ", outfile))
      }
      file.remove(outfile)

      } else{
      stop(paste("outfile", outfile, "already exists!"))
    }
    
  }
  
  ## check sample annotation
  if (!("subjectID" %in% names(sample.annot))) stop("sample.annot must have subjectID")
  
  ## get overlap between anomalies and segments
  # this does some checks too
  anom.seg <- getAnomSegOverlap(anomalies, segments)
  
  # granges object for segments
  grSeg <- .getSegmentGRanges(segments)
  
  
  # check columns of input
  Cnt <- count.fields(infile)
  col.cnt <- max(Cnt)
  nsnp <- length(Cnt)
  nsamp <- (col.cnt - n.col.skip)/3
  
  # impute 2 checks
  if (nsamp != nrow(sample.annot)) stop("number of samples do not match")
  
  # list for keeping track of filtered samples
  log.list <- list()
  
  # gprobs columns for easier subsetting later
  gprobsCols <- c("snp", "rsID", "pos", "alleleA", "alleleB", rep(sample.annot$subjectID, each=3))
  
  cnt <- 1
  i.log <- 1
  # loop through file to read
  opfile <- file(infile, "r")
  while (length(ss <- scan(opfile, what=character(), quiet=TRUE, nlines=block.size)) > 0){
    
    # keep track of time for rate reporting
    startTime <- Sys.time()
    
    
    dat <- matrix(ss, ncol=col.cnt, byrow=TRUE)
    # make a copy for filtering later -- don't need to convert to numbers since we are just replacing numbers
    dat.filtered <- dat
    
    # get snp position
    grSnp <- GRanges(seqnames=paste0("chr", chromosome),
                     ranges=IRanges(start=as.numeric(dat[,3]), width=1)) # width=1?
    
    # get snp overlaps with segments
    # CHECK HERE what happens for a snp on the border between segments
    overlaps.snp <- findOverlaps(grSnp, grSeg, type="within")
    # CHECK HERE what happens if snp is not found in segment map
    snp.segments <- subjectHits(overlaps.snp)
    
    # loop segment by segment
    for (seg in unique(snp.segments)){
      i.snp <- snp.segments %in% seg
      samp.filt <- anom.seg$subjectID[anom.seg$segment %in% seg]
      if (length(samp.filt) == 0) next
      # columns to set to missing
      i.scan <- gprobsCols %in% samp.filt
      
      dat.filtered[i.snp, i.scan] <- filteredValue
      
      log.list[[i.log]] <- data.frame(subjectID=samp.filt, chromosome=chromosome, segment=seg)
      
      i.log <- i.log + 1
    }
    
    # write dat.filtered
    write.table(dat.filtered, file=outfile, append=TRUE, row.names=F, col.names=F, quote=F)
   
    
    endTime <- Sys.time()
    rate <- format(endTime - startTime, digits=4)
    
    if(verbose) message(paste("Block", ceiling(cnt/block.size), "of", ceiling(nsnp/block.size), "Completed -", rate))
    cnt <- cnt + block.size
    
  }
  
  close(opfile)
  
  if (checkLineCount){
    in.lines <- countLines(infile)
    out.lines <- countLines(outfile)
    if (in.lines != out.lines) stop("line mismatch!")
  }
  
  log <- unique(do.call(rbind, log.list)) # unique to get rid of duplicated rows if looping in blocks
  log
  
}

