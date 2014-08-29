test_getAnomalySegmentMap <- function(){
  
  nsnp <- 100
  min.pos <- 100
  max.pos <- 1000
  chromosome <- 1
  snp.pos <- sort(sample(min.pos:max.pos, nsnp))
  
  nseg <- 3
  segments <- data.frame(chromosome=chromosome, segment=1:3)
  segments$bp.start <- sort(c(min.pos, sample((min.pos+1):(max.pos-1), 2)))
  segments$bp.end <- segments$bp.start + diff(c(segments$bp.start, max.pos+1))-1
  
  #grSeg <- GRanges(seqnames=paste0("chr", chromosome),
  #                 ranges=IRanges(start=segments$bp.start, end=segments$bp.end))
  
  n.anom <- 6
  anomalies <- data.frame(subjectID=letters[1:n.anom], stringsAsFactors=F)
  anomalies$chromosome <- chromosome
  anomalies$left.base <- NA
  anomalies$right.base <- NA
  # 1) outside the range of snps
  anomalies[1, c("left.base", "right.base")] <- c(0, 1)
  # 2) outside the range of snps to within the first segment
  anomalies[2, c("left.base", "right.base")] <- c(0, sample((segments$bp.start[1]+1):(segments$bp.end[1]-1), 1))
  # 3) within one segment
  p1 <- sample((segments$bp.start[1]+1):(segments$bp.end[1]-1), 1)
  p2 <- sample((segments$bp.start[1]+1):(segments$bp.end[1]-1), 1)
  anomalies[3, c("left.base", "right.base")] <- c(min(p1,p2), max(p1, p2))
  # 4) within two segments
  p1 <- sample((segments$bp.start[1]+1):(segments$bp.end[1]-1), 1)
  p2 <- sample((segments$bp.start[2]+1):(segments$bp.end[2]-1), 1)
  anomalies[4, c("left.base", "right.base")] <- c(p1, p2)
  # 5) spanning all segments
  anomalies[5, c("left.base", "right.base")] <- c(0, max.pos+1)
  # 6) on the boundary between segments -- should zero out only segment 2
  anomalies[6, c("left.base", "right.base")] <- c(segments$bp.start[2], segments$bp.end[2])
  
  
  anom.seg <- QCpipeline:::getAnomSegOverlap(anomalies, segments)
  
  # checks
  checkEquals(nrow(anom.seg), 8)
  # 1) a should not be found
  checkEquals(nrow(anom.seg[anom.seg$subjectID %in% anomalies$subjectID[1], ]), 0)
  # 2) should find 1 only
  checkEquals(anom.seg$segment[anom.seg$subjectID %in% anomalies$subjectID[2]], c(1))
  # 3) should find 1 only
  checkEquals(anom.seg$segment[anom.seg$subjectID %in% anomalies$subjectID[3]], c(1))
  # 4) should find 1,2
  checkEquals(anom.seg$segment[anom.seg$subjectID %in% anomalies$subjectID[4]], c(1,2))
  # 5) should find 1,2,3
  checkEquals(anom.seg$segment[anom.seg$subjectID %in% anomalies$subjectID[5]], c(1,2,3))
  # 6) should find 1,2,3
  checkEquals(anom.seg$segment[anom.seg$subjectID %in% anomalies$subjectID[6]], c(2))
  
}



test_filterImputationSegments <- function(){
  
  probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                          package="GWASdata")
  sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                          package="GWASdata")
  
  # read in prob file to get snp info
  gprobs <- read.table(probfile, as.is=T)
  
  # read in sample file
  samp <- read.table(sampfile, as.is=T, skip=2, header=F)
  names(samp) <- c("subjectID", "altID", "missing")
  
  # snp position
  snp.pos <- as.numeric(gprobs[,3])

  # segment mapping file
  chromosome <- 22
  segments <- data.frame(chromosome=chromosome, segment=1:3, mb.start=c(20.3, 20.4, 20.5), mb.end=c(20.4, 20.5, 20.6), stringsAsFactors=F)
  segments$bp.start <- segments$mb.start*1e6 + 1
  segments$bp.end <- segments$mb.end*1e6
  #grSeg <- GRanges(seqnames=paste0("chr", chromosome),
  #                 ranges=IRanges(start=segments$bp.start, end=segments$bp.end))
  
  # anomalies
  # tests here:
  # 1) outside the range of snps
  # 2) outside the range of snps to within the first segment
  # 3) within one segment
  # 4) within two segments
  # 5) spanning all segments
  # 6) on the boundary of segments
  # 7) random -- split into 2 for the same subject
  n.anom <- 7
  anomalies <- data.frame(subjectID=sample(samp$subjectID, n.anom), stringsAsFactors=F)
  anomalies$chromosome <- chromosome
  anomalies$left.base <- NA
  anomalies$right.base <- NA
  anomalies[1, c("left.base", "right.base")] <- c(0, 1) # outside the range of snps
  anomalies[2, c("left.base", "right.base")] <- c(0, 20.35e6) # outside the range of snps to within the first segment
  anomalies[3, c("left.base", "right.base")] <- c(20.42e6, 20.48e6) # within one segment
  anomalies[4, c("left.base", "right.base")] <- c(20.45e6, 20.55e6) # within two segments
  anomalies[5, c("left.base", "right.base")] <- c(0, 100e6) # spanning all segments
  anomalies[6, c("left.base", "right.base")] <- c(segments$bp.start[2], segments$bp.end[2]) # on the boundary between segments -- should zero out only segment 2
  p1 <- runif(1, min=min(snp.pos), max=max(snp.pos))
  p2 <- runif(1, min=min(snp.pos), max=max(snp.pos))
  anomalies[7, c("left.base", "right.base")] <- c(min(p1,p2), max(p1,p2)) # random -- will end up split into two to test having the same subject twice
  
  # duplicate the last subject
  anomalies[n.anom+1,] <- anomalies[n.anom, ]
  
  # subset the last set of anomalies
  anomalies[n.anom+1, "right.base"] <- ceiling(anomalies[n.anom+1, "left.base"] + (anomalies[n.anom+1, "right.base"] - anomalies[n.anom+1, "left.base"])/3)
  anomalies[n.anom, "left.base"] <- floor(anomalies[n.anom, "right.base"] - (anomalies[n.anom, "right.base"] - anomalies[n.anom, "left.base"])/3)
  
  
  # get overlapping segments
  anom.seg <- QCpipeline:::getAnomSegOverlap(anomalies, segments)
  outfile <- tempfile()
  
  filtVal <- -1
  
  for (b in c(5000, 10, 32)){
    log <- filterImputationSegments(infile=probfile,
                                    outfile=outfile,
                                    anomalies=anomalies, # must have subjectID
                                    sample.annot=samp, # must have subjectID
                                    segments=segments, # segments mapping file
                                    chromosome=chromosome,
                                    filteredValue=filtVal,
                                    block.size=b,
                                    overwrite=TRUE,
                                    verbose=TRUE)
    
    # check them - we already have gprobs
    gprobs.filt <- read.table(outfile, as.is=T)
    
    # check that the snp columns are identical
    checkEquals(dim(gprobs), dim(gprobs.filt))
    checkEquals(gprobs[,1], gprobs.filt[,1])
    checkEquals(gprobs[,2], gprobs.filt[,2])
    checkEquals(gprobs[,3], gprobs.filt[,3])
    checkEquals(gprobs[,4], gprobs.filt[,4])
    checkEquals(gprobs[,5], gprobs.filt[,5])
    
    # check that non-anom samples are ok
    cols <- c("snp", "rsID", "pos", "alleleA", "alleleB", rep(samp$subjectID, each=3))
    i.col <- !(cols %in% anom.seg$subjectID)
    checkEquals(gprobs[, i.col], gprobs.filt[, i.col])
    
    snp.pos <- gprobs[, 3]
    
    # check anomaly by anomaly
    for (subj in unique(anomalies$subjectID)){
      i.seg <- anom.seg$segment[anom.seg$subjectID %in% subj]
      if (length(i.seg) == 0) next
      # get the segments associated with this anomaly
      seg.range <- c(min(segments$bp.start[i.seg]), max(segments$bp.end[i.seg]))
      # get the snps in that segment
      snps.anom <- snp.pos >= seg.range[1] & snp.pos <= seg.range[2]
      i.samp <- cols %in% subj
      checkTrue(all(gprobs.filt[snps.anom, i.samp] == filtVal))
      checkEquals(gprobs[!snps.anom, i.samp], gprobs.filt[!snps.anom, i.samp])
    }
    
    # check log file
    checkEquals(nrow(log), nrow(unique(anom.seg[, c("subjectID", "chromosome", "segment")])))
    checkTrue(setequal(unique(paste(anom.seg$subjectID, anom.seg$chromosome, anom.seg$segment)),
                       paste(log$subjectID, log$chromosome, log$segment)))
  }
  unlink(outfile)
  
}


test_filterImputationSegments_snpBoundary <- function(){
  
  probfile <- system.file("extdata", "imputation", "IMPUTE2", "example.chr22.study.gens",
                          package="GWASdata")
  sampfile <- system.file("extdata", "imputation", "IMPUTE2", "example.study.samples",
                          package="GWASdata")
  
  # read in prob file to get snp info
  gprobs <- read.table(probfile, as.is=T)
  
  # read in sample file
  samp <- read.table(sampfile, as.is=T, skip=2, header=F)
  names(samp) <- c("subjectID", "altID", "missing")
  
  # snp position
  snp.pos <- as.numeric(gprobs[,3])
  
  # segment mapping file
  chromosome <- 22
  boundary <- snp.pos[round(length(snp.pos)/2)]
  
  segments1 <- data.frame(chromosome=chromosome, segment=1:2, bp.start=c(0, boundary), bp.end=c(boundary-1, 100e6), stringsAsFactors=F)
  segments2 <- data.frame(chromosome=chromosome, segment=1:2, bp.start=c(0, boundary+1), bp.end=c(boundary, 100e6), stringsAsFactors=F)
  
  # anomalies
  anomalies <- data.frame(subjectID=sample(samp$subjectID, 1), stringsAsFactors=F)
  anomalies$chromosome <- chromosome
  anomalies$left.base <- boundary
  anomalies$right.base <- boundary
  
  # get overlapping segments
  anom.seg1 <- QCpipeline:::getAnomSegOverlap(anomalies, segments1)
  anom.seg2 <- QCpipeline:::getAnomSegOverlap(anomalies, segments2)
  checkEquals(anom.seg1$segment, 2)
  checkEquals(anom.seg2$segment, 1)
  
  outfile <- tempfile()
  
  filtVal <- -1
  
  j <- which(snp.pos %in% boundary)
  for (b in c(5000, j, j-1, j+1)){

    ## FIRST SET
    log <- filterImputationSegments(infile=probfile,
                                    outfile=outfile,
                                    anomalies=anomalies, # must have subjectID
                                    sample.annot=samp, # must have subjectID
                                    segments=segments1, # segments mapping file
                                    chromosome=chromosome,
                                    filteredValue=filtVal,
                                    block.size=b,
                                    overwrite=TRUE,
                                    verbose=TRUE)
    
    # check them - we already have gprobs
    gprobs.filt <- read.table(outfile, as.is=T)
    
    cols <- c("snp", "rsID", "pos", "alleleA", "alleleB", rep(samp$subjectID, each=3))

    i.samp <- cols %in% anomalies$subjectID
    
    # check that the others are ok
    checkEquals(gprobs[, !i.samp], gprobs.filt[, !i.samp])
    
    # check that only the one segment was filtered
    seg.range <- c(segments1$bp.start[segments1$segment %in% 2], segments1$bp.end[segments1$segment %in% 2])
    snps.anom <- snp.pos >= seg.range[1] & snp.pos <= seg.range[2]
    
    # check that only SNPs in that segment are filtered
    checkEquals(gprobs[!snps.anom, i.samp], gprobs.filt[!snps.anom, i.samp])
    checkTrue(all(gprobs.filt[snps.anom, i.samp] == filtVal))
    
    ## SECOND SET
    log <- filterImputationSegments(infile=probfile,
                                    outfile=outfile,
                                    anomalies=anomalies, # must have subjectID
                                    sample.annot=samp, # must have subjectID
                                    segments=segments2, # segments mapping file
                                    chromosome=chromosome,
                                    filteredValue=filtVal,
                                    block.size=b,
                                    overwrite=TRUE,
                                    verbose=TRUE)
    
    # check them - we already have gprobs
    gprobs.filt <- read.table(outfile, as.is=T)
    
    cols <- c("snp", "rsID", "pos", "alleleA", "alleleB", rep(samp$subjectID, each=3))
    
    i.samp <- cols %in% anomalies$subjectID
    
    # check that the others are ok
    checkEquals(gprobs[, !i.samp], gprobs.filt[, !i.samp])
    
    # check that only the one segment was filtered
    seg.range <- c(segments2$bp.start[segments2$segment %in% 1], segments2$bp.end[segments2$segment %in% 1])
    snps.anom <- snp.pos >= seg.range[1] & snp.pos <= seg.range[2]
    
    # check that only SNPs in that segment are filtered
    checkEquals(gprobs[!snps.anom, i.samp], gprobs.filt[!snps.anom, i.samp])
    checkTrue(all(gprobs.filt[snps.anom, i.samp] == filtVal))
    
  }
  
  unlink(outfile)
  
}

