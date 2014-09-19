
## function to combine two (or more) gds files.
## TO-DO: add another option - columns to preserve from original snp/scan annotations?
## TO-DO: log file for mismatched SNPs
gdsMerge <- function(genoDataList, sampleList=NULL, snpList=NULL,
                     match.snps.on="position", snpNameList=NULL,
                     sortByScanID=TRUE, newSnpID=TRUE, outPrefix="new") {
  if (is.null(names(genoDataList))) stop("Please supply names for genoDataList")

  stopifnot(all(match.snps.on %in% c("position", "name")))
  
  ## scanIDs from sampleList, or all
  if (!is.null(sampleList)) {
    stopifnot(all(names(sampleList) == names(genoDataList)))
  } else {
    sampleList <- lapply(genoDataList, getScanID)
  }
  scanID.all <- unlist(sampleList, use.names=FALSE)

  ## TO-DO: allow mapping to new scanIDs???
  stopifnot(sum(duplicated(scanID.all)) == 0)

  if (sortByScanID) scanID.all <- sort(scanID.all)

  for (i in 1:length(sampleList)) {
    message(length(sampleList[[i]]), " samples included from ", names(sampleList)[i])
  }

  ## snpIDs from snpList, or all
  if (!is.null(snpList)) {
    stopifnot(all(names(snpList) == names(genoDataList)))
  } else {
    snpList <- lapply(genoDataList, getSnpID)
  }

  ## align snps on chromosome, position, and alleles
  ## TO-DO: currently this gets only the minimum info - also merge other columns?
  ## TO-DO: make log file with list of snps matching on chrom and position but not alleles.
  message("Matching SNPs...")
  for (x in names(snpList)) {
    snpID <- getSnpID(genoDataList[[x]])
    index <- snpID %in% snpList[[x]]
    alleleA <- getAlleleA(genoDataList[[x]], index=index)
    alleleB <- getAlleleB(genoDataList[[x]], index=index)
    snpList[[x]] <- data.frame(snpID=snpID[index],
                      chromosome=getChromosome(genoDataList[[x]], index=index),
                      position=getPosition(genoDataList[[x]], index=index),
                      alleleA=alleleA,
                      alleleB=alleleB,
                      alleles=pasteSorted(alleleA, alleleB),
                      stringsAsFactors=FALSE)
    if (!is.null(snpNameList))
        snpList[[x]][["name"]] <- getSnpVariable(genoDataList[[x]], snpNameList[[x]],
                                                 index=index)
  }
  snp <- snpList[[1]]
  ## add chromosome if matching on position
  if ("position" %in% match.snps.on) match.snps.on <- c("chromosome", match.snps.on)
  ## always match on alleles
  match.snps.on <- c(match.snps.on, "alleles")
  for (x in names(snpList)[-1]) {
    snp <- merge(snp, snpList[[x]], by=match.snps.on,
                 suffixes=c("", paste0(".", x)), sort=FALSE)
  }
  names(snp)[names(snp) == "snpID"] <- paste0("snpID.", names(snpList)[1])
  snp <- snp[order(snp$chromosome, snp$position),]

  if (newSnpID) {
      snp$snpID <- 1:nrow(snp)
  } else {
      snp$snpID <- snp[[paste0("snpID.", names(snpList)[1])]]
  }
  message(nrow(snp), " SNPs matched on ", paste(match.snps.on, collapse=", "))

  ## use first object as allele coding.
  ## make an index for A/B swaps (C/T in one file, T/C in other).
  for (x in names(snpList)[-1]) {
    snp[[paste("swap", x, sep=".")]] <- snp[[paste("alleleA", x, sep=".")]] != snp$alleleA
    message(sum(snp[[paste("swap", x, sep=".")]]), " alleles swapped for ", x)
  }

  ## create new gds file
  message("Creating new GDS file with ", length(scanID.all), " samples and ", nrow(snp), " SNPs")

  gdsfile <- paste0(outPrefix, ".gds")
  gds <- createfn.gds(gdsfile)

  add.gdsn(gds, "sample.id", scanID.all, compress="ZIP.max", closezip=TRUE)
  add.gdsn(gds, "snp.id", snp$snpID, compress="ZIP.max", closezip=TRUE)
  add.gdsn(gds, "snp.chromosome", snp$chromosome, storage="uint8",
           compress="ZIP.max", closezip=TRUE)
  add.gdsn(gds, "snp.position", snp$position, compress="ZIP.max", closezip=TRUE)
  add.gdsn(gds, "snp.allele", paste(snp$alleleA, snp$alleleB, sep="/"),
           compress="ZIP.max", closezip=TRUE)
  ## rsID???
  sync.gds(gds)

  geno.node <- add.gdsn(gds, "genotype", storage="bit2",
                        valdim=c(nrow(snp), length(scanID.all)))
  put.attr.gdsn(geno.node, "snp.order")

  for (i in 1:length(genoDataList)) {
    set <- names(genoDataList)[i]
    message("Reading genotypes from ", set)
    snpID.col <- paste0("snpID.", set)
    snp.index <- match(snp[[snpID.col]], getSnpID(genoDataList[[i]]))
    swap <- snp[[paste0("swap.", set)]]
    for (s in sampleList[[i]]) {
      samp.index <- which(getScanID(genoDataList[[i]]) == s)
      geno <- getGenotype(genoDataList[[i]], scan=c(samp.index, 1), snp=c(1,-1))[snp.index]
      ## if alleles are swapped, new genotype is 2-genotype
      if (i > 1) geno[swap] <- 2 - geno[swap]
      geno[is.na(geno)] <- 3
      write.gdsn(geno.node, geno, start=c(1, which(scanID.all == s)), count=c(-1,1))
    }
  }
  message("Cleaning up")
  closefn.gds(gds)
  cleanup.gds(gdsfile)

  message("Saving annotation")
  cols <- c("snpID", "chromosome", "position", "alleleA", "alleleB", "name",
            names(snp)[grep("snpID.", names(snp), fixed=TRUE)],
            names(snp)[grep("chromosome.", names(snp), fixed=TRUE)],
            names(snp)[grep("position.", names(snp), fixed=TRUE)],
            names(snp)[grep("name.", names(snp), fixed=TRUE)])
  cols <- intersect(cols, names(snp))
  snpAnnot <- SnpAnnotationDataFrame(snp[,cols])
  save(snpAnnot, file=paste0(outPrefix, "_snpAnnot.RData"))

  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=scanID.all,
                                                  stringsAsFactors=FALSE))
  save(scanAnnot, file=paste0(outPrefix, "_scanAnnot.RData"))
}


gdsMergeCheck <- function(genoDataList, outPrefix="new") {
  snpAnnot <- getobj(paste0(outPrefix, "_snpAnnot.RData"))
  genoData <- GenotypeData(GdsGenotypeReader(paste0(outPrefix, ".gds")),
                           snpAnnot=snpAnnot)
  scanID.all <- getScanID(genoData)

  for (i in 1:length(genoDataList)) {
    set <- names(genoDataList)[i]
    message("Reading genotypes from ", set)
    snpID.col <- paste0("snpID.", set)
    snp.index <- match(snpAnnot[[snpID.col]], getSnpID(genoDataList[[i]]))
    scanID.set <- getScanID(genoDataList[[i]])
    for (s in intersect(scanID.set, scanID.all)) {
      geno <- getGenotype(genoDataList[[i]], scan=c(which(scanID.set == s), 1),
                          snp=c(1,-1), char=TRUE, sort=TRUE)[snp.index]
      geno.new <- getGenotype(genoData, scan=c(which(scanID.all == s), 1),
                              snp=c(1,-1), char=TRUE, sort=TRUE)
      stopifnot(allequal(geno, geno.new))
    }
  }
  close(genoData)
}
