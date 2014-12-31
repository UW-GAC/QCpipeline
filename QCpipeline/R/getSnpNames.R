readSnpTable <- function(gds, chromosome, variables, verbose=TRUE) {
  
  nodes <- ls.gdsn(gds)
  if (!all(variables %in% nodes)) stop(paste("not all requested variables in gds file. Available variables are:", paste(variables, sep=", ")))
  
  # read chromosome
  if (verbose) message("Reading chromosome info...")
  chrom <- read.gdsn(index.gdsn(gds, "chrom"))
  sel <- chrom %in% paste("chr", chromosome, sep="")
  
  if (sum(sel) == 0) stop(paste("Chromosome", chromosome, "not found in snptab"))
  
  # read required nodes
  snptab <- data.frame(chrom=chrom[sel], stringsAsFactors=F)
  
  if (verbose) message(paste("Reading SNP info for chromosome ", chromosome, "...", sep=""))
  for (variable in setdiff(variables, "chrom")) {
    if (verbose) message(paste(" ", variable))
    snptab[[variable]] <- readex.gdsn(index.gdsn(gds, variable), sel=sel)
  }
  
  snptab
  
}

# assume gds file has been opened
getSnpNames <- function(snpAnnot, gds, chromosome, verbose=TRUE, extraSnptabVars=NULL){
  
  # for now, one chromosome only
  if (length(unique(snpAnnot$chromosome)) > 1) stop("only one chromosome allowed")
  
  
  # check for required snpAnnot names
  required.snp <- c("chromosome", "position", "alleleA", "alleleB")
  stopifnot(all(required.snp %in% varLabels(snpAnnot)))
  snp <- pData(snpAnnot)
  
  # read in the snp table
  # get chromosome to figure out which SNPs to index
  chrom <- read.gdsn(index.gdsn(gds, "chrom"))
  if (XchromCode(snpAnnot) == chromosome) {
    chromCode <- "X"
  } else if (YchromCode(snpAnnot) == chromosome) {
    chromCode <- "Y"
  } else if (MchromCode(snpAnnot) == chromosome) {
    chromCode <- "M"
  } else {
    chromCode <- chromosome
  }
  required <- c("name", "chrom", "chromStart", "chromEnd", "observed", "strand", "class", "exceptions")
  variables <- union(required, extraSnptabVars)
  snptab <- readSnpTable(gds, chromCode, variables, verbose=verbose)
  
  
  if (verbose) message("Selecting SNPs...")
  
  # select only certain rows
  # we only match on class=single, snps that don't have multiple alignments, and snps were chromEnd-chromStart =1 (not sure why this would be bigger for class="single" but it happens)
  sel.class <- snptab$class == "single"
  sel.exception <- !grepl("MultipleAlignment", snptab$exceptions)
  sel.bases <- snptab$chromEnd - snptab$chromStart == 1
  sel <- sel.class & sel.exception & sel.bases
  snptab$snptab.valid <- sel

  snptab.valid <- snptab[sel, ]
  
  # now the input snp annotation to match
  # must be on the same chromosome
  # must have proper alleles
  alleles <- c("A", "C", "G", "T", "0") # the 0 is for study-only monomorphic SNPs
  sel.chromosome <- snp$chromosome %in% chromosome
  sel.alleles <- snp$alleleA %in% alleles & snp$alleleB %in% alleles
  sel <- sel.chromosome & sel.alleles

  snp.valid <- snp[sel, ]
    
  # only merge snp.valid on snpID and position, since we will merge again with the full snp annotation by snpID.
  merged <- merge(snp.valid[, c("snpID", "position")], snptab.valid, by.x="position", by.y="chromEnd", all.x=T)
  merged$snpAnnot.valid <- TRUE

  # do not include position in "merged", since it was already merged with the valid snps by position.
  merged2 <- merge(snp, merged[, !(names(merged) %in% "position")], , by="snpID", all.x=T)
  merged2$snpAnnot.valid[is.na(merged2$snpAnnot.valid)] <- FALSE
  sel <- is.na(merged2$snptab.valid) & merged2$position %in% snptab$chromEnd[!snptab$snptab.valid]
  merged2$snptab.valid[sel] <- FALSE
  sel <- is.na(merged2$snptab.valid) & merged2$position %in% snptab$chromEnd[snptab$snptab.valid]
  merged2$snptab.valid[sel] <- TRUE  

  if (any(duplicated(merged2$snpID))) warning("duplicated snpIDs in merged output! check them and determine which rsID you want to keep.")
  
  list("merged"=merged2, "snptab"=snptab)

}

## conclusions:
# - work on one chromosome at a time
# - from the snp table:
#   - only match class="single"
#   - remove any with MultipleAlignments in "exceptions"
#   - remove any with chromEnd - chromStart != 0
#   - make sure there are no duplicated positions
# - from SNP annotation
#   - remove indels: anything with alleleA and alleleB not part of c("A", "C", "G", "T", "0")
#   - merge based on position
#   - check alleles
#   - make sure that all inconsistent alleles have "InconsistentAllele" in exceptions.




