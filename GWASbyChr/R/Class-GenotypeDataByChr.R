setClass("GenotypeDataByChr",
         representation(directory = "character",
                        base = "character",
                        chromSep = "character",
                        suffix = "character",
                        annotation = "AnnotationByChr"),
         prototype(directory = "",
                   base = "",
                   chromSep="_chr-",
                   suffix = ".gds"))


# assumes gds
GenotypeDataByChr <- function(directory, base, annotation, chromSep, suffix, ...) {
  if (!file_test("-d", directory)) stop(directory, " is not a directory")
  if (missing(chromSep)) chromSep <- "_chr-"
  if (missing(suffix)) suffix <- ".gds"
  # autodetect base?
  if (missing(base)) {
    files <- list.files(directory, pattern=paste0(suffix, "$"))
    # probably not the best way to do this, but..
    base <- unique(matrix(unlist((strsplit(files, chromSep))), ncol=2, byrow=TRUE)[,1])
    if (length(base) > 1) stop("ambiguous base name in directory. provide base.") # maybe this belongs in validity method
  }
  if (missing(annotation)) {
    annot <- AnnotationByChr(directory=directory, base=base, chromSep=chromSep, suffix="_snpAnnot.RData")
  }
  new("GenotypeDataByChr", directory=directory, base=base, chromSep=chromSep, suffix=suffix, annotation=annot, ...)
}

# to write -- what does it check?
# - if the directory exists.
# - if the directory has gds files
# - if the directory has a snp_segment mapping file
# - if the directory has a scan annotation
# - if the directory has a snp annotation for every file
setValidity("GenotypeDataByChr", 
            function(object) {
              
              # check that the directory exists
              if (!file.exists(object@directory)) return(paste(object@directory, "does not exist"))
              
              # check that the directory is actually a directory
              if (!file_test("-d", object@directory)) return(paste(object@directory, "is not a directory"))
              
              # check that a snp-segment mapping file exists
              #if (!file.exists(file.path(object@directory, paste(object@base, "_snp_segment_map.csv", sep="")))) return("snp segment mapping file does not exist")
              if (!file.exists(getSnpSegmentMap(object))) stop("snp segment mapping file does not exist")
              
              # check that all gds files have a snp annotation
              ## x <- getValidChromosomes(object)
              ## snp.files <- file.path(object@directory, paste(object@base, object@chromSep, x, "_snpAnnot.RData", sep=""))
              ## if (!all(file.exists(snp.files))) stop("not all gds files have associated snp annotations.")
              
              # check scan annotation
              #scan.file <- file.path(object@directory, paste(object@base, "_scanAnnot.Rdata", sep=""))
              #if (!file.exists(scan.file)) stop("scan annotation is missing.")
              TRUE
            })


setMethod("getValidChromosomes",
          signature(object="GenotypeDataByChr"),
          function(object) {
            file.pattern <- paste0("^", object@base, object@chromSep, ".+", object@suffix, "$")
            files <- list.files(object@directory, pattern=file.pattern)
            chroms <- matrix(unlist(strsplit(sub(paste0(object@suffix, "$"), "", files), object@chromSep)),
                             ncol=2, byrow=TRUE)[,2]
            mixedsort(chroms)
          })


setGeneric("getGenotypeFileName", function(object, ...) standardGeneric("getGenotypeFileName"))
setMethod("getGenotypeFileName",
          signature(object="GenotypeDataByChr"),
          function(object, chromosome) {
            tmp <- paste0(object@base, object@chromSep, chromosome, object@suffix)
            file.path(object@directory, tmp)
          })


setGeneric("getGenoData", function(object, ...) standardGeneric("getGenoData"))
setMethod("getGenoData",
          signature(object = "GenotypeDataByChr"),
          function(object, chromosome, snpAnnot=FALSE, scanAnnot=FALSE) {
            # check for valid chromosomes
            if (!(chromosome %in% getValidChromosomes(object))) stop (paste(chromosome, "is not a valid chromosome"))
            
            filename <- getGenotypeFileName(object, chromosome)
            gds <- GdsGenotypeReader(filename)
            
            if (class(snpAnnot) == "logical") {
              if (snpAnnot) snpAnnot <- getSnpAnnotation(object, chromosome) else snpAnnot <- NULL
            } else if (class(snpAnnot) != "SnpAnnotationDataFrame") snpAnnot <- NULL
            
            if (class(scanAnnot) == "logical") {
              if (scanAnnot) scanAnnot <- getScanAnnotation(object) else scanAnnot <- NULL
            } else if (class(scanAnnot) != "ScanAnnotationDataFrame") scanAnnot <- NULL
            
            
            GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
          })



# returns the scan annotation for that chromosome
# generic now defined in GWASTools
#setGeneric("getScanAnnotation", function(object, ...) standardGeneric("getScanAnnotation"))
setMethod("getScanAnnotation",
          signature(object = "GenotypeDataByChr"),
          function(object) {
            getobj(file.path(object@directory, paste0(object@base, "_scanAnnot.RData")))  
          })


# returns the snp annotation for the specified chromosome
#setGeneric("getSnpAnnotation", function(object, ...) standardGeneric("getSnpAnnotation"))
setMethod("getSnpAnnotation",
          signature(object = "GenotypeDataByChr"),
          function(object, chromosome) {
            getSnpAnnotation(object@annotation, chromosome)
          })


setGeneric("getSnpSegmentMap", function(object, ...) standardGeneric("getSnpSegmentMap"))
setMethod("getSnpSegmentMap",
          signature(object="GenotypeDataByChr"),
          function(object) {
            file.path(object@directory, paste0(object@base, "_snp_segment_map.csv"))
          })

# getScanID method, which uses the first gds file
setMethod("getScanID",
          signature(object="GenotypeDataByChr"),
          function(object) {
            gds <- getGenoData(object, chromosome=getValidChromosomes(object)[1])
            scanID <- getScanID(gds)
            close(gds)
            scanID
          })

# hasSnpID, to check if a snpID is in olgaData
setGeneric("hasSnpID", function(object, ...) standardGeneric("hasSnpID"))
setMethod("hasSnpID",
          signature(object="GdsGenotypeReader"),
          function(object, snpID, ...) {
            snpID %in% getSnpID(object)
          })
setMethod("hasSnpID",
          signature(object="GenotypeData"),
          function(object, snpID, ...) {
            hasSnpID(object@data, snpID, ...)
          })
setMethod("hasSnpID",
          signature(object="GenotypeDataByChr"),
          function(object, snpID, ...) {
            found <- rep(FALSE, length(snpID))
            for (chromosome in getValidChromosomes(object)) {
              gds <- getGenoData(object, chromosome=chromosome)
              found <- found | hasSnpID(gds, snpID, ...)
              close(gds)
              if (all(found)) break
            }
            found
          })

# getGenotypeFromSnpID method, which finds the right gds file to open first
setGeneric("getGenotypeFromSnpID", function(object, ...) standardGeneric("getGenotypeFromSnpID"))
setMethod("getGenotypeFromSnpID",
          signature(object="GdsGenotypeReader"),
          function(object, snpID, ...) {
            geno <- NULL           
            snp.ids <- getSnpID(object)
            ind <- na.omit(match(snpID, snp.ids))
            if (length(ind) > 0) {
              geno <- getGenotypeSelection(object, snp=ind, ...)
            }
            geno
          })
setMethod("getGenotypeFromSnpID",
          signature(object="GenotypeData"),
          function(object, snpID, ...) {
            getGenotypeFromSnpID(object@data, snpID, ...)
          })
setMethod("getGenotypeFromSnpID",
          signature(object="GenotypeDataByChr"),
          function(object, snpID, order=c("selection", "file"), transpose=FALSE, ...) {
            order <- match.arg(order)
            geno <- NULL
            for (chromosome in getValidChromosomes(object)) {
              gds <- getGenoData(object, chromosome=chromosome)
              snpID.chr <- intersect(snpID, getSnpID(gds))
              geno.chr <- getGenotypeFromSnpID(gds, snpID.chr, drop=FALSE, use.names=TRUE,
                                               order=order, transpose=transpose, ...)
              close(gds)
              if (transpose) {
                  geno <- cbind(geno, geno.chr)
                  if (all(paste0("snp", snpID) %in% colnames(geno))) break
              } else {
                  geno <- rbind(geno, geno.chr)
                  if (all(paste0("snp", snpID) %in% rownames(geno))) break
              }
            }
            geno
          })

# see what else is needed.

# write some unit tests.. somehow.
