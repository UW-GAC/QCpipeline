setClass("AssocResultsByChr",
         representation(directory = "character",
                        base = "character",
                        chromSep = "character"),
         prototype(directory = "",
                   base = "",
                   chromSep="_chr"))


AssocResultsByChr <- function(directory, base, chromSep, ...) {
  # autodetect base?
  if (missing(chromSep)) chromSep <- "_chr"
  if (missing(base)) {
    files <- list.files(directory, pattern="assoc.+RData")
    if (length(files) == 0) files <- list.files(directory, pattern="meta.+RData")
    if (length(files) == 0) files <- list.files(directory, pattern="RData$")
    # probably not the best way to do this, but..
    base <- unique(matrix(unlist((strsplit(files, chromSep))), ncol=2, byrow=TRUE)[,1])
    if (length(base) > 1) stop("ambiguous base name in directory. provide base.") # maybe this belongs in validity method
  }
  new("AssocResultsByChr", directory=directory, base=base, chromSep=chromSep, ...)
}

# to write -- what does it check?
# - if the directory exists.
# - if the directory has RData files
setValidity("AssocResultsByChr", 
            function(object) {
              
              # check that the directory exists
              if (!file.exists(object@directory)) return(paste(object@directory, " does not exist"))
              
              # check that the directory is actually a directory
              if (!file_test("-d", object@directory)) return(paste(object@directory, " does not exist"))

              # check that there is at least one RData file
              file.pattern <- paste0("^", object@base, object@chromSep, ".+RData", sep="")
              if (length(list.files(object@directory, pattern=file.pattern)) == 0)
                  return(paste(object@directory, " has no .RData files"))
              TRUE
            })


setMethod("getValidChromosomes",
          signature(object="AssocResultsByChr"),
          function(object) {
            file.pattern <- paste0(object@base, object@chromSep, ".+RData", sep="")
            files <- list.files(object@directory, pattern=file.pattern)
            chroms <- matrix(unlist(strsplit(sub("[.][^.]*$", "", files), object@chromSep)), ncol=2, byrow=TRUE)[,2]
            mixedsort(chroms)
          })


setGeneric("getAssocResults", function(object, ...) standardGeneric("getAssocResults"))
setMethod("getAssocResults",
          signature(object = "AssocResultsByChr"),
          function(object, chromosome=getValidChromosomes(object)) {
              chk <- setdiff(chromosome, getValidChromosomes(object))
              if (length(chk) > 0) stop (paste(paste(chk, collapse=","), "is not a valid chromosome"))
              
              assoc.list <- lapply(unique(chromosome), function(c) {
                  getobj(file.path(object@directory,
                                   paste(object@base, object@chromSep, c, ".RData", sep="")))
              })
              x <- do.call(rbind, assoc.list)
              
              ## select non-duplicated snpIDs -- necessary for haplotype2 tests
              sel <- any(duplicated(x$snpID))
              if (sum(sel) > 0) warning("duplicated snpIDs detected. selecting the first non-duplicated snpID.")
              x <- x[!duplicated(x$snpID), ]
              x
          })

## see lookUpSnps.R for definition of lookUpSnps method for AssocResultsByChr
