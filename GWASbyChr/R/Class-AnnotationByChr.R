setClass("AnnotationByChr",
         representation(directory = "character",
                        base = "character",
                        chromSep = "character",
                        suffix = "character"),
         prototype(directory = "",
                   base = "",
                   chromSep="_chr-",
                   suffix=".RData"))


AnnotationByChr <- function(directory, base, chromSep, suffix, ...) {
  if (missing(chromSep)) chromSep <- "_chr-"
  if (missing(suffix)) suffix <- ".RData"
  # autodetect base?
  if (missing(base)) {
    files <- list.files(directory, pattern=paste0(suffix, "$"))
    # probably not the best way to do this, but..
    base <- unique(matrix(unlist((strsplit(files, chromSep))), ncol=2, byrow=TRUE)[,1])
    if (length(base) > 1) stop("ambiguous base name in directory. provide base.") # maybe this belongs in validity method
  }
  new("AnnotationByChr", directory=directory, base=base, chromSep=chromSep, suffix=suffix, ...)
}

# to write -- what does it check?
# - if the directory exists.
# - if the directory has RData files
setValidity("AnnotationByChr", 
            function(object) {
              
              # check that the directory exists
              if (!file.exists(object@directory)) return(paste(object@directory, " does not exist"))
              
              # check that the directory is actually a directory
              if (!file_test("-d", object@directory)) return(paste(object@directory, " does not exist"))
              
              # check that there is at least one RData file
              file.pattern <- paste0("^", object@base, object@chromSep, ".+", object@suffix, "$")
              if (length(list.files(object@directory, pattern=file.pattern)) == 0)
                  return(paste(object@directory, " has no .RData files"))
              
              TRUE
            })


setGeneric("getValidChromosomes", function(object, ...) standardGeneric("getValidChromosomes"))
setMethod("getValidChromosomes",
          signature(object="AnnotationByChr"),
          function(object) {
            file.pattern <- paste0("^", object@base, object@chromSep, ".+", object@suffix, "$")
            files <- list.files(object@directory, pattern=file.pattern)
            chroms <- matrix(unlist(strsplit(sub(paste0(object@suffix, "$"), "", files), object@chromSep)),
                             ncol=2, byrow=TRUE)[,2]
            mixedsort(chroms)
          })


# returns the snp annotation for the specified chromosome
setMethod("getSnpAnnotation",
          signature(object = "AnnotationByChr"),
          function(object, chromosome=getValidChromosomes(object), returnColumns=NULL) {
            if (!all(chromosome %in% getValidChromosomes(object))) stop (paste(chromosome, "is not a valid chromosome"))
            
            snpAnnot.list <- list()
            for (chr in unique(chromosome)){
              snpAnnot <- getobj(file.path(object@directory,
                                           paste0(object@base, object@chromSep, chr, object@suffix)))
              if (!is.null(returnColumns)) {
                keep <- intersect(c("snpID", "chromosome", "position", returnColumns), varLabels(snpAnnot))
                snpAnnot <- snpAnnot[, keep]
              }
              snpAnnot.list[[as.character(chr)]] <- snpAnnot
            }

            snp.list <- lapply(snpAnnot.list, pData)
            meta.list <- lapply(snpAnnot.list, varMetadata)
            
            # rbind snp annotations together
            snp <- do.call(rbind.fill, snp.list)
            snp <- snp[order(snp$snpID), ]
            rownames(snp) <- NULL
            
            # extract out unique meta elements.. this is more complicated
            meta <- meta.list[[1]]
            if (length(meta.list) > 1){
              for (i in 2:length(meta.list)){
                x <- meta.list[[i]]
                x <- x[!(rownames(x) %in% rownames(meta)), , drop=F]
                if (nrow(x) > 1) meta <- cbind(meta, x)
              }
            }
            
            snpAnnot <- SnpAnnotationDataFrame(snp, meta)
            snpAnnot
            
          })

## see lookUpSnps.R for definition of lookUpSnps method for AnnotationByChr
