
setGeneric("lookUpSnps", function(object1, object2, ...) standardGeneric("lookUpSnps"))
setMethod("lookUpSnps",
          signature(object1 = "AnnotationByChr", object2="missing"),
          function(object1, object2, snps=NULL, column="snpID", chromosome=getValidChromosomes(object1), extraCols=NULL, ...) {
            
            ## add option to look up by multiple columns?  (chromosome, position)
            
            # for memory considerations, loop over chromosomes here even though getSnpAnnotation can do it
            annot.list <- list()
            for (c in unique(chromosome)){
              annot <- getSnpAnnotation(object1, c, ...)
              cols <- intersect(c("snpID", "snpName", "chromosome", "position", "alleleA", "alleleB", "type", "rsID", "dataset", "composite.filter", "oevar", extraCols), varLabels(annot))
              
              if (!is.null(snps)){
                index <- annot[[column]] %in% snps
              } else {
                index <- 1:nrow(annot)
              }
              annot.list[[c]] <- pData(annot)[index, cols]
            }
            
            meta <- varMetadata(annot)[cols, , drop=FALSE]
            
            annot <- do.call(rbind, annot.list)
            
            if (nrow(annot) > 0) {
                row.names(annot) <- 1:nrow(annot)
            } else {
                warning(paste("No matching SNPs found in column", column))
            }
            SnpAnnotationDataFrame(annot, meta)
          })



setMethod("lookUpSnps",
          signature(object1 = "AssocResultsByChr", object2 = "missing"),
          function(object1, snps=NULL, column="snpID", chromosome=getValidChromosomes(object1)) {
            
            assoc.list <- list()
            for (c in unique(chromosome)){
              assoc <- getAssocResults(object1, c)
              if (!is.null(snps)) {
                index <- assoc[[column]] %in% snps
              } else {
                index <- 1:nrow(assoc)
              }
              assoc.list[[c]] <- assoc[index, ]
            }

            do.call(rbind, assoc.list)
          })


setMethod("lookUpSnps",
          signature(object = "AnnotationByChr", object2 = "AssocResultsByChr"),
          function(object1, object2, snps=NULL, column="snpID", chromosome=getValidChromosomes(object2), ...) {
            
            annot <- lookUpSnps(object1, snps=snps, column=column, chromosome=chromosome, ...)
            annot <- pData(annot)
            
            assoc <- lookUpSnps(object2, snps=annot$snpID, column="snpID", chromosome=unique(annot$chromosome))
            
            annot$chromosome <- chromToChar(annot$chromosome)
            
            res <- merge(annot, assoc, by="snpID", suffixes=c(".annot", ".assoc"))
            # remove columns that are the same for both
            # annotation takes precedence over association
            for (n in setdiff(intersect(names(assoc), names(annot)), "snpID")){
              n.annot <- paste0(n, ".annot")
              n.assoc <- paste0(n, ".assoc")
              if (allequal(res[[n.annot]], res[[n.assoc]])) res[[n.assoc]] <- NULL
            }
            names(res) <- gsub(".annot", "", names(res), fixed=TRUE)
            
            res$snpID <- as.integer(res$snpID)
            
            if (any(duplicated(res$snpID))) stop("duplicated snpIDs!")
            
            res[order(res$snpID),]
            
          })

