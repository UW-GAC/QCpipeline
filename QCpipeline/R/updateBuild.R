
updateBuild <- function(genoData, map, outPrefix="new",
                        update.alleles=FALSE, remove.unmapped=TRUE,
                        block.size=100) {
         
    gdsfile <- paste0(outPrefix, ".gds")
    gds <- createfn.gds(gdsfile)

    ## add sample data
    varnames <- getVariableNames(genoData@data)
    sample.vars <- varnames[grep("sample", varnames, fixed=TRUE)]
    for (v in sample.vars) {
        add.gdsn(gds, v, getVariable(genoData, v), compress="ZIP.max", closezip=TRUE)
    }
    sync.gds(gds)
        
    ## integer chromosomes
    for (v in c("old.chromosome", "new.chromosome")) {
        map[[v]][map[[v]] %in% "X"] <- 23
        map[[v]][map[[v]] %in% "XY"] <- 24
        map[[v]][map[[v]] %in% "Y"] <- 25
        map[[v]][map[[v]] %in% c("M", "MT", "Mt")] <- 26
        map[[v]][map[[v]] %in% "0"] <- 27
        map[[v]] <- as.integer(map[[v]])
    }
    
    snp <- getSnpVariable(genoData, c("snpID", "chromosome", "position", "rsID"))
    if (update.alleles) {
        stopifnot(all(c("alleleA", "alleleB") %in% names(map)))
    } else {
        snp <- cbind(snp, getSnpVariable(genoData, c("alleleA", "alleleB")))
    }
    names(snp) <- paste0("old.", names(snp))
    snp <- merge(snp, map, all.x=TRUE)
    snp$new.chromosome[is.na(snp$new.chromosome)] <- 27L
    snp$new.position[is.na(snp$new.position)] <- 0L
    snp$new.rsID[is.na(snp$new.rsID)] <- "."
    
    if (remove.unmapped) snp <- snp[snp$new.chromosome != 27,]
    snp <- snp[order(snp$new.chromosome, snp$new.position),]
    snp$new.snpID <- 1:nrow(snp)

    ## index that maps old snp order to new chrom-position order
    ## WARNING: this assumes that rsIDs are unique
    ## will fail if removeUnmapped=FALSE and some old SNPs are not in the map file
    ## TO-DO: fix this if removeUnmapped=FALSE or there are duplicate rsIDs
    snp.index <- match(snp$old.rsID, getSnpVariable(genoData, "rsID"))
    stopifnot(allequal(snp$old.snpID, getSnpID(genoData)[snp.index]))

    ## add snp data
    add.gdsn(gds, "snp.id", snp$new.snpID, compress="ZIP.max", closezip=TRUE)
    add.gdsn(gds, "snp.chromosome", snp$new.chromosome, storage="uint8",
             compress="ZIP.max", closezip=TRUE)
    add.gdsn(gds, "snp.position", snp$new.position, compress="ZIP.max", closezip=TRUE)
    add.gdsn(gds, "snp.allele", paste(snp$alleleA, snp$alleleB, sep="/"),
             compress="ZIP.max", closezip=TRUE)
    add.gdsn(gds, "snp.rs.id", snp$new.rsID, compress="ZIP.max", closezip=TRUE)
    sync.gds(gds)

    nsamp <- nscan(genoData)
    nsnp <- nrow(snp)
    geno.node <- add.gdsn(gds, "genotype", storage="bit2", valdim=c(nsnp, nsamp))
    put.attr.gdsn(geno.node, "snp.order")

    nblocks <- ceiling(nsamp / block.size)
    for (i in 1:nblocks) {
        start <- (i-1)*block.size + 1
        end <- min(i*(block.size), nsamp)
        count <- end - start + 1
        geno <- getGenotype(genoData, scan=c(start, count), snp=c(1,-1))
        geno <- geno[snp.index,]
        geno[is.na(geno)] <- 3
        write.gdsn(geno.node, geno, start=c(1,start), count=c(-1,count))
    }

    closefn.gds(gds)
    cleanup.gds(gdsfile)
    
    cols <- c("snpID", "chromosome", "position", "rsID")
    snp <- snp[,c(paste0("new.", cols), "alleleA", "alleleB")]
    names(snp)[1:4] <- cols
    snpAnnot <- SnpAnnotationDataFrame(snp)
    save(snpAnnot, file=paste0(outPrefix, "_snpAnnot.RData"))
}


