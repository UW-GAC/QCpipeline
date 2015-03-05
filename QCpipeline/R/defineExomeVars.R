# Function to annotate variants, including exomic/non-exomic
# Following VariantAnnotation workflow - http://www.bioconductor.org/help/workflows/annotation/AnnotatingRanges/#annotate-by-position
# SN, UW GAC, Feb 2015

## input:
# 1) snpAnnot - SNP annotation (GWASTools "SnpAnnotationDataFrame" object)
# 2) txdb - a TxDb object with the desired database table annotation - see getTxDb function; alternatively, if the user has already run 'getTxDb' and saved an .sqlite object with the TxDb ojbect, that .sqlite file can be loaded with 'loadDb' and passed to this function argument
# 3) collapsed - T/F argument where T means to return one row per input SNP, with a collapsed list of location types (e.g., "intron; fiveUTR; promoter") and geneIDs (e.g., "126789; 116983"), F means return one row per SNP/location type combination

## output:
# data frame with input SNPs (snpID, chromosome, position), with additional columns:
# 1) excl.chr: logical vector where T = variant was excluded from query because it was located on a chromosome other than autosomes, X, Y, or mitochondrial; F if not excluded based on this criterion
# 2) excl.pos: logical vector where T = variant was excluded from query becuase position=0; F if not excluded based on this criterion
# 3) loctype: Values returned from 'locateVariants'. Possible values are 'coding', 'intron', 'threeUTR', 'fiveUTR', 'intergenic', 'spliceSite', or 'promoter'. If collapsed=TRUE, multiple records are concatenated by ";". 
# 4) geneID: Values returned from locateVariants. For genic SNPs, one or more Entrez GeneIDs; NA otherwise. If collapsed=TRUE, multiple records are concatenated by ";". 
# 5) exomic: logical vector, TRUE where variant is located in exomic region (at least one loctype value of coding); FALSE if not located in exomic region; NA if not included in query (either because excl.chr=TRUE, excl.pos=TRUE, or no match found by locateVariants query)

## # values for testing:
## # rm(list=objects())
## snpAnnot <- getobj("/projects/cidr/Marazita_ofc/sample_snp_annot/Marazita_ofc_HumanCoreExomePlusCustom_Marazita_15050181_A_V13_scn.RData")


################# Start function definition

defineExomeVars <- function(snpAnnot, txdb, collapsed=TRUE) {

  # check required packages
  require(VariantAnnotation)
  require(GenomicFeatures) # probably already loaded in user exection of getTxDb
  require(AnnotationDbi) # probably already loaded in user exection of getTxDb

  ## reshape2 only needed if collapsing results to one row per SNP
  # if(collapsed) {require(reshape2)}

    # check primary SNP annotation object
    stopifnot(class(snpAnnot) == "SnpAnnotationDataFrame")
    nvarAll <- nrow(snpAnnot)
    message(paste("Read in", prettyNum(nvarAll, big.mark=","), "variants from primary annotation", "\n"))

    # getVariables from snpAnnot
    snp.dat <- getVariable(snpAnnot, c("snpID", "chromosome", "position"))

    # convert integer chromosome codes to characters
    snp.dat$chrom.char <- snp.dat$chromosome
    snp.dat$chrom.char[snp.dat$chromosome %in% XchromCode(snpAnnot)] <- "X"
    snp.dat$chrom.char[snp.dat$chromosome %in% XYchromCode(snpAnnot)] <- "X"
    snp.dat$chrom.char[snp.dat$chromosome %in% YchromCode(snpAnnot)] <- "Y"
    snp.dat$chrom.char[snp.dat$chromosome %in% MchromCode(snpAnnot)] <- "M"

    # exclude SNPs not on autosomes, X, Y, M
    snp.dat$excl.chr <- !is.element(snp.dat$chrom.char,
                                    c(autosomeCode(snpAnnot),"X","Y","M"))
    # exclude unmapped SNPs
    snp.dat$excl.pos <- snp.dat$position %in% 0
    message("Excluding ", sum(snp.dat$excl.chr) + sum(snp.dat$excl.pos)," SNPs with either unknown position and/or not mapped to autosomes, chrX, chrY, or chrM")

    snp.use <- snp.dat[!snp.dat$excl.pos & !snp.dat$excl.chr,]
    # to test
    # snp.use <- snp.use[1:10000,]

    snp.use$chrom.code <- paste0("chr",snp.use$chrom.char)
    
    # create a GRanges object for all the SNP positions
    gr.query <- GRanges(snp.use$chrom.code, IRanges(snp.use$position, snp.use$position))

    message("Printing information about specified 'txdb' object:")
    print(txdb)

    # to look at chromosome names
    # intersect(seqlevels(txdb), seqlevels(gr.query))

    # assign the genome definition
    # genome(gr.query) <- unique(genome(txdb))

    # run locateVariants
    message("\nRunning locateVariants on ", prettyNum(nrow(snp.use),big.mark=","), " variants")
    loc <- locateVariants(gr.query, txdb, AllVariants())

    # create data frame of results
    # convert factor to character
    loctypes <- as.character(mcols(loc)$LOCATION)
    loc.rslt <- unique(data.frame(queryID=mcols(loc)$QUERYID,
                                  loctype=loctypes,
                                  geneID=mcols(loc)$GENEID,
                                  stringsAsFactors=FALSE))
    message("More than one result for ", prettyNum(sum(duplicated(loc.rslt$queryID)), big.mark=","), " variants - processing...")

    # queryID provides map back to row in original query
    # assign snpID to results df
    snp.use$rowID <- 1:nrow(snp.use)
    loc.rslt$snpID <- snp.use$snpID[match(loc.rslt$queryID, snp.use$rowID)]

    # if not collapsing results, match loc.rslt to snp.dat
    if (!collapsed) {
      out <- merge(snp.dat, loc.rslt[,c("snpID","loctype","geneID")],
                   by="snpID", all.x=TRUE, all.y=TRUE, sort=FALSE)
    }


    # if collapsing results, rehape the data frame
    # account for where there is more than one result per SNP    
    if (collapsed) {

      # for testing:
      # dim(loc.rslt); class(loc.rslt$queryID); class(loc.rslt$loctype); class(loc.rslt$geneID)
      # loc.rslt <- loc.rslt[is.element(loc.rslt$snpID, 1:100),]

      # melt the data frame
      loc.melt <- melt(loc.rslt, id=(c("queryID", "snpID")), na.rm=TRUE)

      # cast the data frame
      collapsed.rslt <- dcast(loc.melt, queryID + snpID ~ variable,
                    fun.aggregate = function(x) paste(x, collapse = "; "), fill = "NA")

      # merge back in with snp.dat
      out <- merge(snp.dat, collapsed.rslt[,c("snpID","loctype","geneID")],
                   by="snpID", all.x=TRUE, all.y=TRUE, sort=FALSE)

    } # end 'if' on collapse=FALSE

    # post-processing for both collapsed and uncollapsed data frames

    # 1) sanity check via position
    ir.rslt <- ranges(loc) #IRanges object
    pos.rslt <- unique(data.frame(queryID=loc$QUERYID,
                                  position=start(ir.rslt)))

    # line up query ID, snpID, position
    chk.pos <- merge(snp.use, pos.rslt, by.x="rowID", by.y="queryID")
    if(!allequal(chk.pos$position.x, chk.pos$position.y))
      stop("Problem! Row number of input SNPs does not match QUERYID in locateVariants result. Report to code developer") # need to find out something useful to do here

    # 2) annotated "exomic" where any of the loctypes are "coding"
    out$exomic <- grepl("coding",out$loctype)

    # set to 'NA' where loctype is NA
    out$exomic[is.na(out$loctype)] <- NA
    
    # 3) sort by snpID
    out.srt <- out[order(out$snpID),setdiff(names(out),"chrom.char")]

    message("Finished!")

    return(out.srt)

} # end function definition
