# Function to take a SNP annotation object and identify duplicate variants
# Assess overlap in chromosome, position, and - optionally - alleles
# Allows for strand flips - assigned duplicate status and flagged with strand discrepancy
# Written for construction of SoL freeze 2: combining imputed 1, imputed 2, and Metabochip imputation
# Likely to go into QCpipeline
# Can be extended to check for triallelics and assign a "triallelic.var.id"
# SN, 10/23/14

# Allows for supplementary variant-level annotation containing alternate information (i.e., position and alleles) to use for the look-up. Useful for insertion/deletion variants (indels) - e.g., if some indels in primary annotation are annotated using Illumina's representation of position and alleles, whereas other - such as 1000G imputed - indels in primary annotation are annotated with VCF convention for position and alleles. Could also be used to incorporate converted build information 

#################
# Input (required):
# snpAnnot - SNP annotation data frame object. 
# alleleAcol - first allele column to match on, defaults to "alleleA"
# alleleBcol - second allele column to match on, defaults to "alleleB"
# (Note that order of alleles does not matter)
# match.alleles - if FALSE, match only on chromosome and position. if TRUE, additionally match on alleles (allowing for strand discrepancies). Defaults to FALSE (data cleaning standards; for SoL we know we want =TRUE)

# Input (optional, for use with supplementary annotation object):
# supp_snpAnnot - supplementary SNP annotation data frame object. Assumes that "chromosome" and "position" fields contain the position information to use for matching
#  **supplemental snpID needs to match primary snpID
# supp_alleleAcol - first allele column to match on in supplementary SNP annotation, defaults to "alleleA"
# supp_alleleBcol - first allele column to match on in supplementary SNP annotation, defaults to "alleleB"
# supp_snpIDs - list of snpIDs for which to use position and allele information in the supplementary annotation when determining overlap. Defaults to use supplementary annotation for all snpIDs found in common between primary and supplementary annotationDefaults to use supplementary annotation for all snpIDs found in common between primary and supplementary annotation

# Output:
# annotation data frame with additional columns:
# dup.pos.id - an integer (1:n) assigned to pairs (or multiples) of variants deemed to be equivalent based on position alone
# dup.var.id - returned if match.alleles=TRUE, an integer (1:n) assigned to pairs (or multiples) variants deemed to be equivalent based on position and alleles
# dupDiscStrand - returned if match.alleles=TRUE, logical vector flagging where overlapping variants appeared to be represented on a different strand (e.g., indicating difference in plus strand alleles between 1000G and internal reference panel imputations)
# suppAnnotUsed - returned if supp_snpAnnot!=NULL, a logical vector flagging where supplementary annotation was used to determine overlap

################# 

## ## values for testing: rm(list=objects()) options(stringsAsFactors = FALSE)
## library(GWASTools) # will be loaded via NAMESPACE library(Biostrings) #
## snpAnnot should be object already loaded in snpAnnot <-
## getobj('/projects/geneva/gcc-fs2/OLGA/genotype/freeze2/amstilp/results/gds/freeze2_v01/freeze2_snpAnnot_chr22_tmp.RData')
## # supp_snpAnnot should be object already loaded in supp_snpAnnot <-
## getobj('/projects/geneva/gcc-fs2/OLGA/genotype/freeze2/sarahcn/results/alleleMappings/SoL_freeze2_indels_VCFconventions.RData')

## #### adjust snpID to match primary annot snpID (temporary kluge) supp.dat <-
## pData(supp_snpAnnot) supp.dat$snpID <- NA # get indel positions on chr22
## snp.dat <- pData(snpAnnot) snp.dat$indel <- nchar(snp.dat$alleleA)>1 |
## nchar(snp.dat$alleleB)>1 indels <- snp.dat[snp.dat$indel,] supp.dat <-
## supp.dat[supp.dat$chromosome==22,] supp.dat$snpID <-
## indels$snpID[match(supp.dat$position, indels$position)]
## sum(!is.na(supp.dat$snpID)) # filled in 18 - reduce to just those supp.dat <-
## supp.dat[!is.na(supp.dat$snpID),]; dim(supp.dat) supp_snpAnnot <-
## SnpAnnotationDataFrame(supp.dat, YchromCode=as.integer(24))
## YchromCode(supp_snpAnnot) # I assigned Y chrom variants to chromosome==24 ###

## supp_snpIDs <- supp.dat$snpID[1:5] supp_alleleAcol <- alleleAcol <- 'alleleA'
## supp_alleleBcol <- alleleBcol <- 'alleleB'

################# Start function definition

defineDupVars <- function(snpAnnot, alleleAcol = "alleleA", alleleBcol = "alleleB", 
    supp_snpAnnot = NULL, supp_alleleAcol = "alleleA", supp_alleleBcol = "alleleB", 
    supp_snpIDs = NULL, match.alleles=FALSE, start_dup.pos.id=1, start_dup.var.id=1) {
    
    # check primary SNP annotation object
    stopifnot(class(snpAnnot) == "SnpAnnotationDataFrame")
    nvarAll <- nrow(snpAnnot)
    message(paste("Read in", nvarAll, "variants from primary annotation", "\n"))
    
    # getVariables from snpAnnot
    snp.dat <- getVariable(snpAnnot, c("snpID", "chromosome", "position", alleleAcol, 
        alleleBcol))
    
    # standardize column names
    names(snp.dat)[4:5] <- c("alleleA", "alleleB")
    
    # if provided, read in supplementary annotation object
    if (!is.null(supp_snpAnnot)) 
        {
            stopifnot(class(supp_snpAnnot) == "SnpAnnotationDataFrame")
            nsuppvar <- nrow(supp_snpAnnot)
            message(paste("Read in", nsuppvar, "variants from supplementary annotation\n"))
            supp.dat <- getVariable(supp_snpAnnot, c("snpID", "chromosome", "position", 
                supp_alleleAcol, supp_alleleBcol))
            
            # standardize column names
            names(supp.dat)[4:5] <- c("alleleA", "alleleB")
            
            ov <- intersect(supp.dat$snpID, snp.dat$snpID)
            message(paste("\tFound", length(ov), "variants overlapping with primary annotation, matching on snpIDs\n"))
            
            # reduce to overlapping variants
            supp.dat <- supp.dat[is.element(supp.dat$snpID, ov), ]

            # check if we're further subsetting the supplementary annotation based on user input
            if (!is.null(supp_snpIDs)) 
                {
                  
                  # reduce specified list of variants to those in primary annotation
                  ov <- intersect(supp_snpIDs, snp.dat$snpID)
                  supp.dat <- supp.dat[is.element(supp.dat$snpID, ov), ]
                  message("\tUsing supplementary annotation for subset of ", length(ov), 
                    " specified variants \n")
                }  # close if on specifying select snpIDs from supplementary annot
            
            # replace information in primary annotation with supplementary annotation
            snp.keep <- snp.dat[!is.element(snp.dat$snpID, supp.dat$snpID), ]
            
            # first check non-autosome chrom coding - conform to primary annotation
            if (XchromCode(snpAnnot) != XchromCode(supp_snpAnnot)) 
                {
                  message("\tharmonizing chrX integer code")
                  supp.dat$chromosome[is.element(supp.dat$chromosome, XchromCode(supp_snpAnnot))] <- XchromCode(snpAnnot)
                }  # end X chrom check
            
            if (YchromCode(snpAnnot) != YchromCode(supp_snpAnnot)) 
                {
                  message("\tharmonizing chrY integer code")
                  supp.dat$chromosome[is.element(supp.dat$chromosome, YchromCode(supp_snpAnnot))] <- YchromCode(snpAnnot)
                }  # end Y chrom check
            
            if (XYchromCode(snpAnnot) != XYchromCode(supp_snpAnnot)) 
                {
                  message("\tharmonizing chrXY integer code")
                  supp.dat$chromosome[is.element(supp.dat$chromosome, XYchromCode(supp_snpAnnot))] <- XYchromCode(snpAnnot)
                }  # end XY chrom check
            
            if (MchromCode(snpAnnot) != MchromCode(supp_snpAnnot)) 
                {
                  message("\tharmonizing chrM integer code")
                  supp.dat$chromosome[is.element(supp.dat$chromosome, MchromCode(supp_snpAnnot))] <- MchromCode(snpAnnot)
                }  # end M chrom check
            
            # with chromosome codes aligned, merge the two dataframes flag where positional
            # information from the supplementary annotation was used
            snp.keep$suppAnnotUsed <- FALSE
            supp.dat$suppAnnotUsed <- TRUE
            snp.comb <- rbind(snp.keep, supp.dat)
            
            # re-order by snpID
            snp.dat <- snp.comb[order(snp.comb$snpID), ]
        }  # close if on supplementary annot
    
    #### preliminary duplicate definition: matching on position alone
    
    # recode XY variants as X. as long as position refers to chrX position, the
    # distinction between non-PAR/XTR and PAR/XTR regions of the X doesn't affect
    # detection of duplicate positions. In fact it may thwart it, if annotations
    # haven't been consistently annotated wrt to PAR/XTR chrX locations
    snp.dat$chromosome[is.element(snp.dat$chromosome, XYchromCode(snpAnnot))] <- XchromCode(snpAnnot)
    
    # concatenate chrom and bp position into variable 'map'
    snp.dat$map <- paste(snp.dat$chromosome, snp.dat$position)
    # exclude unmapped variants (where either chrom or position is unknown).  assumes
    # 27 as unknown chrom code
    dup.pos <- unique(snp.dat$map[duplicated(snp.dat$map) & !is.element(snp.dat$map, 
        "27 0")])
    ndup <- length(dup.pos)
    nvar <- sum(is.element(snp.dat$map, dup.pos))
    message("Based on position alone, detected ", prettyNum(ndup, big.mark = ","), 
        " duplicated positions involving ", prettyNum(nvar, big.mark = ","), " variants\n")
    
    # assign preliminary dup.var.id based on position alone
    dup.pos.ids <- start_dup.pos.id:(length(dup.pos) + start_dup.pos.id - 1)
    dup.dat <- data.frame(map.pos = dup.pos, id = dup.pos.ids)
    snp.dat$dup.pos.id <- dup.dat$id[match(snp.dat$map, dup.dat$map.pos)]
    
    # refine matching based on alleles, if match.alleles=TRUE
    if (match.alleles) 
        {
            message("\nContinuing to refine position-based matching based on alleles...\n")
            
            # subset snp.dat to the duplicated positions
            snp.dup <- snp.dat[!is.na(snp.dat$dup.pos.id), ]
            
            # sort alleles alphabetically
            snp.dup$alleles <- pasteSorted(snp.dup$alleleA, snp.dup$alleleB)
            
            ## identify where sorted alleles match with consistent strand designations
            snp.dup$pos.alleles <- paste(snp.dup$map, snp.dup$alleles)
            dup.pos.alleles <- unique(snp.dup$pos.alleles[duplicated(snp.dup$pos.alleles)])
            
            ndup <- length(dup.pos.alleles)
            nvar <- sum(is.element(snp.dup$pos.alleles, dup.pos.alleles))
            message("\tMatching on alleles in original strand orientations, detected ", 
                prettyNum(ndup, big.mark = ","), " duplicated position-alleles involving ", 
                prettyNum(nvar, big.mark = ","), " variants\n")
            
            # assign dup.origalles.id so we can later identify strand dicrepancies
            dup.dat <- data.frame(map.pos = dup.pos.alleles, id = 1:length(dup.pos.alleles))
            snp.dup$dup.origalleles.id <- dup.dat$id[match(snp.dup$pos.alleles, dup.dat$map.pos)]
            
            ## identify where sorted alleles appear to be strand flip using Biostrings
            ## package, turn alleles into DNAStrings to get reverse complements only convert
            ## to DNAStrings those alleles that are actually DNA nucelotides becuase of
            ## multi-nucleotide alleles, can't select alleles based on being mamber of
            ## DNA_ALPHABET for now, manually exluding non-DNA allele codes that we're aware
            ## of (not an ideal fix for durability and portability) when one of the alleles
            ## is 0, variant is monomorphic and we can't determine what the minor allele is.
            ## conservative approach is to NOT match up on major allele alone
            
            # don't try to take reverse complement of alleles annotated as 'I','D', or '0'
            take.comp <- !is.element(snp.dup$alleleA, c("I", "D", "0", "<DEL>")) & 
                !is.element(snp.dup$alleleB, c("I", "D", "0", "<DEL>"))
            
            not.dna <- sum(!take.comp)
            message("\tNot taking reverse complement of alleles at ", not.dna, " variants that have I, D, 0, or <DEL> as alleles\n")
            
            # eventually wrap conversion to DNAStringSet with 'try' or 'tryCatch' so that an
            # error will cause script to exit yet still return SNP annotation with dup.var.id
            
            a.str <- DNAStringSet(snp.dup$alleleA[take.comp])
            b.str <- DNAStringSet(snp.dup$alleleB[take.comp])
            
            alleleA.flip <- reverseComplement(a.str)
            alleleB.flip <- reverseComplement(b.str)
            
            # only define flipped alleles where both alleles could be processed as a
            # DNAString
            snp.dup$alleleA.flip <- snp.dup$alleleB.flip <- NA
            snp.dup$alleleA.flip[take.comp] <- as.character(alleleA.flip)
            snp.dup$alleleB.flip[take.comp] <- as.character(alleleB.flip)
            
            # create a string of position, alleles (original), and alleles (reverse
            # complemented)
            snp.dup$alleles.revcomp <- pasteSorted(snp.dup$alleleA.flip, snp.dup$alleleB.flip)
            snp.dup$pos.alleles.revcomp <- paste(snp.dup$map, pasteSorted(snp.dup$alleles, snp.dup$alleles.revcomp, sep=" "))
            
            # check for duplicates within position, alleles, and reverse complemented alleles
            dup.pos.alleles <- unique(snp.dup$pos.alleles.revcomp[duplicated(snp.dup$pos.alleles.revcomp)])
            
            ndup <- length(dup.pos.alleles)
            nvar <- sum(is.element(snp.dup$pos.alleles.revcomp, dup.pos.alleles))
            message("\tFurther allowing for reverse-complemented alleles, detected ", 
                prettyNum(ndup, big.mark = ","), " duplicated position-allele combinations involving ", 
                prettyNum(nvar, big.mark = ","), " variants\n")
            
            ## assign dup.var.id based on positions and alleles
            dup.var.ids <- start_dup.var.id:(length(dup.pos.alleles) + start_dup.var.id - 1)
            dup.dat <- data.frame(map.pos = dup.pos.alleles, id = dup.var.ids)
            snp.dup$dup.var.id <- NA
            snp.dup$dup.var.id <- dup.dat$id[match(snp.dup$pos.alleles.revcomp, dup.dat$map.pos)]
            
            # detect where allele matched but with strand discrepancy
            snp.dup$dupDiscStrand <- is.na(snp.dup$dup.origalleles.id) & !is.na(snp.dup$dup.var.id)
            
            # return to snp.dat and add dup.var.id and dupDiscStrand
            snp.dat$dup.var.id <- snp.dup$dup.var.id[match(snp.dat$snpID, snp.dup$snpID)]
            snp.dat$dupDiscStrand <- snp.dup$dupDiscStrand[match(snp.dat$snpID, snp.dup$snpID)]
            snp.dat$dupDiscStrand[is.na(snp.dat$dupDiscStrand)] <- FALSE
            
        }  # end matching on alleles
    
    # return annotated data frame
    cols.exclude <- c("map")
    snp.return <- snp.dat[, setdiff(names(snp.dat), cols.exclude)]
    return(snp.return)
}  # end function definition 
