# Function to check consistency between probe sequences and A/B alleles in an Illumina manifest
# Written after discovering inconsistencies in HumanOmni2.5-4v1_D.csv array manifest used in HRS1-2 project, which caused many errors (1) when merging with HRS3 and (2) imputing HRS1-2 to 1000 Genomes

# SN 5/30/2014

# arguments:
# 1: data frame of Illumina SNP annotation, including columns an identifier,"SNP","AlleleA_ProbeSeq," and "AlleleB_ProbeSeq"
# 2: name of field to use as identifier (defaults to "Name")
# 3: logical for whether to report on function progress
# returns a subset of input data frame with discrepancies 

# example usage:
## illumina <-  read.csv(file="/projects/geneva/gcc-fs2/SNP_annotation/Illumina/HumanOmni2.5_4v1/HumanOmni2.5-4v1_D.csv", skip=7, nrows=50000)
## out <- checkABalleles(illumina)

##############

#require(stringr)

checkABalleles <- function(snp,identifier="Name",verbose=TRUE) {

  # check for presence of required fields
  if(!all(c(identifier,"SNP","AlleleA_ProbeSeq","AlleleB_ProbeSeq") %in% names(snp)))
    stop("Check required fields for input file")

  # stop if identifier is not unique, as we plan to use it as identifier in output
  if(sum(duplicated(snp[,names(snp)==identifier]))>0)
    stop("'Name' field is not unique as expected")
    
  nsnps.init <- nrow(snp)
  
  # reduce to variants with two probe types (i.e., InfiniumI assays)
  snp <- snp[!is.element(snp$AlleleB_ProbeSeq,"") &
             !is.na(snp$AlleleB_ProbeSeq),]

  # check that Infinium I assays are all strand-ambiguous SNPs
  if(sum(!is.element(snp$SNP,c("[A/T]","[C/G]", "[G/C]","[T/A]")))>0)
    warning("Not all two-probe (InfiniumI) assays appear to be strand ambiguous, which is not as expected")

  # report out initial and reduced dimensions
  if(verbose){
    nsnps <- nrow(snp)
    message(nsnps," SNPs with two probe sequences, out of ",nsnps.init," in input")
    message("Determining consistency between A/B alleles in 'SNP' column and terminal nucleotide of the probe sequences")
  }

  # get last nucleotide of ProbeSeqs
  snp$probeA <- str_sub(snp$AlleleA_ProbeSeq,start=-1)
  snp$probeB <- str_sub(snp$AlleleB_ProbeSeq,start=-1)

  # get design alleles
  snp$designA <- substr(snp$SNP, 2,2)
  snp$designB <- substr(snp$SNP, 4,4)

  # flag discrepancies
  snp$probes.disc <- (snp$designA!=snp$probeA) | (snp$designB!=snp$probeB)

  if (verbose){
    ndisc <- sum(snp$probes.disc)
    message(ndisc," SNPs with discrepancies out of ",nsnps," checked")
  }
  
  # return reduced data frame with discrepant SNPs
  results <- snp[snp$probes.disc,
                 c(identifier,"SNP","AlleleA_ProbeSeq","AlleleB_ProbeSeq","probeA","probeB","designA","designB")]

  return(results)
}

