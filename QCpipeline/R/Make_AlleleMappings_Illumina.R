## Function to take Illumina build 37 SNP annotation data file and make an Allele Mappings table
## Sarah Nelson, UWGCC, May 25, 2011
## Updated March 7, 2012, to allow for indel records where IlmnStrand != SourceStrand

################################# CONTENT
## Example allele mappings table:
#                 snp alle.AB alle.design alle.top alle.fwd alle.plus
#  rs1000000       A           T        A        T         A
#  rs1000000       B           C        G        C         G
#  rs1000002       A           A        A        A         T
#  rs1000002       B           G        G        G         C
#  rs10000023       A           T        A        T         T
#  rs10000023       B           G        C        G         G

## With the following data dictionary:
#  variable        type    description
#  snp     text    rs id and other comparable identifiers for a snp
#  alle.AB text    A or B, per the Illumina genotyping system nomenclature
#  alle.design     text    nucleotide(s) corresponding to A or B allele for design strand
#  alle.top        text    nucleotide(s) corresponding to A or B allele for Illumina TOP strand
#  alle.fwd        text    nucleotide(s) corresponding to A or B allele for FORWARD strand, with respect to dbSNP refSNP exemplar
#  alle.plus       text    nucleotide(s) corresponding to A or B allele for PLUS(+) strand, relative to the forward direction in the human reference genome sequence

## Note that "RefStrand" in build 37 annotations from Illumina is result of BLAST search of DESIGN strand
## In previous strand annotations from Illumina for build 36 arrays, "BlastStrand" held +/- result of  BLAST search of the SOURCE sequence

################################# USAGE
## Requires as input:
# "snp.dat", the name of an data frame made from SNP Illumina annotation file (i.e., "HumanOmni2.5-4v1_D.csv"), with following fields:
#  "IlmnID"                "Name"                  "IlmnStrand"           
#  "SNP"                   "GenomeBuild"          "SourceStrand"                      
#  "SourceSeq" "TopGenomicSeq"         "RefStrand"             

## Returns as output: R data frame object with columns ("snp", "alle.AB", "alle.design", "alle.top", "alle.fwd", "alle.plus"

#Example of using:
#  source("/projects/geneva/geneva_sata/GCC_code/Imputation/R_functions/Make_AlleleMappings_Illumina.R")

##  Load in Illumina annotation
#  column.select=rep("NULL",times=21)
#  column.select[c(1:4,9:11,16:18,21)] <- NA  ## only read in select columns    
#  snp.dat <- read.csv(file="/projects/geneva/geneva_sata/SNP_annotation/Illumina/HumanOmni2.5_4v1/HumanOmni2.5-4v1-Multi-B.csv"),
#                   skip=7,colClasses=column.select,nrows=2450000)   ## last rows are for control probes
#  map.final <- make.allele.mappings(snp.dat)
#	 write.csv(map.final, file="/projects/geneva/geneva_sata/SNP_annotation/Illumina/HumanOmni2.5_4v1/testfn.csv", row.names=FALSE,quote=FALSE))

################################# FUNCTION

make.allele.mappings <- function(snp.dat)
{

  options(stringsAsFactors = FALSE)

  ## Report total number of SNPs and indels
  cat("\tTotal number of probes:", length(snp.dat$Name),"\n")
  cat("\tTotal number of indels:", length(snp.dat$Name[is.element(snp.dat$SNP, c("[D/I]","[I/D]"))]),"\n")
  
  ## Check that SNP name is unique and appropriate to use as the identifier in allele mappings
  if (any(duplicated(snp.dat$Name)))
  {
  	stop(sum(duplicated(snp.dat$Name)), " SNP names are duplicated; need alternate identifier to make Allele Mappings file")
  }
  
  ## set up table template - 2 rows per SNP, first "A", second "B"
  map <- as.data.frame(c(snp.dat$Name, snp.dat$Name))
  names(map) <- "snp"
  snp.dat$order <- 1:length(snp.dat$Name)
  
  ## reorder data frame to match SNP annotation
  map$order <- snp.dat$order[match(map$snp, snp.dat$Name)]
  map <- map[order(map$order),]
  
  ## make first instance of SNP name allele A; second instance allele B
  map$alle.AB[!duplicated(map$snp)] <- "A"
  map$alle.AB[duplicated(map$snp)] <- "B"
  
  ## parse A/B alleles from "SNP" column
  snp.dat$design.A <-  substr(snp.dat$SNP,2,2)
  snp.dat$design.B <-  substr(snp.dat$SNP,4,4)
  
  ## if CNV, change design.A and design.B to NA
  snp.dat$design.A[snp.dat$SNP=="[N/A]"] <- NA
  snp.dat$design.B[snp.dat$SNP=="[N/A]"] <- NA
  
  ## extract the Illumina info for the indels, if present
  indels.dat <- snp.dat[is.element(snp.dat$design.A,c("I","D")),]
  
  indels <- indels.dat$Name
  
  ### Make conditional on presence of indels
  if(length(indels)>0){
  cat("\tCreating mappings for insertion/deletion probes\n")
  	## extract the map data frame
	indels.map <- map[is.element(map$snp, indels),]

	## get alleles from the SourceSeq column; makes a list
	sourceseq.start <- strsplit(indels.dat$SourceSeq,"[",fixed=TRUE)

	## get just the 2nd element of each list
	tmp <- sapply(sourceseq.start, "[[", 2)

	sourceseq.end <- strsplit(tmp,"]",fixed=TRUE)

	## indels.snp holds a vector like: [1] "-/GA"     "-/CTCT"   "-/GAAGT"
	indels.snp <- sapply(sourceseq.end, "[[", 1)

	## ready to define possible alleles parsed from SourceSeq column
	## SourceSeq USUALLY presents deletion first; order of A/B assignment indicated by SNP column D/I vs I/D
	alles <- strsplit(indels.snp,"/",fixed=TRUE)
	ins <- sapply(alles, "[[", 2)
	dels <- sapply(alles, "[[", 1)

	## bind "ins" and "del" columns together
	a <- as.data.frame(cbind(ins, dels))
	a$var <- ifelse(a$ins=="-", a$dels, a$ins)

	## Get reverse complement of the variation
	## make DNAstrings object
	var.dstr <- DNAStringSet(a$var)
	## get reverse complement of var
	a$var.revcomp <- as.character(reverseComplement(var.dstr))

	## map the indel names back on
	a$id <- indels

	## identify the records where "ins" is actually the deletion character
	a$deletion.first <- indels.dat$SNP=="[D/I]"

	## identify records where SourceStrand==IlmnStrand
	a$strand.diff <- indels.dat$IlmnStrand!=indels.dat$SourceStrand

	## Assign design alleles, taking into account both "[D/I]" vs "[I/D]" and whether IllmnStrand equals SourceStrand

	indels.dat$alle1[!a$deletion.first & !a$strand.diff] <- a$var[match(indels.dat$Name[!a$deletion.first & !a$strand.diff], a$id)]
	indels.dat$alle1[!a$deletion.first & a$strand.diff] <- a$var.revcomp[match(indels.dat$Name[!a$deletion.first & a$strand.diff], a$id)]

	indels.dat$alle1[a$deletion.first] <- "-"
	indels.dat$alle2[!a$deletion.first] <- "-"  

	indels.dat$alle2[a$deletion.first & !a$strand.diff] <- a$var[match(indels.dat$Name[a$deletion.first & !a$strand.diff], a$id)]
	indels.dat$alle2[a$deletion.first & a$strand.diff] <- a$var.revcomp[match(indels.dat$Name[a$deletion.first & a$strand.diff], a$id)]  

	## MINUS equiv to Illumina BOT; PLUS equiv to TOP

	## Where design strand is PLUS, TOP alles A and B are alles 1, 2 
	indels.dat$top.A[is.element(indels.dat$IlmnStrand,c("P","PLUS"))] <- indels.dat$alle1[is.element(indels.dat$IlmnStrand,c("P","PLUS"))]
	indels.dat$top.B[is.element(indels.dat$IlmnStrand,c("P","PLUS"))] <-indels.dat$alle2[is.element(indels.dat$IlmnStrand,c("P","PLUS"))]

	## make DNAstrings objects for alle1, alle2
	alle1.dstr <- DNAStringSet(indels.dat$alle1)
	alle2.dstr <- DNAStringSet(indels.dat$alle2)

	## get reverse complement of alle1, alle2
	indels.dat$alle1.revcomp <- as.character(reverseComplement(alle1.dstr))
	indels.dat$alle2.revcomp <- as.character(reverseComplement(alle2.dstr))

	## where design sequence is MINUS (i.e., BOT) top alleles are reverse complement of alle1, 2
	indels.dat$top.A[is.element(indels.dat$IlmnStrand,c("M","MINUS"))]<- indels.dat$alle1.revcomp[is.element(indels.dat$IlmnStrand,c("M","MINUS"))]
	indels.dat$top.B[is.element(indels.dat$IlmnStrand,c("M","MINUS"))] <- indels.dat$alle2.revcomp[is.element(indels.dat$IlmnStrand,c("M","MINUS"))]

	## identify the DESIGN alleles 
	indels.dat$design.A <- indels.dat$alle1
	indels.dat$design.B <- indels.dat$alle2

  }
  
  #### exit indels, return to snp.dat
  ## Allow for no indels in SNP.dat
  if (length(indels)>0)   {
  snp.dat <- snp.dat[!is.element(snp.dat$Name, indels),] }
  
  ## for SNPs, make lookup to table to define reverse complement of alleles A, B
  lup <- as.data.frame(c("A","C","G","T"))
  names(lup) <- "b1"
  lup$b2 <- c("T","G","C","A") ## reverse complement
  
  snp.dat$alle1.revcomp <- lup$b2[match(snp.dat$design.A, lup$b1)]
  snp.dat$alle2.revcomp <- lup$b2[match(snp.dat$design.B, lup$b1)]
  
  cat("\tDefining TOP alleles\n")
  
  ## get TOP alleles.  Where IlmnStrand=TOP (IlmnID contains "_T_"), TOP alleles are same as design alleles
  ## Where IlmnStrand=BOT (IlmnID contains "_B_"), TOP alleles are reverse complement of design alleles
  snp.dat$top.A[is.element(snp.dat$IlmnStrand,c("TOP","P"))] <- snp.dat$design.A[is.element(snp.dat$IlmnStrand,c("TOP","P"))]
  snp.dat$top.B[is.element(snp.dat$IlmnStrand,c("TOP","P"))] <- snp.dat$design.B[is.element(snp.dat$IlmnStrand,c("TOP","P"))]
  
  snp.dat$top.A[is.element(snp.dat$IlmnStrand,c("BOT","M"))] <- snp.dat$alle1.revcomp[is.element(snp.dat$IlmnStrand,c("BOT","M"))]
  snp.dat$top.B[is.element(snp.dat$IlmnStrand,c("BOT","M"))] <- snp.dat$alle2.revcomp[is.element(snp.dat$IlmnStrand,c("BOT","M"))]
  
  ## now that design and TOP alleles are defined, rejoin indels with main SNP annotation, keeping select columns
  cols.keep <- c("IlmnID","Name","IlmnStrand","SNP","SourceStrand","RefStrand",
                 "design.A","design.B","alle1.revcomp","alle2.revcomp","top.A","top.B","order")
  
  ## allow for arrays with no indels - in those cases, just keep main snp.dat object with selected columns
  if (length(indels)>0)  {
    comb <- rbind(snp.dat[,cols.keep], indels.dat[,cols.keep]) } 
  else   {comb <- snp.dat[,cols.keep]  }

  ## restore original order, rename to snp.dat
  snp.dat <- comb[order(comb$order),]
  
  cat("\tDefining FORWARD alleles\n")
  
  ## get FORWARD alleles: make column to indicate dbSNP orientation of DESIGN strand
  snp.dat$dbSNPStrand.fordesign[is.element(snp.dat$IlmnID,grep("_F_",snp.dat$IlmnID, value=TRUE))] <- "FWD"
  snp.dat$dbSNPStrand.fordesign[is.element(snp.dat$IlmnID,grep("_R_",snp.dat$IlmnID, value=TRUE))] <- "REV"
  
  snp.dat$fwd.A[snp.dat$dbSNPStrand.fordesign=="FWD"] <- snp.dat$design.A[snp.dat$dbSNPStrand.fordesign=="FWD"]
  snp.dat$fwd.B[snp.dat$dbSNPStrand.fordesign=="FWD"] <- snp.dat$design.B[snp.dat$dbSNPStrand.fordesign=="FWD"]
  
  snp.dat$fwd.A[snp.dat$dbSNPStrand.fordesign=="REV"] <- snp.dat$alle1.revcomp[snp.dat$dbSNPStrand.fordesign=="REV"]
  snp.dat$fwd.B[snp.dat$dbSNPStrand.fordesign=="REV"] <- snp.dat$alle2.revcomp[snp.dat$dbSNPStrand.fordesign=="REV"]
  
  cat("\tDefining PLUS(+) alleles\n")
  
  ## get + alleles, wrt human genome reference sequence
  ## In Illumina b37, RefStrand is BLAST result of Design Strand
  ## where RefStrand="+", same as DESIGN alleles
  ## where RefStrand="-", reverse complement of DESIGN alleles
  
  snp.dat$plus.A[is.element(snp.dat$RefStrand,"+")] <- snp.dat$design.A[is.element(snp.dat$RefStrand,"+")]
  snp.dat$plus.B[is.element(snp.dat$RefStrand,"+")] <- snp.dat$design.B[is.element(snp.dat$RefStrand,"+")]
  
  snp.dat$plus.A[is.element(snp.dat$RefStrand,"-")] <-  snp.dat$alle1.revcomp[is.element(snp.dat$RefStrand,"-")]
  snp.dat$plus.B[is.element(snp.dat$RefStrand,"-")] <-  snp.dat$alle2.revcomp[is.element(snp.dat$RefStrand,"-")]
  
  ## Add allele definitions to the Allele mappings data frame
  
  ## Map design alleles
  map$alle.design[map$alle.AB=="A"] <- snp.dat$design.A[match(map$snp[map$alle.AB=="A"], snp.dat$Name)]
  map$alle.design[map$alle.AB=="B"] <- snp.dat$design.B[match(map$snp[map$alle.AB=="B"], snp.dat$Name)]
  
  ## Map top alleles
  map$alle.top[map$alle.AB=="A"] <- snp.dat$top.A[match(map$snp[map$alle.AB=="A"], snp.dat$Name)]
  map$alle.top[map$alle.AB=="B"] <- snp.dat$top.B[match(map$snp[map$alle.AB=="B"], snp.dat$Name)]
  
  ## Map fwd alleles
  map$alle.fwd[map$alle.AB=="A"] <- snp.dat$fwd.A[match(map$snp[map$alle.AB=="A"], snp.dat$Name)]
  map$alle.fwd[map$alle.AB=="B"] <- snp.dat$fwd.B[match(map$snp[map$alle.AB=="B"], snp.dat$Name)]
  
  ## Map + alleles
  map$alle.plus[map$alle.AB=="A"] <- snp.dat$plus.A[match(map$snp[map$alle.AB=="A"], snp.dat$Name)]
  map$alle.plus[map$alle.AB=="B"] <- snp.dat$plus.B[match(map$snp[map$alle.AB=="B"], snp.dat$Name)]
  
  dim(map);head(map)

  ## Where design alleles are "NA," change all mapping columns to NA (for cnv probes, e.g.)
  na.snps <- snp.dat$Name[is.element(snp.dat$SNP, c("[N/A]","NA"))]
  cat("\tTotal number for probes with NA alleles (probably CNVs): ", length(na.snps), "\n")
  map$alle.top[is.element(map$snp, na.snps)] <- NA
  map$alle.fwd[is.element(map$snp, na.snps)] <- NA
  map$alle.plus[is.element(map$snp, na.snps)] <- NA
  map$alle.design[is.element(map$snp, na.snps)] <- NA
  
  ## make sure map file is sorted by original order (same as input annotation), delete "order" column, and save to R object
  map.final <- map[order(map$order),-2]
  
  cat("\tMapping file created; do some spot checking with original annotation to confirm mappings!\n")
  
  return(map.final)
  }


make.allele.annotation <- function(map, alleles=c("top", "design", "fwd", "plus")) {
  stopifnot(all(names(map) %in% c("snp", "alle.AB", "alle.design", "alle.top", "alle.fwd", "alle.plus")))
  
  # select which alleles we want
  alleles <- match.arg(alleles)
  col <- paste("alle", alleles, sep=".")

  snp <- map$snp[c(TRUE, FALSE)]
  alleleA <- map[map$alle.AB %in% "A", col]
  alleleB <- map[map$alle.AB %in% "B", col]

  # indels
  indelA <- alleleA %in% "-"
  alleleA[indelA] <- "D"
  alleleB[indelA] <- "I"
  indelB <- alleleB %in% "-"
  alleleA[indelB] <- "I"
  alleleB[indelB] <- "D"
  
  snp.dat <- data.frame(snp, alleleA, alleleB, stringsAsFactors=FALSE)
  names(snp.dat) <- c("snp", paste(c("alleleA", "alleleB"), alleles, sep="."))
  return(snp.dat)
}
