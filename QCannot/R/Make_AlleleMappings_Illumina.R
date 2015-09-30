## Function to take an Illumina SNP annotation data file and make an Allele Mappings table
## Sarah Nelson, UWGCC, May 25, 2011

## Updated March 7, 2012, to allow for indel records where IlmnStrand != SourceStrand
## Updated January 31, 2013 to allow for non-verbose indel mapping (where SourceStrand is not available, print out "I" and "D," respectively, to represent insertion and deletion allele
## Updated October 10, 2013, to allow for Illumina manifests files lacking the "RefStrand" column
## Updated March 22, 2014 to allow for one indel record
## Updated July 10, 2014 to rename some variables, re-organize some code (Tin Louie)
## Updated July 14, 2014 to output indels vcf file as side effect (Tin Louie)
## Updated Aug 1, 2014 to fix bugs/enhance previous update
## Updated Sep 7, 2014 to fix bugs/enhance indels.vcfout feature
## Updated Oct 14, 2014 to report missing columns early on in script (S Nelson)
## Updated Oct 24, 2014 better message about the number of indels with chr0 chrMT pos0
## Updated Jun 18, 2015 handle NA in Refstrand column of input, enhance indels.vcfout feature 
## Updated Sep 29, 2015 preface Hsapiens with BSgenome.Hsapiens.UCSC.hg19:: and preface getSeq, DNAStringSet, reverseComplement with Biostrings::
##
################################# CONTENT
## Example allele mappings table:
#                 snp alle.AB alle.design alle.top alle.fwd alle.plus
#  rs1000000       A           T        A        T         A
#  rs1000000       B           C        G        C         G
#  rs1000002       A           A        A        A         T
#  rs1000002       B           G        G        G         C

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
#  options(stringsAsFactors = FALSE)
#  read SNP Illumina annotation file (i.e., "HumanOmni2.5-4v1_D.csv"), with following required fields:
#  "IlmnID", "Name", "IlmnStrand", "SNP", "SourceSeq"         
#  optional fields:
#  "SourceStrand", "RefStrand"

#  args:
#  indels.verbose - defaults to TRUE, will parse full nucleotide sequence of alleles from the SourceSeq column  
#                   when FALSE, prints I" and "D" characters, respectively, to represent insertion and deletion allele
#               
#  indels.vcfout - default value FALSE; when TRUE, will create indels vcf file
#                  REF and ALT columns can be used to replace Ds and Is in .bim files
#                  no entry for ambiguous indels (could be either insertion or deletion, given the Illumina file)
#                  indels not yet in standardized, left-aligned format 
#                  (which could be the format used by 1000Genomes; 
#                  see Broad's GATK LeftAlignAndTrimVariants option)
#
#  indels.vcfout.filename - default value "indels.needLeftAlign.vcf"
#                           disregard warning "appending column names to file"
#                           

#  returns:
#  data frame object with columns "snp", "alle.AB", "alle.design", "alle.top", "alle.fwd", "alle.plus"
#  *Note "alle.plus" will only be returned if the snp.dat input contains "RefStrand" column

## #Example of using:
## source("/projects/geneva/geneva_sata/GCC_code/Imputation/R_functions/Make_AlleleMappings_Illumina.R")
## source("/projects/geneva/geneva_sata/sarahcn/svn_revisions/trunk/QCpipeline/R/Make_AlleleMappings_Illumina.R")

## ## Load in Illumina annotation
##  column.select=rep("NULL",times=21)
##  column.select[c(1:4,9:11,16:18,21)] <- NA  ## only read in select columns    
##  snp.dat <- read.csv(file="/projects/geneva/gcc-fs2/SNP_annotation/Illumina/HumanOmniExpress_12v1_H/HumanOmniExpress-12v1_H.csv", skip=7,colClasses=column.select,nrows=7000)   ## last rows are for control probes
##  map.final <- make.allele.mappings(snp.dat)
## 	 write.csv(map.final, file="/projects/geneva/geneva_sata/SNP_annotation/Illumina/HumanOmni2.5_4v1/testfn.csv", row.names=FALSE,quote=FALSE))

################################# FUNCTIONS

# getSequence() is a wrapper for getSeq (Biostrings)
# and is invoked by make.allele.mappings(indels.vcfout=TRUE) 
# 
# a$var is the sequence within square brackets [ ], 
# flanked by a$seq1 and a$seq2, 
# from column SourceSeq in the Illumina file
#
# nchar(a$var)  == a$var.n
# nchar(a$seq1) == a$seq1.n
# nchar(a$seq2) == a$seq2.n
#
# return value is an hg19 sequence of length equal to:
# seq1 + seq2       if insertion, or
# seq1 + var + seq2 if deletion
#
# return value *should* match seq1 + seq2 (or seq1 + var + seq2)
# if caller uses column MapInfo from Illumina file
# to correctly guess (without using BLAT) where seq1 and seq2 are located in hg19;
# this guess is represented by a$posI (or a$posD), see I and D below
#
# seq1   var   seq2
# AAAAAA TT    GGGGGG
#      I
#        D
#
getSequence <- function(a, insertion=FALSE) {
	if (insertion) {
		toupper(as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, a$chr.name, a$posI - a$seq1.n +1, a$posI + a$seq2.n)))
	} else {
		toupper(as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, a$chr.name, a$posD - a$seq1.n,    a$posD + a$var.n + a$seq2.n - 1)))
	}
}


# nchar.suffix() is invoked by make.allele.mappings(indels.vcfout=TRUE) 
# to adjust column MapInfo from Illumina file, while attempting to locate seq1 and seq2
# 
# if a$seq1 is at least as long as a$var, then a$var may appear one or more times 
# at the end of a$seq1; this function calculates number of chars of this "suffix"
# by comparing length of a$seq1 (= a$seq1.n), versus the length sans suffix
#   
# the goal is to "right-align" MapInfo values (values that 
# have been left-aligned because the indel appears in a "repeat region")
nchar.suffix <- function(a) {
	ifelse(a$seq1.n >= a$var.n, mapply(function(x,y,z) x - nchar(sub( paste0("(",  y, ")+$"), "", z )), a$seq1.n, a$var, a$seq1), 0)
}


# equalACTGN() is invoked by equal.hg19()
# to test for string equality, 
# skipping over any non-ACGT characters (e.g. N, which is IUPAC wildcard) appearing in arg1
#
equalACGTN <- function(flank, hg19) {
	ACGT = c("A", "C", "G", "T")
	return(0 == mapply(function(c1, c2) sum(c1 != c2 & c1 %in% ACGT), strsplit(flank, ''), strsplit(hg19, '')))
	
	# if side (arg3) == 1, flank is to the "left" (upstream) of the indel,
	# otherwise flank is to the "right" (downstream);
	# test only LEN bases adjacent to the indel
	# if (side == 1) {
		# lastN  = nchar(flank)
		# startN = lastN -LEN +1
		# return(0 == mapply(function(c1, c2) sum(c1 != c2 & c1 %in% ACGT), strsplit(substr(flank, startN, lastN), ''), strsplit(substr(hg19, startN, lastN), '')))
	# } else {
		# return(0 == mapply(function(c1, c2) sum(c1 != c2 & c1 %in% ACGT), strsplit(substr(flank, 1, LEN), ''), strsplit(substr(hg19, 1, LEN), '')))
	# }
	
}

# equal.hg19() is invoked by make.allele.mappings(indels.vcfout=TRUE) 
#
# a$var is the sequence within square brackets [ ], 
# flanked by a$seq1 and a$seq2, 
# from column SourceSeq in the Illumina file
#
# pre-req: 
# nchar(a$var)  == a$var.n
# nchar(a$seq1) == a$seq1.n
# nchar(a$seq2) == a$seq2.n
# a$seq1.n >= LEN (a constant defined below; it's smaller than the Illumina 50-mer probe)
# a$seq2.n >= LEN
#
# if insertion is TRUE:
# disregard a$var, and
# return TRUE if one flank matches the hg19 sequence stored in a$hg19I
#                and the other flank matches at least the LEN characters adjacent to the insertion  
#
# if insertion is FALSE:
# return TRUE if the same conditions above are met, 
#                and a$var matches a$hg19D exactly
#
equal.hg19 <- function(a, insertion=FALSE) {
	LEN = 30
	probeLEN = 50
	stopifnot(all(a$seq1.n >= LEN))
	stopifnot(all(a$seq2.n >= LEN))
	
	if (insertion) {
		seq2.near = equalACGTN(substr(a$seq2, 1, LEN), substr(a$hg19I, a$seq1.n +1, a$seq1.n + LEN))
		seq2.far  = equalACGTN(substr(a$seq2, LEN +1, probeLEN), substr(a$hg19I, a$seq1.n + LEN +1, a$seq1.n + probeLEN))
		seq1.far  = equalACGTN(substr(a$seq1, a$seq1.n -probeLEN +1, a$seq1.n -LEN), substr(a$hg19I, a$seq1.n -probeLEN +1, a$seq1.n -LEN)) 
		seq1.near = equalACGTN(substr(a$seq1, a$seq1.n -LEN +1, a$seq1.n), substr(a$hg19I, a$seq1.n -LEN +1, a$seq1.n)) 
		
		((a$seq1.n >= probeLEN) & seq1.far & seq1.near & seq2.near) | (seq1.near & seq2.near & seq2.far & (a$seq2.n >= probeLEN))
		# equalACGTN(a$seq1, substr(a$hg19I, 1, a$seq1.n), side=1) & 
		# equalACGTN(a$seq2, substr(a$hg19I, a$seq1.n +1, a$seq1.n + a$seq2.n), side=2)
	} else {
		seq2.near = equalACGTN(substr(a$seq2, 1, LEN), substr(a$hg19D, a$seq1.n +a$var.n +1, a$seq1.n +a$var.n +LEN))
		seq2.far  = equalACGTN(substr(a$seq2, LEN +1, probeLEN), substr(a$hg19D, a$seq1.n + a$var.n +LEN +1, a$seq1.n + a$var.n +probeLEN))
		seq1.far  = equalACGTN(substr(a$seq1, a$seq1.n -probeLEN +1, a$seq1.n -LEN), substr(a$hg19D, a$seq1.n -probeLEN +1, a$seq1.n -LEN)) 
		seq1.near = equalACGTN(substr(a$seq1, a$seq1.n -LEN +1, a$seq1.n), substr(a$hg19D, a$seq1.n -LEN +1, a$seq1.n)) 
	
		(a$var == substr(a$hg19D, a$seq1.n +1, a$seq1.n +a$var.n)) & 
		(((a$seq1.n >= probeLEN) & seq1.far & seq1.near & seq2.near) | (seq1.near & seq2.near & seq2.far & (a$seq2.n >= probeLEN)))
		         
	
		# equalACGTN(a$seq1, substr(a$hg19D, 1, a$seq1.n), side=1) & 
		# (a$var == substr(a$hg19D, a$seq1.n +1, a$seq1.n +a$var.n)) & 
		# equalACGTN(a$seq2, substr(a$hg19D, a$seq1.n +a$var.n +1, a$seq1.n +a$var.n +a$seq2.n), side=2)
	}
	
}

# approx.hg19() is more complex than equal.hg19()
# returns a numeric score instead of TRUE/FALSE, where 
# higher scores = more similarity to hg19
# (except when one.sided (arg3) is TRUE, see below)
# 
# insertions:
# each flank is scored zero if flank is too dissimilar to hg19, 
# otherwise score = [length of flank] minus [edit distance between flank & hg19]
#
# deletions:
# score is adjusted by:
# [length of deletion] minus 6 * [number of mismatches between deletion & hg19]
#
# n.b. the deletion, as it appears in Illumina file, is expected to almost always 
# match something in hg19
#
# one.sided:
# used by pseudoBLAT() to see if seq1 or seq2 has score zero
#
approx.hg19 <- function(b, insertion=FALSE, one.sided=FALSE) {
	var.adjustment = 0
	distance1 = NULL
	distance2 = NULL
	
	if (insertion) {
		distance1 = mapply(function(x, y) adist(x,y)[1,1], b$seq1, substr(b$hg19I, 1, b$seq1.n))
		distance2 = mapply(function(x, y) adist(x,y)[1,1], b$seq2, substr(b$hg19I, b$seq1.n +1, b$seq1.n + b$seq2.n))

		#distance1 = mapply(function(x, y) adist(x,y)[1,1], substr(b$seq1, b$seq1.n - LEN + 1, b$seq1.n), substr(b$hg19I, pmax(1, b$seq1.n - LEN +1), b$seq1.n))
		#distance2 = mapply(function(x, y) adist(x,y)[1,1], substr(b$seq2, 1, LEN), substr(b$hg19I, b$seq1.n +1, pmin(b$seq1.n + LEN, b$seq1.n + b$seq2.n)))
		
		
	} else {
		distance1 = mapply(function(x, y) adist(x,y)[1,1], b$seq1, substr(b$hg19D, 1, b$seq1.n))
		distance2 = mapply(function(x, y) adist(x,y)[1,1], b$seq2, substr(b$hg19D, b$seq1.n +1 +b$var.n, b$seq1.n + b$var.n + b$seq2.n))
		var.adjustment = b$var.n  - (6 * mapply(function(c1, c2) sum(c1 != c2), strsplit(b$var, ''), strsplit(substr(b$hg19D, b$seq1.n +1, b$seq1.n +b$var.n), '')))

		# distance1 = mapply(function(x, y) adist(x,y)[1,1], substr(b$seq1, b$seq1.n - LEN + 1, b$seq1.n), substr(b$hg19D, pmax(1, b$seq1.n - LEN +1), b$seq1.n))
		# distance2 = mapply(function(x, y) adist(x,y)[1,1], substr(b$seq2, 1, LEN), substr(b$hg19D, b$seq1.n +1 +b$var.n, pmin(b$seq1.n + b$var.n + LEN, b$seq1.n + b$var.n + b$seq2.n)))
		
	}
	
	# score is 0 if too dissimilar (distance > 6)
	seq1.score = ifelse(distance1 > 6, 0, b$seq1.n - distance1)
	seq2.score = ifelse(distance2 > 6, 0, b$seq2.n - distance2)
	
	if (one.sided) {
		# return 1 if seq1.score is zero, 
		#        2 if seq2.score is zero, 
		#        0 if neither is zero, or both are zero
		ifelse((seq1.score == 0) & (seq2.score > 0), 
		        1, 
		        ifelse((seq1.score > 0) & (seq2.score == 0), 2, 0) ) 
		
	} else {
		return(seq1.score + seq2.score + var.adjustment)
	}
	
}


# pseudoBLAT() is like approx.hg19(), except the focus is on deletions only
# return value has same number of rows (i.e. same number of indels) as b (arg1) 
# but only 3 columns:
# 1) scoreD
# 2) hg19D
# 3) posD

# 1) scoreD will be same as b$scoreD.old (the score of the best deletion scenario so far)
#    unless this function can obtain a higher score by finding an exact match to 
#    seq1 or seq2 in hg19 
#
# 2) hg19D may be the same as b$hg19D.old; if it is different, then:
#    a) hg19D will be *longer* than seq1 + var + seq2, implying that the deletion is 
#       longer than what is shown in column ScoreSeq from Illumina file
#    b) one end of hg19D will be an exact match to seq1 or seq2
#    c) scoreD will be higher than scoreD.old
#
# 3) posD may be the same as b$posD.old; in any event, it can be used 
#    in conjunction with hg19D to print the variant to the vcf output file,
#    (assuming the variant is not ruled an insertion, of course)
#
#
pseudoBLAT <- function(b) {
	# score of the best deletion scenario so far, and its associated variables
	b$posD   = b$posD.old
	b$hg19D  = b$hg19D.old
	b$scoreD = b$scoreD.old
	
	# process only those variants where we haven't found a good match
	# for seq1 or seq2 in hg19
	b$one.sided = approx.hg19(b, one.sided=TRUE)
	b$originalOrder = 1:nrow(b)
	upstream   = b[b$one.sided == 1,]
	downstream = b[b$one.sided == 2,]
	neither    = b[b$one.sided == 0,]
	
	# avoid R CMD check notes
	chr.name = posD = var.n = seq1 = hg19L = seq1 = seq2 = hit = hg19D = seq1.n = seq2.n = scoreD = NULL
	
	LEN = 9999
	if (nrow(downstream) > 0) {
		# grab a long hg19 sequence (LEN bases) upstream or downstream of deletion
		# (deletion begins at posD according to best deletion scenario so far) 
		downstream$hg19L = with(downstream, toupper(as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chr.name, posD + var.n, posD + var.n + LEN -1))))

		# look for exact match of seq1 or seq2 in long hg19 sequence
		# (if there're multiple matches, get the closest one to the deletion)
		downstream$hit = with(downstream, mapply(function(x, y) regexpr(x, y, fixed=T), seq2, hg19L))
		
		# if no match, then hg19D remains the same as hg19D.old, 
		# else extend hg19D by that portion of hg19L up to the hit 
		downstream$hg19D = with(downstream, ifelse(hit == -1, hg19D, paste0(substr(hg19D, 1, seq1.n + var.n),
																			substr(hg19L, 1, hit + seq2.n - 1)) ))
	

		# if no match, then scoreD remains the same as scoreD.old
		# else increase the score by the length of seq1 or seq2 
		# (see the scoring system in approx.hg19();
		#  remember that processed variants are one-sided, i.e. zero contribution to the score from seq1 or seq2)
		downstream$scoreD = with(downstream, ifelse(hit == -1, scoreD, scoreD + seq2.n))
	}
	
	if (nrow(upstream) > 0) {
		upstream$hg19L = with(upstream, toupper(as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chr.name, posD - LEN,   posD -1))))
		upstream$hit   = with(upstream, mapply(function(x, y) rev(gregexpr(x, y, fixed=T)[[1]])[1], seq1, hg19L))
		upstream$hg19D = with(upstream, ifelse(hit == -1, hg19D, paste0(substr(hg19L, hit, LEN), 
	                                                                    substr(hg19D, seq1.n +1, seq1.n + var.n + seq2.n)) ))
		upstream$scoreD = with(upstream, ifelse(hit == -1, scoreD, scoreD + seq1.n))

		# if no match, then posD remains the same as posD.old
		# else for upstream match to seq1, shift posD upstream
		# (deletion now extends all the way up to end of seq1)
		# else for downstream match to seq2, posD remains the same as posD.old 
		# (because even though the deletion has been extended, its starting location has not changed)
		upstream$posD = with(upstream, ifelse(hit == -1, posD, posD - LEN + hit + seq1.n -1))
	}
	
	
	# combine processed variants with un-processed variants, restore original ordering of rows
	neither    =    neither[,c("posD", "hg19D", "scoreD", "originalOrder")]
	upstream   =   upstream[,c("posD", "hg19D", "scoreD", "originalOrder")]
	downstream = downstream[,c("posD", "hg19D", "scoreD", "originalOrder")]
	b = rbind(neither, upstream, downstream)
	b = b[order(b$originalOrder),]
	return(b[,c("posD", "hg19D", "scoreD")])
}




make.allele.mappings <- function(snp.dat, indels.verbose = TRUE, indels.vcfout = FALSE, indels.vcfout.filename = "indels.needLeftAlign.vcf") {
  options(stringsAsFactors = FALSE)

  # avoid R CMD check notes
  scoreD = scoreI = scoreD.old = scoreI.old = hg19D.old = hg19I.old = posD.old = posI.old = NULL
  
  if (indels.vcfout && (!requireNamespace("Biostrings", quietly=T) || !requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly=T)) ) {
	indels.vcfout <- FALSE
	warning("\tno vcf file will be created")
  }
  
  ## Report total number of SNPs and indels
  cat("\tTotal number of probes:", length(snp.dat$Name),"\n")
  cat("\tTotal number of indels:", length(snp.dat$Name[is.element(snp.dat$SNP, c("[D/I]","[I/D]"))]),"\n")
  
  ## Check that SNP name is unique and appropriate to use as the identifier in allele mappings
  if (any(duplicated(snp.dat$Name))) {
  	stop(sum(duplicated(snp.dat$Name)), " SNP names are duplicated; need alternate identifier to make Allele Mappings file")
  }

  ## Check for presence of minimum required columns
  cols.needed <- c("IlmnID","Name","SNP","IlmnStrand")
  cols.missing <- setdiff(cols.needed, names(snp.dat))
  if(length(cols.missing > 0)) {
    stop("SNP input is missing columns(s) ",paste(cols.missing, collapse=" "),"\n")
  }

  # If indels.vcfout=TRUE, also need Chr and Mapinfo
  cols.vcf <- c("Chr","MapInfo")
  cols.missing <- setdiff(cols.vcf, names(snp.dat))
  if(length(cols.missing > 0) & indels.vcfout ) {
    stop("Indel VCF requested, but SNP input is missing columns(s) ",paste(cols.missing,collapse=" "),"\n")
  }

  ## detect if snp.dat contains "RefStrand" - if not, print warning message that plus(+) alleles will not be written out
  print.plus <- TRUE
  if(!is.element("RefStrand", names(snp.dat))) {
       print.plus <- FALSE
       warning("\tSNP manifest lacks RefStrand column; plus(+) alleles NOT written to output")
  }
  
  ## set up table template - 2 rows per SNP, first "A", second "B"
  map <- as.data.frame(c(snp.dat$Name, snp.dat$Name))
  names(map) <- "snp"
  snp.dat$order <- 1:length(snp.dat$Name)
  
  ## re-order data frame to match SNP annotation
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
  
	if (length(indels) > 0) {
		cat("\tCreating mappings for insertion/deletion probes\n")
        if (indels.verbose || indels.vcfout) {
			# extract 4 portions from each indel's SourceSeq:
			# sourceseq.start = portion before the brackets []
			# sourceseq.end   = portion after the brackets []
			# left            = inside the brackets & before the slash /
			# right           = inside the brackets & after the slash /
            tmp <- strsplit(indels.dat$SourceSeq, "[", fixed = TRUE)
			sourceseq.start <- toupper(sapply(tmp, "[[", 1))
            tmp <- sapply(tmp, "[[", 2)
			
			tmp <- strsplit(tmp, "]", fixed = TRUE)
            sourceseq.end <- toupper(sapply(tmp, "[[", 2))
            tmp <- sapply(tmp, "[[", 1)
			
            tmp <- strsplit(tmp, "/", fixed = TRUE)
            left  <- toupper(sapply(tmp, "[[", 1))
            right <- toupper(sapply(tmp, "[[", 2))

			# gather info about each indel:
			# var = the ACGTs inside the brackets
			# strand.diff = TRUE if IlmnStrand is different from SourceStrand
            a <- as.data.frame(cbind(left, right))
            a$var <- ifelse(a$left == "-", a$right, a$left)
            tmp <- Biostrings::DNAStringSet(a$var)
            a$var.revcomp <- as.character(Biostrings::reverseComplement(tmp))
            a$Name <- indels
            a$DI <- indels.dat$SNP == "[D/I]"
            a$strand.diff <- indels.dat$IlmnStrand != indels.dat$SourceStrand

			# design.A and .B depend on DI and strand.diff
            indels.dat$design.A[a$DI] <- "-"
            indels.dat$design.B[a$DI & !a$strand.diff] <- a$var[match(indels.dat$Name[a$DI & !a$strand.diff], a$Name)]
            indels.dat$design.B[a$DI & a$strand.diff] <- a$var.revcomp[match(indels.dat$Name[a$DI & a$strand.diff], a$Name)]

            indels.dat$design.B[!a$DI] <- "-"
            indels.dat$design.A[!a$DI & !a$strand.diff] <- a$var[match(indels.dat$Name[!a$DI & !a$strand.diff], a$Name)]
            indels.dat$design.A[!a$DI & a$strand.diff] <- a$var.revcomp[match(indels.dat$Name[!a$DI & a$strand.diff], a$Name)]

			# top.A and .B depend on IlmnStrand and design.A and .B (or revcomp)
			# PLUS is equivalent to TOP, MINUS is equivalent to BOT
            tmp <- Biostrings::DNAStringSet(indels.dat$design.A)
            indels.dat$alle1.revcomp <- as.character(Biostrings::reverseComplement(tmp))
            tmp <- Biostrings::DNAStringSet(indels.dat$design.B)
            indels.dat$alle2.revcomp <- as.character(Biostrings::reverseComplement(tmp))

            indels.dat$top.A[is.element(indels.dat$IlmnStrand,
                c("P", "PLUS"))] <- indels.dat$design.A[is.element(indels.dat$IlmnStrand, c("P", "PLUS"))]
            indels.dat$top.B[is.element(indels.dat$IlmnStrand,
                c("P", "PLUS"))] <- indels.dat$design.B[is.element(indels.dat$IlmnStrand, c("P", "PLUS"))]
			
            indels.dat$top.A[is.element(indels.dat$IlmnStrand,
                c("M", "MINUS"))] <- indels.dat$alle1.revcomp[is.element(indels.dat$IlmnStrand, c("M", "MINUS"))]
            indels.dat$top.B[is.element(indels.dat$IlmnStrand,
                c("M", "MINUS"))] <- indels.dat$alle2.revcomp[is.element(indels.dat$IlmnStrand, c("M", "MINUS"))]
				
			### make vcf file ###
			if (indels.vcfout) {
				a$var.original    = a$var
				a$sourceseq.start = sourceseq.start
				a$sourceseq.end   = sourceseq.end
				
				# choose between using a$var or a$var.revcomp, but try both if warranted
				# (run code twice using for loop (below), unless break at bottom of loop)
				if (print.plus) { # is there a RefStrand column ? 
					a$shouldflip = xor(a$strand.diff, indels.dat$RefStrand == "-")

					# ensure shouldflip is not NA
					# (strand.diff could be NA, RefStrand could be NA)
					a$shouldflip = ifelse(is.na(a$shouldflip), FALSE, a$shouldflip)

				} else {
					a$shouldflip = FALSE
				}
				
				a$chr     = indels.dat$Chr
				a$MapInfo = indels.dat$MapInfo
				# cannot evaluate chr0 or chrMT or position 0, so drop them
				# (hg19 MT inaccessible via R library currently used) 
				num.indels = nrow(a)
				a = a[!is.na(a$chr) & !is.na(a$MapInfo) & is.element(a$chr, c(1:22, "X", "Y", "XY")) & a$MapInfo != 0,] 
				num.eliminated = num.indels - nrow(a)
				if (num.eliminated > 0) {
					message("\tn.b. ", num.eliminated, " indels won't appear in output due to chrMT chr0 position0")
				}

				a$chr      = ifelse(a$chr == "XY", "X", a$chr)
				a$chr.name = paste0("chr", a$chr)

				iteration1.results = NULL
				for (iteration in 1:2) {
					if (iteration == 2) {
						message("\tre-try ", nrow(a), " indels by using reverse complement")
					}
					
					# gather more info per indel
					# seq1      var   seq2
					# A.....A   TT    G.....G
					a$var  = ifelse(a$shouldflip, a$var.revcomp, a$var.original)
					a$seq1 = ifelse(a$shouldflip, as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(a$sourceseq.end))),   a$sourceseq.start)
					a$seq2 = ifelse(a$shouldflip, as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(a$sourceseq.start))), a$sourceseq.end)
					
					a$var.n  = nchar(a$var)
					a$seq1.n = nchar(a$seq1)
					a$seq2.n = nchar(a$seq2)
					

					# do NOT trim flanking seqs (let equal.hg19() compare up to LEN chars) 
					# LEN = ?
					# a$seq1   = ifelse(a$seq1.n > LEN, substr(a$seq1, a$seq1.n -LEN +1, a$seq1.n), a$seq1)
					# a$seq2   = ifelse(a$seq2.n > LEN, substr(a$seq2, 1, LEN), a$seq2)
					# a$seq1.n = ifelse(a$seq1.n > LEN, LEN, a$seq1.n)
					# a$seq2.n = ifelse(a$seq2.n > LEN, LEN, a$seq2.n)
					
					
					# column MapInfo of Illumina file contains genomic position 
					# which usually corresponds to * or @ below
					# so that BLAT is not needed to find seq1/var/seq2 in hg19
					#
					# seq1   var   seq2
					# AAAAAA TT    GGGGGG
					#      *
					#        @
					#
					# position = * if variant is an insertion (e.g. 1KG_2_71060821)
					#          = @ if               deletion  (e.g. exm-IND2-70914329)
					# 
					# MapInfo is not always at * or @
					# so we test different scenarios below
					# (hope the adjusted position lands at * or @)
					
					# initialize
					a$posD  = 0
					a$hg19D = "ACGT"
					a$del   = FALSE # mark variants as possible deletions using a$del = TRUE
									
					# loop through scenarios
					# skip variants already determined as possible deletions 
					scenarios = cbind(a$MapInfo,  a$MapInfo +1,  a$MapInfo - a$var.n +1,  a$MapInfo - a$var.n +2,  
									  a$MapInfo +nchar.suffix(a), 
									  a$MapInfo +nchar.suffix(a) +1, 
									  a$MapInfo +nchar.suffix(a) -a$var.n +1, 
									  a$MapInfo +nchar.suffix(a) -a$var.n +2)
					for (i in 1:ncol(scenarios)) {
						a$posD  = ifelse(a$del, a$posD,  scenarios[,i])
						a$hg19D = ifelse(a$del, a$hg19D, getSequence(a))
						a$del   = ifelse(a$del, TRUE,    equal.hg19(a))
					}

					# initialize for insertion scenarios:
					a$posI  = 0
					a$hg19I = "ACGT"
					a$ins   = FALSE
					scenarios = cbind(a$MapInfo, a$MapInfo -1, a$MapInfo +nchar.suffix(a),  a$MapInfo -1 +nchar.suffix(a))
					for (i in 1:ncol(scenarios)) {
						a$posI  = ifelse(a$ins, a$posI,  scenarios[,i])
						a$hg19I = ifelse(a$ins, a$hg19I, getSequence(a, insertion=TRUE))
						a$ins   = ifelse(a$ins, TRUE,    equal.hg19(a, insertion=TRUE))
					}
					
					# if any variant has not yet been tagged as insertion or deletion
					# (because of equal.hg19() test results)
					# then go through the above scenarios again (plus one new scenario), 
					# but instead of equal.hg19(), use approx.hg19() which returns 
					# a numeric score per scenario, then pick scenario with highest score
					a$approx.hg = "." # to be put into INFO column of vcf  
					b = a[!a$ins & !a$del,]
					if (nrow(b) > 0) {
						message("\tnumber of indels that will undergo approximate matching is ", nrow(b))
						b$approx.hg = "APPROX"
						
						b$scoreD     = 0
						b$scoreD.old = -10 # dummy value (best score so far)
						b$posD.old = 0
						b$hg19D.old = "ACGT"
						
						scenarios = cbind(b$MapInfo, b$MapInfo +1,  b$MapInfo - b$var.n +1,  b$MapInfo - b$var.n +2,  
										  b$MapInfo +nchar.suffix(b), 
										  b$MapInfo +nchar.suffix(b) +1,
										  b$MapInfo +nchar.suffix(b) -b$var.n +1,
										  b$MapInfo +nchar.suffix(b) -b$var.n +2)
						for (i in 1:ncol(scenarios)) {
							b$posD       = scenarios[,i]
							b$hg19D      = getSequence(b)
							b$scoreD     = approx.hg19(b) 
						
							# compare scenarios, keep the better score & its associated variables in .old
							b$posD.old    = ifelse(b$scoreD.old >= b$scoreD, b$posD.old,   b$posD)
							b$hg19D.old   = ifelse(b$scoreD.old >= b$scoreD, b$hg19D.old,  b$hg19D)
							b$scoreD.old  = ifelse(b$scoreD.old >= b$scoreD, b$scoreD.old, b$scoreD)
						}
						
						# new deletion scenario 
						results = pseudoBLAT(b) 
						
						# compare & keep the better score 
						# n.b. pseudoBLAT results imply a new b$var without changing b$var or b$var.n, 
						#      so don't run more deletion scenarios
						b$posD    = ifelse(b$scoreD.old >= results$scoreD, b$posD.old,   results$posD)
						b$hg19D   = ifelse(b$scoreD.old >= results$scoreD, b$hg19D.old,  results$hg19D)
						b$scoreD  = ifelse(b$scoreD.old >= results$scoreD, b$scoreD.old, results$scoreD)
						

						# insertion scenarios
						b$scoreI     = 0
						b$scoreI.old = -10
						b$posI.old = 0
						b$hg19I.old = "ACGT"

						scenarios = cbind(b$MapInfo, b$MapInfo -1, b$MapInfo +nchar.suffix(b), b$MapInfo -1 +nchar.suffix(b))
						for (i in 1:ncol(scenarios)) {
							b$posI   = scenarios[,i]
							b$hg19I  = getSequence(b, insertion=TRUE)
							b$scoreI = approx.hg19(b, insertion=TRUE)
							
							b$posI.old    = ifelse(b$scoreI.old >= b$scoreI, b$posI.old,   b$posI)
							b$hg19I.old   = ifelse(b$scoreI.old >= b$scoreI, b$hg19I.old,  b$hg19I)
							b$scoreI.old  = ifelse(b$scoreI.old >= b$scoreI, b$scoreI.old, b$scoreI)
						}
						b$posI    = b$posI.old
						b$hg19I   = b$hg19I.old
						b$scoreI  = b$scoreI.old

						
						# compare best insertion scenario vs best deletion scenario
						# also, the deletion score should be high enough to imply a perfect var match plus one matching flank
						# while the insertion score should be high enough to imply one matching flank plus one ok flank
						# n.b. both b$del and b$ins may remain FALSE 
						b$del = (b$scoreD > b$scoreI) & (b$scoreD >= (pmin(b$seq1.n, b$seq2.n) + b$var.n))
						b$ins = (b$scoreI > b$scoreD) & (b$scoreI >  (pmin(b$seq1.n, b$seq2.n) + 10))
						# b$del  = ((b$scoreD > (b$scoreI + 0.02)) & (b$scoreD > 0.6)) | ((b$scoreD > b$scoreI) & (b$scoreD > 0.8))
						# b$ins  = ((b$scoreI > (b$scoreD + 0.02)) & (b$scoreI > 0.6)) | ((b$scoreI > b$scoreD) & (b$scoreI > 0.8))
						
						# remove extraneous columns
						# merge b into a
						b = subset(b, select = -c(scoreD, scoreI, scoreD.old, scoreI.old, hg19D.old, hg19I.old, posD.old, posI.old))
						a = a[a$ins | a$del,]
						a = rbind(a, b)
					}
				
					# is iteration 2 warranted?
					# not if all variants have been tagged as insertion, deletion, or both
					if (iteration == 1) {
						b = a[!a$ins & !a$del,]
						if (nrow(b) > 0) { 
							iteration1.results = a[a$ins | a$del,]
							a = b
							a$shouldflip = ! a$shouldflip
						} else {
							break 
						}
					} else { # else this is iteration 2
						a = rbind(a, iteration1.results)
					}
					
				} # end of for loop with iteration 1 and 2
				
				# keep only indels that are unambiguously an insertion or deletion 
				# i.e. toss indels that could be either, and toss indels that appear to be neither
				num.attempted = nrow(a)
				a = a[xor(a$ins, a$del),]
				num.unambiguous = nrow(a)
				if (num.unambiguous > 0) {
					if (num.attempted != num.unambiguous) {
						warning("\tunable to classify ", num.attempted - num.unambiguous, " indels (will not appear in vcf)")
					}
					# define pos, ref, alt columns of vcf output
					a$pos = ifelse(a$ins, a$posI, a$posD -1) # reconcile vcf requirements & posI or posD (see * @ above)
					
					# n.b. pseudoBLAT results, if any were kept, invalidate var.n, so don't use it
					# a$ref = ifelse(a$ins, substr(a$hg19I, a$seq1.n, a$seq1.n), substr(a$hg19D, a$seq1.n, a$seq1.n +a$var.n))
					a$ref = ifelse(a$ins, substr(a$hg19I, a$seq1.n, a$seq1.n), substr(a$hg19D, a$seq1.n, nchar(a$hg19D) - a$seq2.n))
					
					a$alt = ifelse(a$ins, paste0(a$ref, a$var), substr(a$ref, 1, 1))
					
					# sort by chr, pos 
					vcf = a[,c("chr", "pos", "Name", "ref", "alt", "approx.hg")]
					vcf$chr = ifelse(vcf$chr == "X", "23", vcf$chr)
					vcf$chr = ifelse(vcf$chr == "Y", "24", vcf$chr)
					vcf$chr = as.integer(vcf$chr)
					vcf = vcf[order(vcf$chr, vcf$pos),]
					vcf$chr = as.character(vcf$chr)
					vcf$chr = ifelse(vcf$chr == "23", "X", vcf$chr)
					vcf$chr = ifelse(vcf$chr == "24", "Y", vcf$chr)
					
					# create dummy values to comply with vcf format
					vcf$column6 = 100
					vcf$column7 = "PASS"
					vcf$column8 = vcf$approx.hg
					vcf$approx.hg <- NULL
					vcf$column9 = "GT"
					vcf$column10 = "0/1"
					
					# write header rows of vcf
					# disregard warning "appending column names to file"
					writeLines(c("##fileformat=VCFv4.1", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">", "##INFO=<ID=APPROX,Number=0,Type=Flag,Description=\"underwent approximate matching\">"), indels.vcfout.filename)
					colnames(vcf) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "dummySubject")
					write.table(vcf, indels.vcfout.filename, sep="\t", quote=F, row.names=F, col.names=T, append=T)
				} else {
					warning("\tno vcf written (no indel is unambiguously an insertion or deletion)")
				}
				
			} # end of if indels.vcfout
        } # end of if indels.verbose or indels.vcfout 
		
		if (! indels.verbose) {  
			## write out simpler version of indel mapping
            del.first <- indels.dat[indels.dat$SNP == "[D/I]",]
            ins.first <- indels.dat[indels.dat$SNP == "[I/D]",]
            if (nrow(del.first) > 0) {
                del.first$top.A <- del.first$design.A <- "D"
                del.first$top.B <- del.first$design.B <- "I"
            }
            if (nrow(ins.first) > 0) {
                ins.first$top.A <- ins.first$design.A <- "I"
                ins.first$top.B <- ins.first$design.B <- "D"
            }
            indels.dat <- rbind(del.first, ins.first)
            indels.dat$alle1.revcomp <- indels.dat$design.A
            indels.dat$alle2.revcomp <- indels.dat$design.B
        }

		# temporarily remove indels from snp.dat
		snp.dat <- snp.dat[!is.element(snp.dat$Name, indels),] 
		
	} # end of if length(indels) > 0

  
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
  
  ## now that design and TOP alleles are defined, 
  ## rejoin indels with main SNP annotation, keeping select columns
  cols.keep <- c("IlmnID","Name","IlmnStrand","SNP",
                 "design.A","design.B","alle1.revcomp","alle2.revcomp","top.A","top.B","order")

  message("snp.dat has columns", paste(names(snp.dat),collapse=" "))
  
 ## for debugging - find which of "cols.keep" may be missing from snp.dat
  cols.missing <- setdiff(cols.keep, names(snp.dat))
  if(length(cols.missing > 0)) {
    warning("snp.dat is missing columns(s) ",paste(cols.missing, collapse=" "),"\n")
  }
  
  if (print.plus) cols.keep <- c("RefStrand",cols.keep)
  
  if (length(indels) > 0) {
    comb <- rbind(snp.dat[,cols.keep], indels.dat[,cols.keep])
  } else { 
    comb <- snp.dat[,cols.keep] 
  }
  
  ## restore original order, rename to snp.dat
  snp.dat <- comb[order(comb$order),]
  
  ## get FORWARD alleles: make column to indicate dbSNP orientation of DESIGN strand
  cat("\tDefining FORWARD alleles\n")
  snp.dat$dbSNPStrand.fordesign[is.element(snp.dat$IlmnID,grep("_F_",snp.dat$IlmnID, value=TRUE))] <- "FWD"
  snp.dat$dbSNPStrand.fordesign[is.element(snp.dat$IlmnID,grep("_R_",snp.dat$IlmnID, value=TRUE))] <- "REV"
  
  snp.dat$fwd.A[snp.dat$dbSNPStrand.fordesign=="FWD"] <- snp.dat$design.A[snp.dat$dbSNPStrand.fordesign=="FWD"]
  snp.dat$fwd.B[snp.dat$dbSNPStrand.fordesign=="FWD"] <- snp.dat$design.B[snp.dat$dbSNPStrand.fordesign=="FWD"]
  
  snp.dat$fwd.A[snp.dat$dbSNPStrand.fordesign=="REV"] <- snp.dat$alle1.revcomp[snp.dat$dbSNPStrand.fordesign=="REV"]
  snp.dat$fwd.B[snp.dat$dbSNPStrand.fordesign=="REV"] <- snp.dat$alle2.revcomp[snp.dat$dbSNPStrand.fordesign=="REV"]

  if (print.plus){
      cat("\tDefining PLUS(+) alleles\n")

      ## get + alleles, wrt human genome reference sequence
      ## In Illumina b37, RefStrand is BLAST result of Design Strand
      ## where RefStrand="+", same as DESIGN alleles
      ## where RefStrand="-", reverse complement of DESIGN alleles

	## if RefStrand is NA, plus alleles will be NA
	
      snp.dat$plus.A[is.element(snp.dat$RefStrand,"+")] <- snp.dat$design.A[is.element(snp.dat$RefStrand,"+")]
      snp.dat$plus.B[is.element(snp.dat$RefStrand,"+")] <- snp.dat$design.B[is.element(snp.dat$RefStrand,"+")]

      snp.dat$plus.A[is.element(snp.dat$RefStrand,"-")] <-  snp.dat$alle1.revcomp[is.element(snp.dat$RefStrand,"-")]
      snp.dat$plus.B[is.element(snp.dat$RefStrand,"-")] <-  snp.dat$alle2.revcomp[is.element(snp.dat$RefStrand,"-")]
    }

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

  if(print.plus) {
      ## Map + alleles
      map$alle.plus[map$alle.AB=="A"] <- snp.dat$plus.A[match(map$snp[map$alle.AB=="A"], snp.dat$Name)]
      map$alle.plus[map$alle.AB=="B"] <- snp.dat$plus.B[match(map$snp[map$alle.AB=="B"], snp.dat$Name)]
    }

  
  dim(map);head(map)

  ## Where design alleles are "NA," change all mapping columns to NA (for cnv probes, e.g.)
  na.snps <- snp.dat$Name[is.element(snp.dat$SNP, c("[N/A]","NA"))]
  cat("\tTotal number for probes with NA alleles (probably CNVs): ", length(na.snps), "\n")
  map$alle.top[is.element(map$snp, na.snps)] <- NA
  map$alle.fwd[is.element(map$snp, na.snps)] <- NA
  map$alle.design[is.element(map$snp, na.snps)] <- NA
  if(print.plus) {map$alle.plus[is.element(map$snp, na.snps)] <- NA}
  
  ## make sure map file is sorted by original order (same as input annotation), delete "order" column, and save to R object
  map.final <- map[order(map$order),-2]
  
  cat("\tMapping file created; do some spot checking against original annotation to confirm mappings!\n\n")
  
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
