\name{make.allele.mappings}
\alias{make.allele.mappings}
\alias{make.allele.annotation}
\title{Make allele mapping file}
\description{Function to take Illumina SNP annotation data file and make an Allele Mappings table}
\usage{
make.allele.mappings(snp.dat, indels.verbose=TRUE)

make.allele.annotation(map, alleles=c("top", "design", "fwd", "plus"))
}
\arguments{
  \item{snp.dat}{a data frame made from SNP Illumina
    annotation file (i.e., "HumanOmni2.5-4v1_D.csv"), with the following
    required fields: "IlmnID", "Name", "IlmnStrand", "SNP", and "SourceSeq"; 
	and the following optional fields:"SourceStrand" and "RefStrand"   }
  \item{indels.verbose}{if \code{FALSE}, alleles for indels will be written out as 'I' and 'D,' rather than parsing the full nucleotide sequence of indel alleles from the SourceSeq column}
  \item{map}{a data frame with allele mappings (output of
    \code{make.allele.mappings}), with following fields: "snp",
    "alle.AB", "alle.design", "alle.top", "alle.fwd", "alle.plus"}
  \item{alleles}{character string with type of alleles to return in the
    annotation table ("top", "design", "fwd", or "plus")}
}
\details{
  In the case of indels, the \code{make.allele.mappings} codes A or B as
  "-" and the other allele with the insertion/deletion sequence.
  \code{make.allele.annotation} codes indels as "D" and "I" (useful for
  making PLINK files).
  
  Note that "RefStrand" in build 37 annotations from Illumina is result of BLAST search of DESIGN strand.
  In previous strand annotations from Illumina for build 36 arrays, "BlastStrand" held +/- result of  BLAST search of the SOURCE sequence.
  
  If input "snp.dat" object does not contain a "RefStrand" column, then "alle.plus" mappings will not be returned. Thus the function is compatible with Illumina annotations lacking the "RefStrand" column (e.g., build 36 manifests).
}
\value{
  \code{make.allele.mappings} returns a data frame object with columns
  ("snp", "alle.AB", "alle.design", "alle.top", "alle.fwd",
  "alle.plus").

  \code{make.allele.annotation} returns a data frame with columns "snp",
  "alleleA.*", "alleleB.*" where "*"=\code{alleles}.
}
\author{Sarah Nelson}
\seealso{}
\examples{
\dontrun{
##  Load in Illumina annotation
column.select <- rep("NULL",times=21)
column.select[c(1:4,9:11,16:18,21)] <- NA  ## only read in select columns    
snp.dat <- read.csv(file="/projects/geneva/geneva_sata/SNP_annotation/Illumina/HumanOmni2.5_4v1/HumanOmni2.5-4v1-Multi-B.csv"),
                    skip=7,colClasses=column.select,nrows=2450000)   ## last rows are for control probes
map.final <- make.allele.mappings(snp.dat)
write.csv(map.final, file="/projects/geneva/geneva_sata/SNP_annotation/Illumina/HumanOmni2.5_4v1/testfn.csv",
          row.names=FALSE,quote=FALSE)

## Create a data frame with one row per SNP, TOP alleles, and indels coded as D/I
snp.annot <- make.allele.annotation(map.final, alleles="top")

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
#  alle.plus       text    nucleotide(s) corresponding to A or B allele
#  for PLUS(+) strand, relative to the forward direction in the human
#  reference genome sequence
}
}
\keyword{}
