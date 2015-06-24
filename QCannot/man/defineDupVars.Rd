\name{defineDupVars}
\alias{defineDupVars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Define duplicate variants in a SNP annotation based on position
and, optionally, alleles.
}
\description{
Identifies pairs or sets of variants with the same map position
(chromosome and base pair) and, optionally, also with the same alleles.
Returns integer IDs (e.g., 'dup.pos.id' and 'dup.var.id') for sets of
duplicate variants that can be used to select unique positions and/or
unique position-allele combinations from a set of variants.
}
\usage{
defineDupVars(snpAnnot, alleleAcol = "alleleA", alleleBcol = "alleleB",
supp_snpAnnot = NULL, supp_alleleAcol = "alleleA", supp_alleleBcol =
"alleleB", supp_snpIDs = NULL, match.alleles = FALSE,
start_dup.pos.id=1, start_dup.var.id=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpAnnot}{
SNP annotation data frame object (\code{SnpAnnotationDataFrame} as defined in
\pkg{GWASTools} package)
}
\item{alleleAcol}{
first allele column to match on, defaults to 'alleleA'  
}
\item{alleleBcol}{
second allele column to match on, defaults to 'alleleB'  
}
\item{supp_snpAnnot}{
optionally, a SNP annotation data frame object (\code{SnpAnnotationDataFrame} as defined in
\pkg{GWASTools} package) containing supplementary information to use for
aligning position and alleles
}
  \item{supp_alleleAcol}{
first allele column to match on in supplementary SNP annotation, defaults to 'alleleA'
}
\item{supp_alleleBcol}{
first allele column to match on in supplementary SNP annotation, defaults to 'alleleB'  
}
  \item{supp_snpIDs}{
optionally, vector of snpIDs for which to use supplementary SNP
annotation. Defaults to use supplementary annotation for all snpIDs found in common between primary and supplementary annotation.
}
\item{match.alleles}{
if \code{FALSE}, match only on chromosome and position. If \code{TRUE},
additionally match on alleles, allowing for strand discrepancies. Defaults to \code{FALSE}. 
}
\item{start_dup.pos.id}{
  Integer value at which to start dup.pos.id, defaults to 1.
}
\item{start_dup.var.id}{
  Integer value at which to start dup.var.id, defaults to 1.
}
}
\details{
'defineDupVars' takes a primary SNP annotation and determines
overlapping variants based on chromosome, position, and - optionally, if
\code{match.alleles}=\code{TRUE}, also on alleles. When matching on alleles, the order of
the two alleles does not need to be the same -- e.g., two variants at the
same position with alleleA=C, alleleB=T and alleleA=T, alleleB=C will
still be identified as overlapping. Furthemore, where alleles are
consistent yet appear to be annotated on opposite DNA strands,
overlap is still identified by checking the reverse complement of
specified alleles.

A supplementary SNP annotation may be provided when there is alternative
position and/or allele information that should be used to determine
overlap within the primary SNP annotation. For example, this may be
useful when identifying overlap at insertion/deletion (indel) variants
where records in the primary SNP annotation are not consistenly in VCF
convention (i.e., left-aligned indel alleles and
positions). Supplementary annotations may also be desired when
monomorphic variants in the primary SNP annotation have
information about only one of the alleles, or to integrate position
information from an updated build. 
}
\value{
Returns a data frame of the same dimension as \code{snpAnnot}, with the
following fields added:
\item{suppAnnotUsed}{Logical vector indicating where supplementary
  annotation was used to identify overlap in position and/or
  alleles. Returned only if a supplementary annotation is provided.}
\item{dup.pos.id}{Integer ID (defaults to 1:n) for duplicate variants when matching
  on position (chromosome and base pair position) only}
\item{dup.var.id}{Integer ID (defaults to 1:n) for duplicate variants when matching on
  position and alleles, allowing for discrepancies (i.e., 'strand
  flips') in annotated alleles. Returned only if match.alleles=\code{TRUE}}
\item{dupDiscStrand}{Logical vector indicating where plus strand allele
  annotation was discrepant between members of a \code{dup.var.id}
  pair/muliple. Returned only if match.alleles=\code{TRUE}}
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Sarah Nelson
}
\note{
Note that the \code{snpID} integer ID, which is a required field of all
\code{SnpAnnotationDataFrame} objects, is used to match records between the
primary and supplementary SNP annotation objects. Thus users shold make
sure a consistent \code{snpID} is used between the two annotation objects.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[GWASTools]{SnpAnnotationDataFrame}}
  \code{\link[Biostrings]{DNAStringSet}}
  \code{\link[Biostrings]{reverseComplement}}
  
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
# load SNP annotation object
snp.ann <- getobj("SOL_freeze2_chr-22_snpAnnot.RData")
dim(snp.ann)
##     snps variables 
##   381488        22 
# load supplementary annotation object
supp.ann <- getobj("SoL_freeze2_suppAnnot_monomorphs_indels_v0.RData")
out <- defineDupVars(snpAnnot=snp.ann, supp_snpAnnot=supp.ann, match.alleles=TRUE)
dim(out)
# [1] 381488      9
head(out[!is.na(out$dup.pos.id),])
##            snpID chromosome position alleleA alleleB suppAnnotUsed dup.pos.id
## 346     28273766         22 16423023       G       A         FALSE          1
## 347     28273767         22 16423023       G     GTA         FALSE          1
## 2429499 28274647         22 16890136       A       G          TRUE          2
## 1227    28274648         22 16890136       G       A         FALSE          2
## 1368    28274789         22 16904583       C       T         FALSE          3
## 359035  28274790         22 16904583       T       C         FALSE          3
##         dup.var.id dupDiscStrand
## 346             NA         FALSE
## 347             NA         FALSE
## 2429499          1         FALSE
## 1227             1         FALSE
## 1368             2         FALSE
## 359035           2         FALSE

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

