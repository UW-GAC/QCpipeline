\name{checkABalleles}
\alias{checkABalleles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to check consistency between probe sequences and A/B alleles in an Illumina manifest
}
\description{
Checks consistency between terminal nuceleotide in A/B probe sequences and
the A/B design alleles as indicated in the 'SNP' column, for Illumina
SNP array manifests. Only applies to strand ambiguous SNPs, which have
two probe types (i.e. Infinium I assays)
}
\usage{
checkABalleles(snp, identifier = "Name", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp}{
Data frame of Illumina SNP annotation, including the following columns:
an identifier, 'SNP,' 'AlleleA_ProbeSeq,' and 'AlleleB_ProbeSeq'
}
\item{identifier}{
name of field to use as identifier (defaults to 'Name')
}
\item{verbose}{
logical for whether to report on function progress  
}
}
\details{
Written after discovering inconsistencies in HumanOmni2.5-4v1_D.csv array manifest used in HRS1-2 project, which caused many errors (1) when merging with HRS3 and (2) imputing HRS1-2 to 1000 Genomes  
}
\value{
returns a subset of input data frame with discrepancies, as a data frame
}
\references{

Nelson, S.C., Laurie, C.C., Doheny, K.F. & Mirel, D.B. Is 'forward' the same as 'plus'?...and other adventures in SNP allele nomenclature. Trends in Genetics 28, 361-363 (2012).  
  
}
\author{
 S. Nelson
}
\note{
Note the HumanOmni2.5-4v1_D.csv manifest is the only one detected thus
far to have this annotation issue.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

illumina <-  read.csv(file="/projects/geneva/gcc-fs2/SNP_annotation/Illumina/HumanOmni2.5_4v1/HumanOmni2.5-4v1_D.csv", skip=7, nrows=10000)
dim(illumina); names(illumina)
out <- checkABalleles(illumina)
dim(out); head(out)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
