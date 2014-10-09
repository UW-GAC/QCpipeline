\name{convertBuild}
\alias{convertBuild}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert genome build 
}
\description{
Converts genome build in a SNP annotation, using R implementation of UCSC liftOver tool  
}
\usage{
convertBuild(snpAnnot, chain.file, snp.include = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpAnnot}{
A SNP annotation object (object of class SnpAnnotationDataFrame) with
\dQuote{chromosome} (integer chromosome code) and \dQuote{position} columns in original build,
and column "rsID" 
}
  \item{chain.file}{
Path to UCSC chain file, from original build to target build
}
  \item{snp.include}{
Optionally, a integer vector containing the snpIDs of SNPs to be
included; defaults to including all SNPs in snpAnnot
}
  \item{verbose}{
logical flag to report out counts of successful and unsuccessful conversion
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 Returns input SNP annotation as a data frame with added columns:
 chromosome and position in the target build
 (\dQuote{chromosome.converted},\dQuote{position.converted})
}
\note{ 
 Note that pseudoautosomal SNPs are input as either chromosome X or
 Y, depending on initial chromosome codes (i.e. XY SNP with base pair position
 on the X or position on the Y). In the output, converted chromosome is
 given as X or Y, meaning users may need to re-evaluate whether X
 and Y SNPs belong to the pseudoautosomal regions of those chromosomes,
 based on the updated base pair positions, and re-assign chromosome
 codes accordingly.
}
\references{
http://www.bioconductor.org/packages/release/bioc/manuals/rtracklayer/man/rtracklayer.pdf
http://genome.ucsc.edu/cgi-bin/hgLiftOver
}
\author{
Nelson, S.C.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
 library(OLGApipline)
 data.dir <- "/projects/geneva/gcc-fs2/OLGA/genotype/freeze1/gds/freeze1"
 olgaData <- OlgaGenotypeData(data.dir)
 snpAnnot <- getSnpAnnotation(olgaData, 21)
 chain.file <- "/projects/geneva/gcc-fs2/SNP_annotation/UCSC_downloads/hg19/liftOver/hg19ToHg18.over.chain"
 converted <- convertBuild(snpAnnot, chain.file)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{manip}
