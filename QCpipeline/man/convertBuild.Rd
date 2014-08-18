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
convertBuild(snpAnnot, chain.file, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpAnnot}{
Open SNP annotation object (object of class 'SnpAnnotationDataFrame') with "chromosome" and "position" columns in original build
}
  \item{chain.file}{
Path to UCSC chain file, from original build to target build
}
  \item{verbose}{
logical flag to report out counts of successful and unsuccessful conversion
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
 Returns input SNP annotation as a data frame with added columns: chromosome and position in the target build
}
\references{
http://www.bioconductor.org/packages/release/bioc/manuals/rtracklayer/man/rtracklayer.pdf
http://genome.ucsc.edu/cgi-bin/hgLiftOver
}
\author{
Sarah Nelson
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
