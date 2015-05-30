\name{simulateGenotypeData}
\alias{simulateImputedGenotypeData}
\alias{simulateAnnotation}
\alias{simulateAssocResults}

\title{Simulate Genotype Data}

\description{
Simulate genotype data for testing.
}

\usage{
simulateImputedGenotypeData(path, prefix, nsamp=50, nsnp=100, nchromosome=23, nsegments=3, missing=NULL, filterColumn="oevar")
simulateAnnotation(path, prefix, nsnp=100, nchromosome=23)
simulateAssocResults(path, prefix, nsnp=100, nchromosome=23)
}

\arguments{
  \item{path}{directory for output genotype chromosome files for imputed data}
  \item{prefix}{string prefix for genotype chromosome files for imputed data}
  \item{nsamp}{number of samples to simulate}
  \item{nsnp}{number of SNPs per chromosome to simulate}
  \item{nchromosome}{number of chromosomes to simulate}
  \item{nsegments}{number of segments per chromosome}
  \item{missing}{fraction of SNPs to set to missing. If null, no SNPs are set to missing.}
  \item{filterColumn}{name of column in SNP annotation indicating the SNP quality filter}
}

\details{
\code{simulateImputedGenotypeData} makes a directory of simulated gentoype files, one per chromosome, with associated SNP annotations and snp-segment mapping. Files are located in the directory specified by \code{path}, with names like "\code{prefix}_chr-1.gds" and "\code{prefix}_chr-1_snpAnnot.RData".

\code{simulateAnnotation} makes a directory of simulated annotation files, one per chromosome.

\code{simulateAssocResults} makes a directory of simulated association test files, one per chromosome.
}
%\value{
%}
\author{
Adrienne Stilp, Stephanie Gogarten
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
%\examples{
%}

