\name{gdsCombine}
\alias{gdsCombine}
\alias{checkGdsCombine}

\title{gdsCombine
}

\description{
Combine multiple GDS files into one file.
}
\usage{
gdsCombine(gdsGenoList,
                       filename,
                       blockSize=5000,
                       snpExcludeList=NULL,
                       verbose=TRUE,
                       dryRun=FALSE)
checkGdsCombine(genoData, gdsGenoList, snp, snpExcludeList=NULL, verbose=TRUE, bySnp = TRUE, blockSize = 5000)
}
\arguments{
  \item{gdsGenoList}{Named list of \code{\link{GdsGenotypeReader}} objects.}
  \item{filename}{Filename of combined GDS file to create.}
  \item{blockSize}{Number of SNPs to operate on at a time (ignored in \code{checkGdsCombine} if \code{bySnp=FALSE})}
  \item{snpExcludeList}{list with snpIDs to exclude for each dataset, if any}
  \item{verbose}{logical: report progress info?}
  \item{dryRun}{logical indicator for whether to write output file (dryRun=TRUE) or just test (dryRun=FALSE)}
  \item{genoData}{\code{\link{GdsGenotypeReader}} or \code{\link{GenotypeData}} object of combined GDS file.}
  \item{snp}{Output from \code{gdsCombine}: SNP data frame mapping SNPs to original dataset.}
  \item{bySnp}{logical: check by snp?}
}
\details{
\code{gdsCombine} combined multiple GDS files into a single file. The component files in \code{genoDataList} are required to have the same scanIDs. All SNPs from the component files are included in the component file.
}
\note{
This function will need to be updated for getGenotypeSelection once it's implemented (instead of reading directly from the gds genotype node).
}
\value{
If given, snpExcludeList should have a subset of the names of \code{gdsGenoList}.

\code{snp} is returned by \code{gdsCombine} and is a data frame of the SNP annotation with columns:
\itemize{
\item{\code{snpID}:} new snpID
\item{\code{dataset}:} dataset of original SNP, taken from \code{names{genoDataList}}
\item{\code{snpID.original}:} snpID in the original dataset
\item{\code{chromosome}:} chromosome integer
\item{\code{position}:} position integer
\item{\code{alleleA}:} alleleA, taken from the original dataset
\item{\code{alleleB}:} alleleB, taken from the original dataset
\item{\code{id}:} unique identifier created from the SNP data frame using \code{paste(dataset, snpID.original)}
}
}
\author{Adrienne Stilp}

\seealso{
\code{\link{GdsGenotypeReader}}, \code{\link{GenotypeData}}
}
