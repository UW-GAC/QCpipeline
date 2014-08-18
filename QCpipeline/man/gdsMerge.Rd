\name{gdsMerge}
\alias{gdsMerge}
\alias{gdsMergeCheck}
\title{Merge multiple GenotypeData objects into a new file}
\description{Merges multiple GenotypeData objects on common SNPs and
  writes to a new file}
\usage{
gdsMerge(genoDataList, sampleList=NULL, snpList=NULL,
         sortByScanID=TRUE, newSnpID=TRUE, outPrefix="new")

gdsMergeCheck(genoDataList, outPrefix="new")
}
\arguments{
  \item{genoDataList}{named list of GenotypeData objects.
    list names will be used to identify datasets in output snp annotation}
  \item{sampleList}{optional list of scanIDs to include - list names same as genoDataList}
  \item{snpList}{optional list of snpIDs to include - list names same as genoDataList}
  \item{sortByScanID}{\code{TRUE} if scanIDs should be sorted in new GDS file (\code{FALSE} writes scans
    from each file in blocks)}
  \item{newSnpID}{\code{TRUE} if snpIDs in new file should be 1:n.  If
    \code{FALSE}, snpIDs will be taken from the first object in \code{genoDataList}.}
  \item{outPrefix}{output file prefix (will append ".gds,"
    "_snpAnnot.RData," and "_scanAnnot.RData.")}
}
\details{
\code{gdsMerge} merges SNPs that match on chromosome, position, and
alleles (but will flip any snps where A/B is reversed, e.g. C/T on one
array and T/C on the other).  SNPs in the output file are ordered by
chromosome and position.

It saves scan and snp annotation with
basic info.  Allele definitions are taken from the first
dataset.   The output snp annotation includes original snpIDs.

One can provide a list of scanIDs (\code{sampleList}) and snpIDs
(\code{snpList}) to include from each dataset, default
is to include all of them.
}
\author{Stephanie Gogarten}
%\seealso{}
%\examples{}
