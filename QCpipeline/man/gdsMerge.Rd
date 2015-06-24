\name{gdsMerge}
\alias{gdsMerge}
\alias{gdsMergeCheck}
\title{Merge multiple GenotypeData objects into a new file}
\description{Merges multiple GenotypeData objects on common SNPs and
  writes to a new file}
\usage{
gdsMerge(genoDataList, outPrefix="new", sampleList=NULL, snpList=NULL,
         match.snps.on="position", snpNameList=NULL,
         sortByScanID=TRUE, newSnpID=TRUE, dosage=FALSE, verbose=TRUE)

gdsMergeCheck(genoDataList, outPrefix="new", verbose=TRUE)
}
\arguments{
  \item{genoDataList}{named list of GenotypeData objects.
    list names will be used to identify datasets in output snp annotation}
  \item{outPrefix}{output file prefix (will append ".gds,"
    "_snpAnnot.RData," and "_scanAnnot.RData.")}
  \item{sampleList}{optional list of scanIDs to include - list names same as genoDataList}
  \item{snpList}{optional list of snpIDs to include - list names same as genoDataList}
  \item{match.snps.on}{One or more of ("position", "name") indicating how to match SNPs.  "position" will match SNPs on chromosome and position, and "name" will match on the columns give in \code{snpNameList}.  The same alleles are always required for matching (but A/B designations need not be the same).}
  \item{snpNameList}{list of SNP annotation columns to match on if \code{match.snps.on="name"} - list names same as genoDataList}
  \item{sortByScanID}{\code{TRUE} if scanIDs should be sorted in new GDS file (\code{FALSE} writes scans
    from each file in blocks)}
  \item{newSnpID}{\code{TRUE} if snpIDs in new file should be 1:n.  If
    \code{FALSE}, snpIDs will be taken from the first object in \code{genoDataList}.}
  \item{dosage}{\code{TRUE} if any of the items in \code{genoDataList} are dosage files (genotypes stored in 'float32' instead of 'bit2' format). If \code{dosage=TRUE}, the output file will also be a dosage file.}
  \item{verbose}{logical for whether to print messages}
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
\seealso{\code{\link{gdsCombine}}, \code{\link[QCannot]{gdsUpdateBuild}}}
\examples{
\dontrun{
See /projects/geneva/geneva_sata/QC_reports/SOPs/datasetMerging/exampleDatasets/dataMerging_example.pdf
}
}
