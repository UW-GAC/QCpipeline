\name{getSnpNames}
\alias{getSnpNames}

\title{
getSnpNames
}
\description{
Get rsIDs based on SNP chromosome and position.
}
\usage{
getSnpNames(snpAnnot, gds, chromosome, verbose = TRUE, extraSnptabVars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpAnnot}{\code{\link{SnpAnnotationDataFrame}} object to find rsIDs for.}
  \item{gds}{Opened gds file of SNP table.}
  \item{chromosome}{Chromosome to find names for (see details).}
  \item{verbose}{report progress information?}
  \item{extraSnptabVars}{Any additional variables to read in from the SNP table gds file (see details).}
}
\details{
\code{getSnpNames} matches rsIDs from the SNP table to a set of SNPs using chromosome and position. The input gds file must have nodes \code{name}, \code{chrom}, \code{chromStart}, \code{chromEnd}, \code{observed}, \code{strand}, \code{class}, and \code{exceptions}. By default, these columns are returned. If additional columns are desired, specify them with \code{extraSnptabVars}

At this point, only one chromosome can be run at a time. \code{snpAnnot} must have only one unique chromosome. Future development may allow multiple chromosomes to be run together.

\code{getSnpNames} only matches SNPs from the SNP annotation whose alleleA and alleleB are one of \code{A}, \code{C}, \code{G}, \code{T}, or \code{0}, i.e., indels or CNVs will not have a matched name even if they also appear in the SNP table. These SNPs are indicated by the column \code{snpAnnot.valid}. The function only matches SNPs from the SNP table that have \code{class='single'}, SNPs that are not multiply-aligned (indicated by the \code{exception} column containing \code{'MultipleAlignments'}), and SNPs that have a range of only one basepair. These SNP records are indicated with the column \code{snptab.valid} in the output. If the SNP position is not in the SNP table, then \code{snptab.valid} will be \code{NA}.

Only the valid SNPs from both the gds SNP table and from the SNP annotation are merged to find rsIDs. The output of this merge is the re-merged back with the original SNP annotation by snpID. In some cases, there may be duplicate, valid rsIDs for a given snpID. In this case, the output will have a duplicate record for the given snpID so that the user can decide which rsID to keep.

}

\value{
\item{merged}{A data frame containing the SNP annotation merged with the SNP table. The output rsID is in the column \code{name}. In some cases, there may be duplicate records for a given snpID (see details). The merged data frame also contains the \code{snptab.valid} and \code{snpAnnot.valid} columns.}
\item{snptab}{Requested columns from the SNP table for all SNPs on the requested chromosome, stored as a data frame. The function adds a column \code{snptab.valid}. SNP records with \code{snptab.valid=TRUE} indicate that they were considered acceptable SNPs for merging (see details).}
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Adrienne Stilp
}


\seealso{
\code{\link{readSnpTable}}
}
\examples{
\dontrun{
library(gdsfmt)
gds <- openfn.gds("snp_table.gds"")
snpAnnot <- getobj(snpfile) # snp annotation for chromosome 1
out <- getSnpNames(snpAnnot, gds, 1)

out <- getSnpNames(snpAnnot, gds, 1, extraSnptabVars=c("refUCSC"))
}
}
