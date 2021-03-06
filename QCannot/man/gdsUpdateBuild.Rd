\name{gdsUpdateBuild}
\alias{gdsUpdateBuild}
\title{Update build of a GDS file}
\description{Updates the build of a GenotypeData object and writes a new GDS file}
\usage{
gdsUpdateBuild(genoData, map, outPrefix="new", remove.unmapped=TRUE,
               update.alleles=FALSE, block.size=100)
}
\arguments{
\item{genoData}{A \code{\link{GenotypeData}} object}
\item{map}{a data frame with names ("old.rsID", "old.chromosome", "old.position",
                "new.rsID", "new.chromosome", "new.position")}
\item{outPrefix}{output file prefix (will append ".gds," and "_snpAnnot.RData")}
\item{remove.unmapped}{Logical for whether to remove unmapped SNPs.}
\item{update.alleles}{Logical for whether to update the alleles at the same time.  This option requires additional columns in \code{map}: ("old.alleleA", "old.alleleB", "new.alleleA", "new.alleleB")}
\item{block.size}{Number of samples to process at once}
}
\details{
Takes a GenotypeData object and converts from one build to another based on a user-specified map.  Writes a new GDS file and SNP annotation with the updated build.  Alleles can optionally be updated as well.

\code{genoData} must contain a SNP annotation column "rsID" with no duplicate values.

If \code{remove.unmapped=FALSE}, all SNPs in \code{genoData} must have entries in \code{map}.
}
\author{Stephanie Gogarten}
\seealso{\code{\link{convertBuild}}, \code{\link[QCpipeline]{gdsMerge}}}
\examples{
\dontrun{
See /projects/geneva/geneva_sata/QC_reports/SOPs/datasetMerging/exampleDatasets/dataMerging_example.pdf
}
}
