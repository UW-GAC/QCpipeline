\name{defineExomeVars}
\alias{defineExomeVars}

\title{
Annotate variants with respect to genomic location and gene function.
}
\description{
  This function annotates variants with respect to genomic location and gene function:
  intergenic, intronic, exonic, promoter, splice site, or untranslated
  region (5' and 3' UTR). 
}
\usage{
defineExomeVars(snpAnnot, txdb, collapsed = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpAnnot}{\code{\link{SnpAnnotationDataFrame}} object to annotate}
  \item{txdb}{A \code{\link{TxDb}} object to use for the annotation}
  \item{collapsed}{logical, whether to return one row per input variant versus one
    row for each location type}
}
\details{
  This function uses the \code{\link{locateVariants}} function to
  annotate variants. The \code{txdb} argument refers to a
  \code{\link{TxDb}} object, which the user can either:

      - create with the function\code{\link{getTxDb}}, a companion
      function to \code{defineExomeVars}.

  or
  
      - load from a previously created and saved \code{\link{TxDb}}
      object, using the \code{\link{saveDb}} and \code{\link{loadDb}}
      functions.
  
  Users should check that (1) variants in the \code{snpAnnot} are in a
  consistent genome build and that (2) the \code{txdb} object is also in
  the same build. Genome build is an argument to the
  \code{\link{getTxDb}} function.

  When run with \code{collapsed=TRUE}, one row per variant will be
  returned, with different location types and geneIDs collapsed into one string with
  a ';' separator. When run with \code{collapsed=FALSE}, one row will be
  returned for each query result.

  All the variants included in the \code{snpAnnot} input will be
  returned in the output. The returned values from the
  \code{\link{locateVariants}} function will be \code{NA} where
  the variant was either was not included in the query or when no
  information was found for that variant. Reasons for not including a
  variant in the query are (1) unmapped position or (2) chromosome value other
  than an autosome, X, Y, or mitochondrial.
}
\value{
  A data frame is returned with the following columns:
  \item{snpID}{integer snpID, from the input \code{snpAnnot}}
  \item{chromosome}{integer chromosome, from the input \code{snpAnnot}}
  \item{position}{base pair position, from the input \code{snpAnnot}}
  \item{excl.chr}{logical indicator for whether variant was excluded from query due
    to chromosome value other than an autosome, X, Y, or mitochondrial}
  \item{excl.pos}{logical inciator for whether variant was excluded from query due
   to position=0 (unmapped)}
  \item{loctype}{Values returned from \code{\link{locateVariants}}. Possible values
  are 'coding', 'intron', 'threeUTR', 'fiveUTR', 'intergenic',
  'spliceSite', or 'promoter'. If \code{collapsed=TRUE}, multiple
  records are concatenated by ';', with \code{NA} values removed.} 
  \item{geneID}{Values returned from \code{\link{locateVariants}}. For
  genic SNPs, one or more Entrez GeneIDs; \code{NA} otherwise. If
  \code{collapsed=TRUE}, multiple records are concatenated by ';', with
  \code{NAs} removed.} 
  \item{exomic}{logical, \code{TRUE} where variant is located in exomic region (at least one \code{loctype} value of 'coding'); \code{FALSE} if not located in exomic region; \code{NA} if not included in query (because \code{excl.chr=TRUE}, \code{excl.pos=TRUE}, or no match found by \code{\link{locateVariants}} query).} 
}
\references{
  \url{http://www.bioconductor.org/help/workflows/annotation/AnnotatingRanges/}
  \url{http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf}
  \url{http://genome.ucsc.edu/cgi-bin/hgTables}
  \url{http://www.ncbi.nlm.nih.gov/gene}
}
\author{
Sarah Nelson
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{locateVariants}}, \code{\link{TxDb}}, \code{\link{saveDb}},
      \code{\link{loadDb}}
}
\examples{
\dontrun{

library(GenomicFeatures)
library(GWASTools)
library(reshape2)

# load snps to annotate
snpfile <- system.file("data", "illuminaSnpADF.RData",
  package="GWASdata")
snpAnnot <- getobj(snpfile)

# get UCSC knownGene table, build hg19, as TxDb object
txdb <- getTxDb("hg19", "knownGene")

# example run with default, collapse=TRUE
out <- defineExomeVars(snpAnnot, txdb)
dim(out) 
# [1] 3300    8
names(out)
# [1] "snpID"      "chromosome" "position"   "excl.chr"   "excl.pos"  
# [6] "loctype"    "geneID"     "exomic"

# preview variants with multiple loctypes
ids <- head(out$snpID[grep(";",out$loctype)],3)
out[is.element(out$snpID, ids),]
#     snpID chromosome position excl.chr excl.pos              loctype    geneID
# 18  999757         21 14778884    FALSE    FALSE promoter; intergenic 100423018
# 43 1000085         21 15459599    FALSE    FALSE   intron; intergenic        NA
# 44 1000089         21 15479497    FALSE    FALSE   intron; intergenic        NA
#   exomic
# 18  FALSE
# 43  FALSE
# 44  FALSE

# example run with collapse=FALSE
out <- defineExomeVars(snpAnnot, txdb, collapsed=FALSE)
out[is.element(out$snpID, ids),]
#      snpID chromosome position excl.chr excl.pos    loctype    geneID exomic
# 18  999757         21 14778884    FALSE    FALSE   promoter 100423018  FALSE
# 19  999757         21 14778884    FALSE    FALSE intergenic      <NA>  FALSE
# 44 1000085         21 15459599    FALSE    FALSE     intron      <NA>  FALSE
# 45 1000085         21 15459599    FALSE    FALSE intergenic      <NA>  FALSE
# 46 1000089         21 15479497    FALSE    FALSE     intron      <NA>  FALSE
# 47 1000089         21 15479497    FALSE    FALSE intergenic      <NA>  FALSE

# example of creating, saving, and loading TxDb object
txdb_refgene <- getTxDb("hg19", "refGene")

# print information about the TxDb object
txdb_refgene
# TxDb object:
# | Db type: TxDb
# | Supporting package: GenomicFeatures
# | Data source: UCSC 
# | Genome: hg19
# | Organism: Homo sapiens
# | UCSC Table: refGene
# | Resource URL: http://genome.ucsc.edu/
# | Type of Gene ID: Entrez Gene ID
# | Full dataset: yes
# | miRBase build ID: NA
# | transcript_nrow: 53431
# | exon_nrow: 253728
# | cds_nrow: 212598
# | Db created by: GenomicFeatures package from Bioconductor
# | Creation time: 2015-02-19 15:06:34 -0800 (Thu, 19 Feb 2015)
# | GenomicFeatures version at creation time: 1.18.3
# | RSQLite version at creation time: 1.0.0
# | DBSCHEMAVERSION: 1.0

# save TxDb as a SQLite datbase file
saveDb(txdb,file="TxDb.hg19.refGene.sqllite")

# load SQLite database file
txdb_refgene <- loadDb("TxDb.hg19.refGene.sqllite")

# use as txdb function argument
out <- defineExomeVars(snpAnnot, txdb_refgene)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{attribute}
