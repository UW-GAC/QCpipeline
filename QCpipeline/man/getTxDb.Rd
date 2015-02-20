\name{getTxDb}
\alias{getTxDb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get a TxDb object
}
\description{
Create a transcript database (TxDb) object from UCSC Genome
Bioinformatics resources.
}
\usage{
getTxDb(build, database_table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{build}{
genome build, e.g 'hg18' (build 36) or 'hg19' (build 37)
}
\item{database_table}{
  name of database table to get: e.g.,  'refGene' for  RefSeq refGene
  table, under track Ref Seq Genes; 'knownGene' for knownGene table,
  under track UCSC genes...or other tables listed by \code{\link{supportedUCSCtables}}
}
}
\details{
 Specifying \code{database_table='knownGene'} will make use of the
 pre-existing annotation packages
 \code{\link{TxDb.Hsapiens.UCSC.hg19.knownGene}} (for build hg19) and
 \code{\link{TxDb.Hsapiens.UCSC.hg18.knownGene}} (for build hg18). 
 Other \code{database_table} arguments will use the
 \code{\link{makeTranscriptDbFromUCSC}} function to create a TxDb
 object, pulling the annotation from UCSC Genome Browser.
}
\value{
Returns a \code{\link{TxDb}} object, which can either be (1) saved as an
SQLite file for later use, with
\code{\link{saveDb}} and/or (2) passed as an argument to the companion \code{\link{defineExomeVars}} function
}
\references{
  \url{http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf}
  \url{http://genome.ucsc.edu/cgi-bin/hgTables}
}
\author{
Sarah Nelson
}
% \note{
% %%  ~~further notes~~
% }

\seealso{
  \code{\link{TxDb}}, \code{\link{saveDb}}, \code{\link{loadDb}},
  \code{\link{TxDb.Hsapiens.UCSC.hg19.knownGene}},
  \code{\link{TxDb.Hsapiens.UCSC.hg18.knownGene}}, \code{\link{makeTranscriptDbFromUCSC}}

}
\examples{
\dontrun{
# UCSC genes, build 37
txdb_knownGene <- getTxDb("hg19", "knownGene")

# RefSeq genes, build 36
txdb_refGene <- getTxDb("hg18", "refGene")
}
}

\keyword{manip}
