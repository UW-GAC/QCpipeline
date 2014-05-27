\name{readSnpTable}
\alias{readSnpTable}

\title{
readSnpTable
}
\description{
Read information from the SNP table, stored in gds format.
}
\usage{
readSnpTable(gds, chromosome, variables, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gds}{An opened gds file that stores SNP table information.}
  \item{chromosome}{The chromosome to read. See details.}
  \item{variables}{Variables to read in from the SNP table.}
  \item{verbose}{Logical indicatoring whether to report progress information.}
}
\details{
\code{readSnpTable} will read SNP table information from a gds file for one chromosome at a time. The chromosome specified will typically be 1-22, X, Y, or M; other chromosomes are present in the file but generally will not be used. The \code{variable} argument allows users to specify which columns of the table they would like information from. Each requested variable must be stored as a node in the gds file. Unfortunatley, the available SNP table columns change based on the genome build. Use \code{ls.gdsn(gds)} to see which nodes are present in the gds file.
}
\value{
A data frame with one row per SNP and one column per requested variable.
}
\author{
Adrienne Stilp
}

\references{
http://genome.ucsc.edu/cgi-bin/hgTables
}

\seealso{
\code{getSnpNames}
}
\examples{
\dontrun{
gds <- openfn.gds(filename.gds)
snptab <- readSnpTab(gds, 1, c("name", "chromStart", "chromEnd", "exceptions"))
closefn.gds(gds)
}
}
