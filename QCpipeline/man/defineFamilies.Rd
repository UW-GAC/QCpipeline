\name{defineFamilies}
\alias{defineFamilies}
\title{Define families}
\description{Define families based on kinship coefficients}
\usage{
defineFamilies(ibd, KC.threshold=0.1, start.fam.id=1)
}
\arguments{
  \item{ibd}{Data frame with columns "sample1", "sample2", "KC"}
  \item{KC.threshold}{Minimum KC value for samples to belong to the same
  family}
  \item{start.fam.id}{Starting value for family ID.}
}
\details{}
\value{
  Returns \code{ibd} with an additional "family" column.  Families will
  be numbered sequentially starting with \code{start.fam.id}.
}
\author{Cathy Laurie}
\seealso{}
\examples{
}
\keyword{}
