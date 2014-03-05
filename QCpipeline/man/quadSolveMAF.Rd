\name{quadSolveMAF}
\alias{quadSolveMAF}
\title{Effective sample size MAF filter}
\description{
  Function to find effective sample size MAF filter using quadratic
  equation
}
\usage{
quadSolveMAF(X, N)
}
\arguments{
  \item{X}{filter threshold}
  \item{N}{sample size}
}
\details{
  Returns the effective MAF filter threshold that satisfies
  2*MAF*(1-MAF)*N >X.
}
\author{Adrienne Stilp}
\seealso{}
\examples{
X <- 30
N <- c(1000, 3000, 10000)
quadSolveMAF(X, N)
}
