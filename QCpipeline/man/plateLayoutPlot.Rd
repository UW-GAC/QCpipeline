\name{plateLayoutPlot}
\alias{plateLayoutPlot}
\title{Plot problematic samples by plate and well}
\description{Plots problematic samples by plate and well}
\usage{
plateLayoutPlot(samples, out.file, nprob)
}
\arguments{
  \item{samples}{data frame with columns "scanID", "plate", "well", "annot.sex", "sex", "obsdup.id", "unexp.dup", "unobsdup.id"}
  \item{out.file}{output plot PDF file}
  \item{nprob}{minimum number of problems required to plot a plate}
  }
\details{\code{plateLayoutPlot} generates a pdf showing samples with identity issues. Currently only works with arrays that have 96 wells (A-H on one axis and 1-12 on the other).}
\value{
  
}
\author{Cathy Laurie}
\seealso{}
\examples{

}
\keyword{}
