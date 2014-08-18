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
\details{\code{plateLayoutPlot} generates a pdf showing samples with identity issues. Currently only works with arrays that have 96 wells (A-H on one axis and 1-12 on the other). The columns of the data frame must include:
\itemize{
\item{\code{scanID}:}{ sample scanID }
\item{\code{plate}:}{ sample plate}
\item{\code{well}:}{ sample well}
\item{\code{annot.sex}:}{ annotated sex, values \code{"M"}, \code{"F"}, \code{"U"}}
\item{\code{sex}:}{ genetic sex, values \code{"M"}, \code{"F"}, \code{"U"}}
\item{\code{obsdup.id}:}{ integer identifier for observed duplicate set (i.e., all members of a set are dups of one another)}
\item{\code{unexp.dup}:}{ logical indicator for dup sets that were not expected to be duplicates}
\item{\code{unobsdup.id}:}{ integer identifier for expected dupliates in which the members were not observed to be duplicates}
}
}
%\value{}
\author{Cathy Laurie, Adrienne Stilp}
%\seealso{}
%\examples{}
