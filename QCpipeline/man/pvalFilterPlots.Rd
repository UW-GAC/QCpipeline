\name{pvalFilterPlots}
\alias{qqPlotPng}
\alias{manhattanPlotPng}
\title{QQ and Manhattan plots with multiple panels and filtering}
\description{
  Plot p values in QQ and manhattan plots in multi-panel PNGs with
  different filter levels
}
\usage{
qqPlotPng(pval, filters, outfile)
manhattanPlotPng(pval, chromosome, filters, outfile)
}
\arguments{
  \item{pval}{p values to plot}
  \item{chromosome}{chromosome titles for manhattan plot, same length as
    \code{pval}}
  \item{filters}{named list of filters for each planel of the plots.
    should have 4 elements for \code{qqPlotPng} and 3 elements for
    \code{manhattanPlotPng}.}
  \item{outfile}{name of output png file}
}
\details{
  Generates plots with 4 panels (\code{qqPlotPng}) or 3 panels
  (\code{manhattanPlotPng})
  where each panel has the filtering provided in the corresponding
  element of \code{filters} and the names of \code{filters} are the plot
  titles.
}
\author{Stephanie Gogarten}
\seealso{\code{\link{qqPlot}}, \code{\link{manhattanPlot}}}
\examples{}
