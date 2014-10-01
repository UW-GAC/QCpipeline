\name{pvalFilterPlots}
\alias{qqPlotPng}
\alias{manhattanPlotPng}
\alias{calculateLambda}

\title{QQ and Manhattan plots with multiple panels and filtering}
\description{
  Plot p values in QQ and manhattan plots in multi-panel PNGs with
  different filter levels
}
\usage{
qqPlotPng(pval, stat, df, filters, outfile, ncol=2, addText="", main=NULL, ...)
manhattanPlotPng(pval, chromosome, filters, outfile, addText="", main=NULL, ...)
calculateLambda(stat, df)
}
\arguments{
  \item{pval}{p values to plot}
  \item{stat}{Test statistic from chisq distribution, e.g., \code{(Beta/SE)^2} for Wald tests}
  \item{df}{degrees of freedom for the test statistic (typically 1 for most GWAS models without interactions)}
  \item{chromosome}{chromosome titles for manhattan plot, same length as
    \code{pval}}
  \item{filters}{named list of filters for each planel of the plots.
    should have 4 elements for \code{qqPlotPng} and 3 elements for
    \code{manhattanPlotPng}.}
  \item{outfile}{name of output png file}
  \item{ncol}{number of columns in QQ plot}
  \item{addText}{string to add to the upper left corner of the plot}
  \item{main}{main title to add above all plots}
  \item{...}{other parameters passed to \code{\link{qqPlot}} or \code{\link{manhattanPlot}}}
}
\details{
  Generates plots with 4 panels (\code{qqPlotPng}) or 3 panels
  (\code{manhattanPlotPng})
  where each panel has the filtering provided in the corresponding
  element of \code{filters} and the names of \code{filters} are the plot
  titles.
  If \code{main} is specified, then the outer margins will be adjusted to fit the title.
  
  \code{calculateLambda} will calculate lambda given the test statistic (e.g., \code{Z^2} or likelihood ratio statistic) and the degrees of freedom of that statistic.
}
\author{Stephanie Gogarten}
\seealso{\code{\link{qqPlot}}, \code{\link{manhattanPlot}}}
%\examples{}
