\name{plot2DwithHist}
\alias{plot2DwithHist}
\title{Scatterplot with density}
\description{
  \code{plot2DwithHist} produces a scatterplot of y vs x, along with histograms of
  the marginal distributions of x and y.}
\usage{
plot2DwithHist(x, y, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
  sublab=NULL, col2D="#0000ff22", mn=NULL, sd=NULL)
}
\arguments{
  \item{x}{vector of x coordinates}
  \item{y}{vector of y coordinates}
  \item{xlab}{x-axis label (defaults to variable name)}
  \item{ylab}{y-axis label (defaults to variable name)}
  \item{xlim}{x-axis limits (defaults to [min,max] of X, plus a bit of space}
  \item{ylim}{y-axis limits (defaults to [min,max] of Y, plus a bit of space}
  \item{sublab}{sub-label (instead of main, since there's no room)}
  \item{col2D}{color for plotting points}
  \item{mn}{2-element vector with mean of x and y}
  \item{sd}{2-element vector with sd of x and y}
}
\details{}
\author{Leila Zelnick}
\seealso{}
\examples{
library(MSBVAR)
# generate some multivariate normal example data
n <- 5000
mu <- c(0, 2)
vmat <- matrix(c(1, 0.7, 0.7, 1), nrow=2)

dat <- rmultnorm(n, mu, vmat) # generates n multivariate normal obs.
x <- dat[,1]
y <- dat[,2]

plot2DwithHist(x, y, xlab="This is the X variable", ylab="This is the Y variable.",
  sub="Example Plot!")
# defining axis limits
plot2DwithHist(x, y, xlab="This is the X variable", ylab="This is the Y variable.",
  sub="Example Plot!", xlim=c(0,4), ylim=c(-2,2))
}
\keyword{}
