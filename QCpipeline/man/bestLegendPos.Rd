\name{bestLegendPos}
\alias{bestLegendPos}
\title{Best legend position}
\description{Find the section of a plot with the fewest number of points
for placing a legend}
\usage{
bestLegendPos(x, y)
}
\arguments{
  \item{x}{Points on the x axis}
  \item{y}{Points on the y axis}
}
%\details{}
\value{
  Returns one of the legend placement strings from the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}.
}
\author{Stephanie Gogarten}
%\seealso{}
\examples{
x <- sample(1:100, 50, replace=TRUE)
y <- sample(1:100, 50, replace=TRUE)
plot(x,y)
legend(bestLegendPos(x,y), "points", pch=1)
}
