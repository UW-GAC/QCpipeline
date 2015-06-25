\name{multiplot}
\alias{multiplot}
\title{Multiple plot function}
\description{Combine multiple ggplot objects into one plot}
\usage{
multiplot(..., plotlist=NULL, cols=1, layout=NULL, letters=FALSE, byrow=FALSE)
}
\arguments{
  \item{...}{ggplot objects}
  \item{plotlist}{list of ggplot objects}
  \item{cols}{Number of columns in layout}
  \item{layout}{A matrix specifying the layout. If present, 'cols' is ignored.}
  \item{letters}{logical for whether to add letters A, B, etc. to the top left of each plot (for multi-panel figures in publications)}
  \item{byrow}{logical for whether the plots should be printed by row or by column}
}
\details{
Taken from 
http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects).

If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
then plot 1 will go in the upper left, 2 will go in the upper right, and
3 will go all the way across the bottom.
}
\author{Stephanie Gogarten}
%\seealso{}
\examples{
library(ggplot2)
dat <- data.frame(a=1:10, b=11:20)
g1 <- ggplot(dat, aes(x=a, y=b)) + geom_point()
g2 <- ggplot(dat, aes(x=a, y=b)) + geom_boxplot()
g3 <- ggplot(dat, aes(x=a)) + geom_bar()
multiplot(g1, g2, g3, layout=matrix(c(1,2,3,3), nrow=2, byrow=TRUE), letters=TRUE)
}
