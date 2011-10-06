# This script is from Leila Zelnick, Feb 2011
# added mean and SD options

# plot.hist produces a scatterplot of y vs x, along with histograms of
# the marginal distributions of x and y
plot2DwithHist <- function(x, y,
  xlab=NULL, # x-axis label (defaults to variable name)
  ylab=NULL, # y-axis label (defaults to variable name)
  xlim=NULL, # x-axis limits (defaults to [min,max] of X, plus a bit of space
  ylim=NULL, # y-axis limits (defaults to [min,max] of Y, plus a bit of space
  sublab=NULL, # sub-label (instead of main, since there's no room)
  col2D="#0000ff22",
  mn=NULL, # 2-element vector with mean of x and y
  sd=NULL  # 2-element vector with sd of x and y
){ # number of bins for histograms: defaults to hist default

  if(is.null(xlim)) {
    x.range <- max(x, na.rm=T) - min(x, na.rm=T)
    xlim <- c(min(x, na.rm=T)-0.1*x.range, max(x, na.rm=T)+0.1*x.range)}
  if(is.null(ylim)) {
    y.range <- max(y, na.rm=T) - min(y, na.rm=T)
    ylim <- c(min(y, na.rm=T)-0.1*y.range, max(y, na.rm=T)+0.1*y.range)}
  
  layout(matrix(c(2,4,1,3), 2, 2, byrow=T), heights=c(1,3), widths=c(3,1)) 
  # scatterplot
  par(mar=c(5,4,1.5,1.5)+0.1, las=1)
  plot(x, y, xlab=xlab, ylab=ylab, sub=sublab,
    cex.sub=1.2, type="n", xlim=xlim, ylim=ylim)
  points(x, y, col=col2D)
  # mean and SD
  if(!is.null(mn)) {
    abline(v=mn[1], h=mn[2])
  }
  if(!is.null(sd)) {
    for(i in 1:10) {
      abline(v=mn[1]+i*sd[1], h=mn[2]+i*sd[2], lty=2, col="gray")
      abline(v=mn[1]-i*sd[1], h=mn[2]-i*sd[2], lty=2, col="gray")
    }
  }
  # x density
  par(mar=c(1.5,4,2,1.5)+0.1)
  max.density.x <- max(density(x)$y, na.rm=T)
  plot(density(x), xlim=xlim, ylim=c(0, 1.1*max.density.x), main="", xlab="", ylab="Density")
  # y density
  par(mar=c(5,1.5,1.5,3)+0.1)
  max.density.y <- max(density(y)$y, na.rm=T)
  plot(density(y)$y, density(y)$x, type="l", xlim=c(0, 1.1*max.density.y), 
    ylim=ylim, main="", xlab="Density", ylab="")
}                              

## Example
#library(MSBVAR)
## generate some multivariate normal example data
#n <- 5000
#mu <- c(0, 2)
#vmat <- matrix(c(1, 0.7, 0.7, 1), nrow=2)
#
#dat <- rmultnorm(n, mu, vmat) # generates n multivariate normal obs.
#x <- dat[,1]
#y <- dat[,2]
#
#
#plot.2D.with.hist(x, y, xlab="This is the X variable", ylab="This is the Y variable.",
#  sub="Example Plot!")
## defining axis limits
#plot.2D.with.hist(x, y, xlab="This is the X variable", ylab="This is the Y variable.",
#  sub="Example Plot!", xlim=c(0,4), ylim=c(-2,2))
