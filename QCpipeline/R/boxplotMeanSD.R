boxplotMeanSD <- function(x, y, data=NULL, xlab=NULL, ylab=NULL, nSD=1, ...)
{
  # get data
  if (!is.null(data))
  {
    if (is.null(xlab)) xlab <- x
    if (is.null(ylab)) ylab <- y
    x <- as.factor(data[, x])
    y <- data[, y]
  } else {
    x <- as.factor(x)
    y <- y
  }

  # parameters for plotting
  nlev <- length(unique(na.omit(x)))
  means <- tapply(y, x, function(x) mean(x, na.rm=T))
  stdev <- nSD * tapply(y, x, function(x) sd(x, na.rm=T))
  miny <- min(c(y, means-stdev), na.rm=T)
  maxy <- max(c(y, means+stdev), na.rm=T)
  maxy <- miny + 1.2 * (maxy-miny)

  # boxplots
  bp <- boxplot(y~x, plot=F)
  bxp(bp, xlim=c(0.5, nlev + 0.5), ylim=c(miny, maxy), xlab=xlab, ylab=ylab, ...)
  # add mean & SD
  par(new=TRUE)
  plotCI(x=means, uiw=stdev, col="red", xlim=c(0.5,nlev + 0.5), ylim=c(miny, maxy), xlab="", ylab="", xaxt="n")
  legend("top", c(paste("mean +/-", nSD, "SD"), "boxplot"), col = c("red", "black"), lty=c(1,1), lwd=1) 
}

