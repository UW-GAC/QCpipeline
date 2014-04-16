## functions to filter pvalues from association tests


# function to find effective sample size MAF filter using quadratic equation
quadSolveMAF <- function(X, N) {
  sq <- sqrt(1 - (2 * X / N)) / 2
  pmin(0.5 + sq, 0.5 - sq)
}



## filters - named list with 4 sets of filters (names are plot titles) (can have < 4)
qqPlotPng <- function(pval, filters, outfile) {
  png(outfile, width=720, height=720)
  par(mfrow=c(2,2), mar=c(5,5,4,2)+0.1, lwd=1.5,
      cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5)    
  for (i in 1:length(filters)) {
    filt <- filters[[i]]
    title <- names(filters)[i]
    lambda <- median(-2*log(pval[filt]), na.rm=TRUE) / 1.39
    subtitle <- paste("lambda =", format(lambda, digits=4, nsmall=3))
    qqPlot(pval[filt], truncate=FALSE, main=title, sub=subtitle)
  }
  dev.off()
}


## filters - named list with 3 sets of filters (names are plot titles)
manhattanPlotPng <- function(pval, chromosome, filters, outfile) {
  png(outfile, width=720, height=720)
  par(mfrow=c(length(filters),1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5)
  for (i in 1:length(filters)) {
    filt <- filters[[i]]
    title <- names(filters)[i]
    manhattanPlot(p=pval[filt], chromosome=chromosome[filt], main=title)
  }
  dev.off()
}
