## functions to filter pvalues from association tests


# function to find effective sample size MAF filter using quadratic equation
quadSolveMAF <- function(X, N) {
  sq <- sqrt(1 - (2 * X / N)) / 2
  pmin(0.5 + sq, 0.5 - sq)
}


calculateLambda <- function(stat, df){
  if (any(sum(stat < 0, na.rm=T))) stop("no negative values allowed in stat (does beta/se need to be squared?)")
  median(stat, na.rm=TRUE) / qchisq(0.5, df=df)
}


## filters - named list with 4 sets of filters (names are plot titles) (can have < 4)
qqPlotPng <- function(pval, stat, df, filters, outfile, ncol=2, addText="", main=NULL, ...) {
  nrow <- ceiling(length(filters) / ncol)
  png(outfile, width=360*ncol, height=360*nrow) # maybe need to calculate padding here too
  
  if (!is.null(main)){
    oma <- c(0,0,length(main)+0.1,0)
  } else{
    oma <- par()$oma
  }
  par(mfrow=c(nrow, ncol), mar=c(5,5,4,2)+0.1, lwd=1.5,
      cex.axis=1.5, cex.lab=1.5, cex.sub=1.5, cex.main=1.5, oma=oma)    
  for (i in 1:length(filters)) {
    filt <- filters[[i]]
    title <- names(filters)[i]
    lambda <- calculateLambda(stat[filt], df)
    subtitle <- paste("lambda =", format(lambda, digits=4, nsmall=3))
   
    if (length(filt) == 0){
      plot(1, type="n", axes=F, xlab="", ylab="", main=title, sub=subtitle)
      next
    }

    qqPlot(pval[filt], main=title, sub=subtitle)

    if (i == 1){
      mtext(side=3, line=-1, text=addText, padj=0.9, adj=0.02, outer=T, cex=1.5) 
    }
  }
  
  if (!is.null(main)){
    title(main, outer=T)
  }
  
  dev.off()
}

## filters - named list with 3 sets of filters (names are plot titles)
manhattanPlotPng <- function(pval, chromosome, filters, outfile, addText="", main=NULL, ...) {
  png(outfile, width=720, height=720)
  
  if (!is.null(main)){
    oma <- c(0,0,length(main)+0.1,0)
  } else{
    oma <- par()$oma
  }
  
  par(mfrow=c(length(filters),1), mar=c(5,5,4,2)+0.1, lwd=1.5, cex.lab=1.5, cex.main=1.5, oma=oma)
  for (i in 1:length(filters)) {
    filt <- filters[[i]]
    title <- names(filters)[i]

     if (length(filt) == 0){
      plot(1, type="n", axes=F, xlab="", ylab="", main=title)
      next
    }

    manhattanPlot(p=pval[filt], chromosome=chromosome[filt], main=title, ...)
    
    if (i == 1){
      mtext(side=3, line=-1, text=addText, padj=0.9, adj=0.02, outer=T, cex=1.5) 
    }
    
    if (!is.null(main)){
      title(main, outer=T)
    }
    
  }

  dev.off()
}
