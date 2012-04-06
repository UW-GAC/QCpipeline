chromIntensityPlotIdeogram <- function(intenData, scan.ids, chrom.ids,
                                       main=NULL, info=NULL, ...,
                                       ideo.zoom=TRUE, ideo.rect=FALSE,
                                       cex.axis=1.7, cex.lab=1.7,
                                       cex.main=1.7, cex.sub=1.5, cex.leg=1.5) {

  if (length(scan.ids) != length(chrom.ids)) {
    stop("scan.ids and chrom.ids must be parallel vectors of the same length")
  }

  snpID <- getSnpID(intenData)
  chrom <- getChromosome(intenData)
  pos <- getPosition(intenData)
  
  chrom.char <- chrom.ids
  chrom.char[chrom.ids == XchromCode(intenData)] <- "X"
  chrom.char[chrom.ids == YchromCode(intenData)] <- "Y"
  
  layout(matrix(c(1,2,3), nrow=3, ncol=1), heights=c(0.4, 0.4, 0.2))
  for (i in 1:length(scan.ids)) {
    par(mar=c(5,4,4,2)+0.1, mgp=c(2.5,0.75,0))
    par(cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, cex.sub=cex.sub)
    chromIntensityPlot(intenData, scan.ids=scan.ids[i], chrom.ids=chrom.ids[i],
                       type="LRR", main=main[i], info=info[i], ...)
    chromIntensityPlot(intenData, scan.ids=scan.ids[i], chrom.ids=chrom.ids[i],
                       type="BAF", main=main[i], info=info[i], ...)
    par(mar=c(1,4,1,2)+0.1)
    posc <- pos[chrom == chrom.ids[i]]
    if (ideo.zoom) {
      ideo.x <- c(min(posc), max(posc))
    } else {
      ideo.x <- c(0,lengthChromosome(chrom.char[i],"bases"))
    }
    plot(ideo.x, c(-2,2),
         type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    paintCytobands(chrom.char[i], units="bases", width=1, cex.leg=cex.leg)
    if (ideo.rect) rect(min(posc), -1.2, max(posc), 0.2, border="red", lwd=2)
  }
  
}


anomStatsPlotIdeogram <- function(intenData, genoData, anom.stats, snp.ineligible,
                                  win=5, main=NULL, info=NULL, ...,
                                  ideo.zoom=FALSE, ideo.rect=TRUE,
                                  cex.axis=1.7, cex.lab=1.7,
                                  cex.main=1.7,  cex.leg=1.5) {
  chrom <- getChromosome(intenData)
  pos <- getPosition(intenData)

  chrom.ids <- anom.stats$chromosome
  chrom.char <- chrom.ids
  chrom.char[chrom.ids == XchromCode(intenData)] <- "X"
  chrom.char[chrom.ids == YchromCode(intenData)] <- "Y"
  
  layout(matrix(c(1,2,3), nrow=3, ncol=1), heights=c(0.4, 0.4, 0.2))
  for (i in 1:nrow(anom.stats)) {
    par(mar=c(5,4,4,2)+0.1, mgp=c(2.5,0.75,0))
    par(cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab)
    anomStatsPlot(intenData, genoData, anom.stats[i,], snp.ineligible=snp.ineligible,
                  type="LRR", win=win, main=main[i], info=info[i], ...)
    anomStatsPlot(intenData, genoData, anom.stats[i,], snp.ineligible=snp.ineligible,
                  type="BAF", win=win, main=main[i], info=info[i], ...)
    
    par(mar=c(1,4,1,2)+0.1)
    posc <- pos[chrom == chrom.ids[i]]
    n <- length(posc)
    xlim.all <- c(posc[1], posc[n])
    leftx <- pos[anom.stats$left.index[i]]
    rightx <- pos[anom.stats$right.index[i]]
    dif <- rightx-leftx
    left <- max(leftx-win*dif, xlim.all[1]); right <- min(rightx+win*dif, xlim.all[2])
    if (ideo.zoom) {
      ideo.x <- c(left, right)
    } else {
      ideo.x <- c(0,lengthChromosome(chrom.char[i],"bases"))
    }
    plot(ideo.x, c(-2,2),
         type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    paintCytobands(chrom.char[i], units="bases", width=1, cex.leg=cex.leg)
    if (ideo.rect) rect(left, -1.2, right, 0.2, border="red", lwd=2)
  }
  
}
