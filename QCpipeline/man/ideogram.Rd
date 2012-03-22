\name{ideogram}
\alias{chromIntensityPlotIdeogram}
\alias{anomStatsPlotIdeogram}
\title{BAF/LRR plots with chromosome ideograms}
\description{Plot BAF/LRR with chromosome ideograms at the bottom.}
\usage{
chromIntensityPlotIdeogram(intenData, scan.ids, chrom.ids, ...,
                           ideo.zoom=TRUE, ideo.rect=FALSE,
                           cex.axis=1.7, cex.lab=1.7, cex.main=1.7,
                           cex.sub=1.5, cex.leg=1.5)

anomStatsPlotIdeogram(intenData, genoData, anom.stats,
                      snp.ineligible, win=5, ...,
                      ideo.zoom=FALSE, ideo.rect=TRUE,
                      cex.axis=1.7, cex.lab=1.7, cex.main=1.7,
                      cex.leg=1.5)
}
\arguments{
  \item{intenData}{\code{\link{IntensityData}} object}
  \item{genoData}{\code{\link{GenotypeData}} object}  
  \item{scan.ids}{ A vector containing the sample indices of the plots.}
  \item{chrom.ids}{ A vector containing the chromosome indices of the plots.}
  \item{anom.stats}{
    data.frame of chromosome anomalies with statistics, usually the output
    of \code{anomSegStats}.  Names must include "anom.id", "scanID", "chromosome",
    "left.index", "right.index", "method", "nmark.all", "nmark.elig", "left.base", "right.base",
    "nbase", "non.anom.baf.med", "non.anom.lrr.med", "anom.baf.dev.med",
    "anom.baf.dev.5", "anom.lrr.med", "nmark.baf", "nmark.lrr". Left and right refer
    to start and end, respectively, of the anomaly, in position order.
  }
  \item{snp.ineligible}{
    vector of ineligible snp ids (e.g., intensity-only, failed snps, XTR and HLA regions).
    See \code{\link{HLA}} and \code{\link{pseudoautosomal}}.
  }
  \item{win}{
    size of the window (a multiple of anomaly length) surrounding the anomaly to plot
  }
  
  \item{...}{additional arguments to \code{\link{chromIntensityPlot}} or
    \code{\link{anomStatsPlot}}}
  \item{ideo.zoom}{logical for whether to zoom in on the ideogram to
    match the range of the BAF/LRR plots}
  \item{ideo.zoom}{logical for whether to draw a rectangle on the
    ideogram indicating the range of the BAF/LRR plots}
  \item{cex.axis}{cex value for the axis}
  \item{cex.lab}{cex value for the labels}
  \item{cex.main}{cex value for the title}
  \item{cex.sub}{cex value for the subtitle}
  \item{cex.leg}{cex value for the ideogram legend}
}
\details{
  \code{chromIntensityPlotIdeogram} is a wrapper for
  \code{\link{chromIntensityPlot}}.
  \code{anomStatsPlotIdeogram} is a wrapper for \code{\link{anomStatsPlot}}.

  These functions call \code{\link{paintCytobands}} to draw ideograms.
}
\value{}
\author{Stephanie Gogarten}
\seealso{\code{\link{chromIntensityPlot}},
  \code{\link{anomStatsPlot}},
  \code{\link{paintCytobands}}
}
\examples{
library(GWASdata)
data(illuminaScanADF)
data(illuminaSnpADF)
blfile <- system.file("extdata", "illumina_bl.nc", package="GWASdata")
blnc <- NcdfIntensityReader(blfile)
intenData <-  IntensityData(blnc, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

genofile <- system.file("extdata", "illumina_geno.nc", package="GWASdata")
genonc <- NcdfGenotypeReader(genofile)
genoData <-  GenotypeData(genonc, scanAnnot=illuminaScanADF, snpAnnot=illuminaSnpADF)

# chromIntensityPlotIdeogram
scanID <- getScanID(illuminaScanADF, index=1:2)
#png("tmp_%03d.png", width=720, height=900)
chromIntensityPlotIdeogram(intenData=intenData, scan.ids=scanID,
                           chrom.ids=c(23,23),
                           colorGenotypes=TRUE, genoData=genoData)
#dev.off()

# anomStatsPlotIdeogram
scan.ids <- illuminaScanADF$scanID[1:2]
chrom.ids <- unique(illuminaSnpADF$chromosome)
snp.ids <- illuminaSnpADF$snpID[illuminaSnpADF$missing.n1 < 1]
snp.failed <- illuminaSnpADF$snpID[illuminaSnpADF$missing.n1 == 1]

# example results from anomDetectBAF
baf.anoms <- data.frame("scanID"=scan.ids[1],"chromosome"=21,
  "left.index"=100,"right.index"=200, sex="M", method="BAF",
  anom.id=1, stringsAsFactors=FALSE)

# example results from anomDetectLOH
loh.anoms <- data.frame("scanID"=scan.ids[2],"chromosome"=22,
  "left.index"=400,"right.index"=500, sex="F", method="LOH",
  anom.id=2, stringsAsFactors=FALSE)

anoms <- rbind(baf.anoms, loh.anoms)
data(centromeres.hg18)
stats <- anomSegStats(intenData, genoData, snp.ids=snp.ids, anom=anoms,
  centromere=centromeres.hg18)

#png("tmp_%03d.png", width=720, height=900)
anomStatsPlotIdeogram(intenData, genoData, anom.stats=stats,
  snp.ineligible=snp.failed, centromere=centromeres.hg18)
#dev.off()

close(genoData)
close(intenData)
}
\keyword{}
