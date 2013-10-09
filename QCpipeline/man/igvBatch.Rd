\name{igvBatch}
\alias{igvBatch}
\alias{igvPlot}
\title{Make IGV plot or batch file}
\description{Generate batch file for IGV to make snapshots of reads from
a BAM file, or capture the snapshots directly from within R.}
\usage{
igvBatch(chromosome, position, window=10, bam, outdir="./", prefix="igv",
         file="batch.igv", genome="hg19", collapse=c("yes", "no", "both"))
igvPlot(sock, chromosome, position, window=10, bam, outdir="./", prefix="igv",
        genome="hg19", collapse=c("yes", "no", "both"))
}
\arguments{
  \item{chromosome}{Vector of chromosomes}
  \item{position}{Vector of positions}
  \item{window}{Integer of number of bases to display on either side of \code{position}}
  \item{bam}{Vector of BAM files to load}
  \item{outdir}{Output directory for snapshot png files}
  \item{prefix}{Prefix for snapshot png files}
  \item{file}{Output batch file}
  \item{genome}{Genome to load}
  \item{collapse}{"yes"=collapse the reads, "no"=expand the reads,
    "both"=make two plots, one collapsed and one expanded.}
  \item{sock}{A socket connection to IGV.}
}
\details{
}
\value{
  \code{igvBatch} creates a text file which be run in IGV with the command "File->Run Batch Script."
  The result of running this file will be .png files for each
  chromosome/position pair.

  \code{igvPlot} takes snapshot .png files for each chromosome pair
  directly on a running IGV session.
}
\author{Stephanie Gogarten}
\seealso{}
\examples{
bamfile <- "/projects/geneva/gcc-fs2/Peters_seq/bam/B00Y2ACXX_1_GGCTAC.bam"
outdir <- tempdir()
outfile <- tempfile()
chrom <- c(1, 1)
pos <- c(155265446, 224804891)
igvBatch(chrom, pos, bam=bamfile, outdir=outdir, file=outfile,
         genome="hg19", collapse="yes")
readLines(outfile)
unlink(outfile)

\dontrun{
library(SRAdb)
sock <- IGVsocket()
igvPlot(sock, chrom, pos, bam=bamfile, outdir=outdir,
        genome="hg19", collapse="yes")
}
unlink(outdir, recursive=TRUE)
}
\keyword{}
