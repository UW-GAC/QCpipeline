\name{igvBatch}
\alias{igvBatch}
\title{Generate batch file for IGV}
\description{Generate batch file for IGV to make snapshots of reads from
a BAM file.}
\usage{
igvBatch(chromosome, position, bam, outdir="./",
         file="batch.igv", genome="hg19", collapse=FALSE)
}
\arguments{
  \item{chromosome}{Vector of chromosomes}
  \item{position}{Vector of positions}
  \item{bam}{BAM file to load}
  \item{outdir}{Output directory for snapshot png files}
  \item{file}{Output batch file}
  \item{genome}{Genome to load}
  \item{collapse}{Logical for whether to collapse the reads}
}
\details{
}
\value{
  A text file which be run in IGV with the command "File->Run Batch Script."
  The result of running this file will be .png files for each
  chromosome/position pair.
}
\author{Stephanie Gogarten}
\seealso{}
\examples{
bamfile <- "/projects/geneva/gcc-fs2/Peters_seq/bam/B00Y2ACXX_1_GGCTAC.bam"
outdir <- tempdir()
outfile <- tempfile()
chrom <- c(1, 1)
pos <- c(155265446, 224804891)
igvBatch(chrom, pos, bamfile, outdir, file=outfile,
         genome="hg19", collapse=TRUE)
readLines(outfile)
file.remove(outfile)
}
\keyword{}
