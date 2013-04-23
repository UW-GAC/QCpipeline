\name{GenotypeReader}
\alias{GenotypeReader}
\title{Return a NcdfGenotypeReader or GdsGenotypeReader}
\description{
  Return a NcdfGenotypeReader or GdsGenotypeReader based on the file extension.
}
\usage{
GenotypeReader(file)
}
\arguments{
  \item{file}{Character string with the path to a NetCDF file (ending in
    ".nc") or a GDS file (ending in ".gds").}    
}
\details{
}
\value{
  A \code{\link{NcdfGenotypeReader}} or \code{\link{GdsGenotypeReader}} object.
}
\author{Stephanie Gogarten}
\seealso{}
\examples{
ncfile <- "tmp.nc"
simulateGenotypeMatrix(ncdf.filename=ncfile)
dat <- GenotypeReader(ncfile)
stopifnot(is(dat, "NcdfGenotypeReader"))
close(dat)
gdsfile <- "tmp.gds"
convertNcdfGds(ncfile, gdsfile)
dat <- GenotypeReader(gdsfile)
stopifnot(is(dat, "GdsGenotypeReader"))
close(dat)
unlink(c(ncfile, gdsfile))
}
\keyword{}
