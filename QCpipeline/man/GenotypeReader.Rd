\name{GenotypeReader}
\alias{GenotypeReader}
\alias{IntensityReader}
\title{Return a NcdfGenotypeReader or GdsGenotypeReader}
\description{
  Return a NcdfGenotypeReader or GdsGenotypeReader based on the file extension.
}
\usage{
GenotypeReader(file)
IntensityReader(file)
}
\arguments{
  \item{file}{Character string with the path to a NetCDF file (ending in
    ".nc") or a GDS file (ending in ".gds").}    
}
%\details{}
\value{
  \code{GenotypeReader}: a \code{\link{NcdfGenotypeReader}} or \code{\link{GdsGenotypeReader}} object.
  \code{IntensityReader}: a \code{\link{NcdfIntensityReader}} or \code{\link{GdsIntensityReader}} object.
}
\author{Stephanie Gogarten}
%\seealso{}
\examples{
ncfile <- "tmp.nc"
simulateGenotypeMatrix(filename=ncfile, file.type="ncdf")
dat <- GenotypeReader(ncfile)
stopifnot(is(dat, "NcdfGenotypeReader"))
close(dat)
gdsfile <- "tmp.gds"
simulateGenotypeMatrix(filename=gdsfile, file.type="gds")
dat <- GenotypeReader(gdsfile)
stopifnot(is(dat, "GdsGenotypeReader"))
close(dat)
unlink(c(ncfile, gdsfile))
}
