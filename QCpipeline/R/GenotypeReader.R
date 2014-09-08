## return either a NcdfGenotypeReader or GdsGenotypeReader based on file extension

GenotypeReader <- function(file) {
  ext <- file_ext(file)
  if (ext == "gds") {
    GdsGenotypeReader(file)
  } else if (ext == "nc") {
    NcdfGenotypeReader(file)
  } else {
    stop("genotype file must end in '.gds' or '.nc'")
  }
}

## return either a NcdfIntensityReader or GdsIntensityReader based on file extension

IntensityReader <- function(file) {
  ext <- file_ext(file)
  if (ext == "gds") {
    GdsIntensityReader(file)
  } else if (ext == "nc") {
    NcdfIntensityReader(file)
  } else {
    stop("genotype file must end in '.gds' or '.nc'")
  }
}
