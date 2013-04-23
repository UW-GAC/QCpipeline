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
