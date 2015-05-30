
# function to return either GenotypeDataByChr or GenotypeData object given path
GenotypeDataFromPath <- function(path) {
    ## check if directory
    if (file_test("-d", path)) {
        GenotypeDataByChr(path)
    } else {
        GenotypeData(GdsGenotypeReader(path))
    }
}
