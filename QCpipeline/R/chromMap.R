chromToChar <- function(chrom.int) {
    chrom.map <- c(setNames(1:22, 1:22), "23"="X", "24"="XY", "25"="Y", "26"="M", "27"="U", "28"="XYY")
    unname(chrom.map[chrom.int])
}

chromToInt <- function(chrom.char) {
    chrom.map <- c(setNames(1:22, 1:22), "X"=23, "XY"=24, "Y"=25, "M"=26, "U"=27, "XYY"=28)
    unname(as.integer(chrom.map[chrom.char]))
}
