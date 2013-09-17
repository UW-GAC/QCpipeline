igvBatch <- function(chromosome, position, window=10,
                     bam, outdir="./", prefix="igv",
                     file="batch.igv", genome="hg19", collapse=FALSE) {
  stopifnot(length(chromosome) == length(position))

  con <- file(file, "w")
  out <-  c("new",
            paste("load", bam),
            paste("snapshotDirectory", outdir),
            paste("genome", genome))
  writeLines(out, con)

  for (i in 1:length(chromosome)) {
    out <- c(paste("goto chr", chromosome[i], ":", position[i]-window, "-",
                   position[i]+window, sep=""),
             "sort position")
    if (collapse) out <- c(out, "collapse")
    out <- c(out, paste("snapshot ", prefix, "_chr", chromosome[i], "_",
                        position[i], ".png", sep=""))
    writeLines(out, con)
  }
  close(con)
}
