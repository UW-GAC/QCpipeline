igvBatch <- function(chromosome, position, bam, outdir="./",
                     file="batch.igv", genome="hg19", collapse=FALSE) {
  stopifnot(length(chromosome) == length(position))

  con <- file(file, "w")
  out <-  c("new",
            paste("load", bam),
            paste("snapshotDirectory", outdir),
            paste("genome", genome))
  writeLines(out, con)

  for (i in 1:length(chromosome)) {
    out <- c(paste("goto chr", chromosome[i], ":", position[i], sep=""),
             "sort position")
    if (collapse) out <- c(out, "collapse")
    out <- c(out, paste("snapshot chr",chromosome[i], "_", position[i], ".png", sep=""))
    writeLines(out, con)
  }
  close(con)
}
