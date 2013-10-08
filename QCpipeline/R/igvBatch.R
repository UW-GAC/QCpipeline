igvBatch <- function(chromosome, position, window=10,
                     bam, outdir="./", prefix="igv",
                     file="batch.igv", genome="hg19",
                     collapse=c("yes", "no", "both")) {
  stopifnot(length(chromosome) == length(position))
  collapse <- match.arg(collapse)

  con <- file(file, "w")
  out <-  c("new",
            paste("genome", genome),
            paste("load", bam),
            paste("snapshotDirectory", outdir))
  writeLines(out, con)

  for (i in 1:length(chromosome)) {
    out <- c(paste("goto chr", chromosome[i], ":", position[i]-window, "-",
                   position[i]+window, sep=""),
             "sort position")
    if (collapse %in% c("no", "both")) {
      out <- c(out, paste("snapshot ", prefix, "_chr", chromosome[i], "_",
                          position[i], "_expand.png", sep=""))
    }
    if (collapse %in% c("yes", "both")) {
      out <- c(out, "collapse")
      out <- c(out, paste("snapshot ", prefix, "_chr", chromosome[i], "_",
                          position[i], "_collapse.png", sep=""))
    }
    writeLines(out, con)
  }
  close(con)
}


igvPlot <- function(sock, chromosome, position, window=10,
                     bam, outdir="./", prefix="igv",
                     genome="hg19",
                     collapse=c("yes", "no", "both")) {
  stopifnot(length(chromosome) == length(position))
  collapse <- match.arg(collapse)

  IGVclear(sock)
  IGVgenome(sock, genome)
  IGVload(sock, bam)
  for (i in 1:length(chromosome)) {
    IGVgoto(sock, paste("chr", chromosome[i], ":", position[i]-window, "-",
                        position[i]+window, sep=""))
    IGVsort(sock, "position")
    if (collapse %in% c("no", "both")) {
      IGVsnapshot(sock,
                  fname=paste(prefix, "_chr", chromosome[i], "_",
                    position[i], "_expand.png", sep=""),
                  dirname=outdir)
    }
    if (collapse %in% c("yes", "both")) {
     IGVcollapse(sock)
     IGVsnapshot(sock,
                 fname=paste(prefix, "_chr", chromosome[i], "_",
                   position[i], "_collapse.png", sep=""),
                 dirname=outdir)
    }
  }
}
