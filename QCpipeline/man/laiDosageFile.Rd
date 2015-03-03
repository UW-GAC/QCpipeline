\name{laiDosageFile}
\alias{laiDosageFile}
\alias{checkLaiDosageFile}

\title{
Create LAI dosage GDS files
}

\description{
These functions create a GDS file from an RFMix local ancestry inference output file.
}
\usage{
laiDosageFile(lafile, outfile, pops, blockSize = 5000, onlyUnique = FALSE, dryRun = FALSE)
checkLaiDosageFile(lafile, genoDataList, pops, snps, blockSize=5000, onlyUnique=FALSE)
}


\arguments{
  \item{lafile}{input LAI file}
  \item{outfile}{output gds file}
  \item{pops}{named vector mapping populations; see details.}
  \item{blockSize}{number of lines to read at a time}
  \item{onlyUnique}{logical indicator to only write unique variants; see details.}
  \item{dryRun}{logical indicating whether to write the gds file or just to process \code{lafile}}
  \item{genoDataList}{named list of genoData objects}
  \item{snps}{data frame of SNP annotation; output from laiDosageFile}
}


\details{
\code{laiDosageFile} creates a single GDS file with the dosages of each population at each SNP/sample combination  stored as nodes named \code{"dosage_<pop>"}. NOTE: Only tested for 3 populations but will probably work with different numbers of populations.

\code{checkLaiDosageFile} checks the output from \code{laiDosageFile} against the RFMix file.

\code{pops} should be a named vector indicating the coding in the file, e.g. \code{pops=c("eur"=0, "afr"=1)}.

If \code{onlyUnique=FALSE}, all ancestries will be written to the gds file. If \code{onlyUnique=TRUE}, only ancestries that are not identical to the previous ancestry will be written to the gds file.

}


\value{
In addition to creating the GDS file, a snp annotation is also returned with columns:
\itemize{
\item{snpID: integer identifier}
\item{chromosome: chromosome}
\item{position: position}
\item{rsID: snp name in \code{laifile}}
}


If \code{onlyUnique=TRUE}, then two additional columns are returned:
\itemize{
\item{pos.start: first position with these ancestry values}
\item{pos.end: last position with these ancestry values}
}
Additionally, \code{position} is the average of \code{pos.start} and \code{pos.end}.

}
\author{Adrienne Stilp}

\examples{
library(gdsfmt)
lafile <- tempfile()

## make some dummy local ancestry data
dat <- data.frame(CHROM=rep(c(1,2), each=10), POS=rep(1:10, times=2), ID=paste0("rs", 1:20), REF="-", ALT="-", QUAL="-", FILTER="-", INFO="-", FORMAT="-")
dat$zz100 <- c(rep("2/0", 5), rep("0/2", 11), rep("1/1", 4))
dat$zz101 <- c(rep("2/0", 19), "1/1")
dat$zz102 <- "1/0"
names(dat) <- gsub("^zz", "", names(dat))
gz1 <- gzfile(lafile, "w")
write.table(dat, file=gz1, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, na="")
close(gz1)
gz <- file(lafile, "r")
lines <- readLines(gz)
lines <- c("##file format: similar to VCFv4.1",
           "##fileDate=2014-11-28",
           "##description: 0=European ancestry, 1=African ancestry, 2=Native American ancestry",
           "#CHROM  POS  ID  REF	ALT	QUAL	FILTER	INFO	FORMAT	100	101	102", lines)
close(gz)
  
gz1 <- gzfile(lafile, "w")
writeLines(lines, gz1)
close(gz1)


## setup for running the function
# For example an entry of 1/0 corresponds to a european allele dosage of 1 and an african allele dosage of 1.
pops <- c("eur"=0, "afr"=1, "amer"=2)
 
gdsfile <- tempfile()
snps <- laiDosageFile(lafile, gdsfile, pops)
snps
 
## reading from the output file
# open for reading
gds <- openfn.gds(gdsfile)
gds.eur <- GdsGenotypeReader(gds, genotypeVar="dosage_eur")
gds.afr <- GdsGenotypeReader(gds, genotypeVar="dosage_afr")
gds.amer <- GdsGenotypeReader(gds, genotypeVar="dosage_amer")

head(getGenotype(gds.eur))
head(getGenotype(gds.afr))
head(getGenotype(gds.amer))

genoDataList <- list(eur=GenotypeData(gds.eur), afr=GenotypeData(gds.afr), amer=GenotypeData(gds.amer))
checkLaiDosageFile(lafile, genoDataList, pops, snps, onlyUnique=FALSE)

close(gds.eur)
unlink(gdsfile)

gdsfile <- tempfile()
snps.unique <- laiDosageFile(lafile, gdsfile, pops, onlyUnique=TRUE)
snps.unique
unlink(gdsfile)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
