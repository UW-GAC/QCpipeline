\name{filterImputationSegments}
\alias{filterImputationSegments}
\alias{getAnomSegOverlap}

\title{Filter imputation segments}

\description{
Filter genotype probabilities in imputation segments
}

\usage{
filterImputationSegments(infile, outfile, anomalies, sample.annot, segments, chromosome, filteredValue = -1, block.size = 5000, overwrite = FALSE, verbose = TRUE)
getAnomSegOverlap(anomalies, segments)
}

\arguments{
  \item{infile}{input gprobs file to filter}
  \item{outfile}{output file (filtered gprobs)}
  \item{anomalies}{data frame of anomalies, which must include subjectID, chromosome, left.base, right.base}
  \item{sample.annot}{data frame of sample file from impute2, which must include subjectID}
  \item{segments}{data frame of imputation segments, which must include segment, chromosome, bp.start, bp.end}
  \item{chromosome}{chromosome to operate on}
  \item{filteredValue}{value to use for filtered genotype probabilities}
  \item{block.size}{number of SNPs to read in per block}
  \item{overwrite}{overwrite existing outfile?}
  \item{verbose}{report progress info?}
}

\value{
\code{filterImputationSegments} returns a data frame of filtered gprobs, with columns subjectID, chromosome, and segment.

\code{getAnomSegOverlap} returns a data frame of overlaps between anomalies and segments (using GRanges)
}

\author{
Adrienne Stilp
}


