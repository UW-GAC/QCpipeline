\name{getRefStrand}
\alias{getRefStrand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine +/- strand orientation for Illumina manifest, using BLAT results
}
\description{
This function integrates an Illumina manifest with BLAT results produced
with the companion function \code{runBlat} to determine +/- strand
orientation.
}
\usage{
getRefStrand(snp.dat, blat.rslt, metric.sel="score")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snp.dat}{An R object created by reading in an Illumina
    manifest. Required columns are \dQuote{Name,} \dQuote{Chr,} \dQuote{MapInfo,} and
    \dQuote{IlmnStrand.} Assumed that the \dQuote{Chr} field contains the
    following character codes for non-autosomes: X, XY, Y, and MT.
 }
 \item{blat.rslt}{An R object created by reading in BLAT search results
   produced by the companion function \code{runBlat}.
 }
 \item{metric.sel}{BLAT metric to use when selecting among multiple BLAT
 matches in the right location. Options are \dQuote{score} or
 \dQuote{pct.identity.} Default is \dQuote{score.}
 }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns a data frame with the same row dimension as the \code{snp.dat}
  input. The following columns are added:
  \item{blat.strand}{+ or - strand orientation of the searched sequence}
  \item{blat.id}{integer ID of the BLAT result}
  \item{RefStrand}{+ or - orientation of the design alleles}
}
\note{
  BLAT results are aligned to the Illumina manifest first by limiting to
  BLAT matches with the same chromosome as the input variant, secondly
  by limiting to BLAT matches that span the input variant's base pair
  position. If multiple matches remain, the BLAT result with the
  highest metric (user-specified with the \code{metric.sel} argument) is chosen.
  Input variants without a BLAT match at the right chromosome and
  position will be returned in the output but with \dQuote{NA} values
  for the added columns.

  RefStrand is determined using IlmnStrand (TOP/BOT designation of
  design alleles) and BLAT strand (+/- designation of searched BLAT
  sequence, assumed to be TopGenomicSeq).
  
  Assumes that genome build is consistent between Illumina manifest and
  BLAT search.

  Assumes that BLAT search was run on TopGenomicSeq.

}
%\references{
% }
\author{
Nelson, S.C.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[QCannot]{runBlat}}
}

% \examples{
% \dontrun{
% library(QCannot)
% fasta.input <- system.file("seq_input_fasta_example.fa", package=QCannot)
% blat.rslt <- run.blat(seq.input=fast.input, make.fasta=FALSE,
% return.result=TRUE)
% dim(blat.rslt); blat.rslt
% # [1]  1 25
% #   match mismatch repmatch Ns Q_gap_count Q_gap_bases T_gap_count T_gap_bases
% # 1   201        0        0  0           0           0           0           0
% #   strand                            Q_name Q_size Q_start Q_end T_name
% # 1      + hg19_ct_CHRNA5track_4165_78858078    201       0   201  chr15
% #      T_size  T_start    T_end block_count blocksizes qStarts   tStarts score
% # 1 102531392 78857977 78858178           1       201,      0, 78857977,   201
% #   pct.identity id top.hit
% # 1          100  1    TRUE
% }
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{manip}
