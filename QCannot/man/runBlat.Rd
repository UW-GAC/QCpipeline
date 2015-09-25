\name{runBlat}
\alias{runBlat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run BLAT search on input DNA sequence(s)
}
\description{
This function runs the command line BLAT tool to find region(s) of the
genome with high similarity to the user's input sequence(s).
}
\usage{
runBlat(seq.input, make.fasta=TRUE, output.fn="blat_query", ooc=TRUE,
tileSize=11, stepSize=tileSize, repMatch=1024, minScore=30, minIdentity=90,
find.best=TRUE, out.type="psl", return.result=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq.input}{One of two objects: (1) a dataframe where the first
    column is the sequence identifier and the second column is a
    character string of the DNA sequence to be searched OR (2) if the user already has a FASTA file of input sequences created, \code{seq.input} can be this FASTA file name (full file path).
 }
  \item{make.fasta}{Logical flag on whether the \code{seq.input} object needs
    to be converted into FASTA format. Default=TRUE
 }
 \item{output.fn}{Base name for all output files. Defaults to \dQuote{blat_query.}
 }
 \item{ooc}{Logical flag for whether to speed up the BLAT search by providing a file of commonly repeated sequences. Default=TRUE
}
 \item{tileSize}{Defaults to 11.
 }
 \item{stepSize}{Defaults to tileSize. To match web BLAT query, set to 5.
 }
 \item{repMatch}{Defaults to 1024. To match web BLAT query, set to 2253.
 }
 \item{minScore}{Defaults to 30. To match web BLAT query, lower to 20. For very permissive searching, can set to 0.
 }
 \item{minIdentity}{Default to 90. For very permissive searching, can set to 0.
 }
 \item{find.best}{Logical flag on whether to apply \dQuote{near best in
     genome} algorithm (pslReps) to full set of BLAT searches. Default=TRUE. Will only work if \code{out.type="psl"}.
 }
 \item{out.type}{BLAT output file type. Defaults to \dQuote{psl.}
 }
 \item{return.result}{Logical flag whether to return BLAT results as R object. Default=FALSE. Will only work if output file type is set to \dQuote{psl}.
} 
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  If return.result=TRUE, returns R data frame of all BLAT results; if
  find.best=TRUE, also returns a logical flag  \dQuote{top.hit} indicating the
  results of the pslReps program's attempt to find
  \dQuote{near best in genome} matches.
  If return.result=FALSE, see output.fn* files written in current
  working directory.
  In the output, \dQuote{Q_name} refers to \dQuote{query name} and
  contains the name for each searched sequence provided by the user.
}
\note{
  This function runs command line tools and is not currently set up to
  be run on Windows machines. If \code{find.best=TRUE}, the pslReps program is
  run after the BLAT tool. 
  The \dQuote{score} and \dQuote{pct.identity} fields are produced when using the
  web-based BLAT tool but not with the command line tool. This function
  calculates these two metrics based on formulae provided by BLAT author
  Jim Kent.
}
\references{
  For more information on the BLAT settings, please see the BLAT
  online documentation:
  \url{http://genome.ucsc.edu/goldenPath/help/blatSpec.html},
  \url{http://genome.ucsc.edu/FAQ/FAQblat.html}

  For more information on the ooc setting, please see:
  \url{http://genome.ucsc.edu/FAQ/FAQblat.html#blat6}
  
  For more information on the pslReps program, please see:
  \url{https://genome.ucsc.edu/goldenPath/help/blatSpec.html#pslRepsUsage}

  For a description of the FASTA file format, please see:
  \url{http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml}

}
\author{
Nelson, S.C.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
library(QCannot)
### for running on an existing FASTA file
fasta.input <- system.file("extdata", "seq_input_fasta_example.fa",
package="QCannot")
# preview fasta file
system(paste("cat", fasta.input))
## >rs1000073
## TCCAGCAGGCACCTAGCACTATTTTCTAGTGCTCCATCTCTTAGCAGGGACTCTGTTCAG[A/G]TCTTTGTTTTCTGTG
## GGTGGCAATGAGGAAATGAGCTCACCTGAGCCATGAGGGAAACTA
## >rs1000313
## TTGCTGTTCCTAAACTGGCCTTTGTGCAATGCATAGTGAAAGCCGCTCACCACTGGCTTA[A/G]GAAACTCCCGGGAGG
## CAGTTGCAAAGCTTTTTGAGGGAGATCGGAGGTGGGACAGTTTGC

# run BLAT function
blat.rslt <- runBlat(seq.input=fasta.input, make.fasta=FALSE, return.result=TRUE)
dim(blat.rslt); blat.rslt
##   match mismatch repmatch Ns Q_gap_count Q_gap_bases T_gap_count T_gap_bases
## 1   121        0        0  0           1           1           0           0
## 2   121        0        0  0           1           1           0           0
##   strand    Q_name Q_size Q_start Q_end T_name    T_size   T_start     T_end
## 1      + rs1000073    122       0   122   chr1 249250621 157255335 157255456
## 2      + rs1000313    122       0   122   chr1 249250621  15405428  15405549
##   block_count blocksizes qStarts              tStarts score pct.identity id
## 1           2     61,60,   0,62, 157255335,157255396,   120         99.2  1
## 2           2     62,59,   0,63,   15405428,15405490,   120         99.2  2
##   top.hit
## 1    TRUE
## 2    TRUE

### for running on non-FASTA input
seq.input <- system.file("extdata", "seq_input.txt", package=QCannot)
# preview input file
seq.input
##         id
## 1 rs1000073
## 2 rs1000313
##                                                                                                                             seq
## 1 TCCAGCAGGCACCTAGCACTATTTTCTAGTGCTCCATCTCTTAGCAGGGACTCTGTTCAG[A/G]TCTTTGTTTTCTGTGGGTGGCAATGAGGAAATGAGCTCACCTGAGCCATGAGGGAAACTA
## 2 TTGCTGTTCCTAAACTGGCCTTTGTGCAATGCATAGTGAAAGCCGCTCACCACTGGCTTA[A/G]GAAACTCCCGGGAGGCAGTTGCAAAGCTTTTTGAGGGAGATCGGAGGTGGGACAGTTTGC

blat.rslt <- runBlat(seq.input=seq.input, return.result=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
