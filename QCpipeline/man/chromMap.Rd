\name{chromMap}
\alias{chromToChar}
\alias{chromToInt}

\title{
Map between integer and character chromosome codes.
}
\description{
Map between integer and character chromosome codes.
}
\usage{
chromToChar(chrom.int)
chromToInt(chrom.char)
}
\arguments{
  \item{chrom.int}{vector of integer chromosome codes (1-28)}
  \item{chrom.char}{vector of character chromosome codes (1-22, X, XY, Y, M, U, XYY)}
}
%\details{
%}
\value{
\code{chromToChar} returns a vector of character chromosome codes.
\code{chromToInt} returns a vector of integer chromosome codes.
}
\author{
Stephanie Gogarten
}
%\seealso{}
\examples{
chromToChar(1:28)
chromToInt(c(1:22, "X", "XY", "Y", "M", "U", "XYY"))
}
