\name{ibdAssignRelatedness}
\alias{ibdAssignRelatednessPCRelate}

\title{
Assign relationships
}

\description{
  \code{ibdAssignRelatednessPCRelate} identifies observed relatives using the kinship coefficients and k0 estimates from the PC-Relate model.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usage{
ibdAssignRelatednessPCRelate(k0, kc, cut.kc.dup=1/(2^(3/2)),
                         cut.kc.fs=1/(2^(5/2)), cut.kc.deg2=1/(2^(7/2)),
                         cut.kc.deg3=1/(2^(9/2)), cut.k0=0.1)
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\arguments{
  
  \item{k0}{A vector of k0 values.}
  
  \item{kc}{A vector of kinship coefficient values.}

\item{cut.kc.dup}{Kinship coefficient threshold for dividing duplicates and first degree relatives.}
\item{cut.kc.fs}{Kinship coefficient threshold for dividing full siblings and second degree relatives.}
\item{cut.kc.deg2}{Kinship coefficient threshold for dividing second and third degree relatives.}
\item{cut.kc.deg3}{Kinship coefficient threshold for dividing third degree relatives and unrelated.}
\item{cut.k0}{k0 threshold for dividing parent-offsprings pairs from other relatives.}


}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%\details{
%}

\value{
  \code{ibdAssignRelatednessPCRelate} and returns a vector of relationships with
  values "Dup"=duplicate, "PO"=parent-offspring, "FS"=full sibling,
  "Deg2"=second degree, "Deg3"=third degree, "U"=unrelated, and
  "Q"=unknown.
}

\author{Cecelia Laurie, Stephanie Gogarten}

\seealso{
  \code{\link[GENESIS]{pcrelate}},
  \code{\link{pedigreePairwiseRelatedness}}
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\examples{
k0 <- c(0, 0, 0.25, 0.5, 0.75, 1)
kc <- c(0.5, 0.25, 0.25, 0.125, 0.063, 0)
obs.rel <- ibdAssignRelatednessPCRelate(k0, kc)
}


\keyword{manip}
