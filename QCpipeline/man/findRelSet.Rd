\name{findRelSet}
\alias{findRelSet}
\title{Partition a Sample into an 'Unrelated Set' and a 'Related Set'}
\description{\code{findRelSet} is used to partition a sample from a genetic study into an 'unrelated set' and a 'related set'.  The 'unrelated set' contains individuals who are all mutually unrelated to each other, and the 'related set' contains individuals who are related to someone in the 'unrealted set'.}
\usage{
findRelSet(kinMat, thresh, unrel.set = NULL)
}
\arguments{
  \item{kinMat}{A symmetric matrix of pairwise kinship coefficients for every pair of individuals in the sample (the values on the diagonal do not matter, but the upper and lower triangles must both be filled) used for partitioning the sample into the 'unrelated' and 'related' sets.  See 'Details' for how this interacts with \code{thresh} and \code{unrel.set}.}
  \item{thresh}{Threshold value on \code{kinMat} used for declaring each pair of individuals as related, unrelated, or divergent.  The default value is 0.025.  See 'Details' for how this interacts with \code{kinMat}.}
  \item{unrel.set}{An optional vector of indices for identifying individuals that are forced into the unrelated set.  See 'Details' for how this interacts with \code{kinMat}.}
}
\details{We recommend using the software KING to estimate pairwise kinship coefficients to be used in \code{kinMat}.  Any pair of individuals with a pairwise kinship greater than \code{thresh} will be declared 'related' and any pair of individuals with a pairwise kinship less than \code{-thresh} will be declared 'divergent'.  Typically, \code{thresh} is set to be the amount of error around 0 expected in the kinship estimate for a pair of truly unrelated individuals.  If \code{unrel.set = NULL}, the PC-AiR algorithm is used to find an 'optimal' partition (see 'References' for a paper describing the algorithm).  If \code{unrel.set} is specified, then all individuals indexed in \code{unrel.set} are forced in the 'unrelated set' and the PC-AiR algorithm is used to partition the rest of the sample; this is especially useful for including reference samples of known ancestry in the 'unrelated set'.
}
\value{A vector of indices identifying which individuals are in the 'related set'.  Note: these are not individual IDs, they are indices for which rows/columns of \code{kinMat} are in the 'related set'.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Matthew P. Conomos}
\note{\code{findRelSet} is called internally in the function \code{pcair} but may also be used on its own to partition the sample into 'unrelated' and 'related' sets without performing PCA.}
%\seealso{
% 	\code{\link{pcair}} which uses this function for finding principal components in the presence of related individuals.
% 	\code{\link{king2mat}} for creating a matrix of pairwise kinship coefficient estimates that can be used as \code{kinMat} from KING output text files.
%}
\examples{
\dontrun{
	# load matrix of estimated kinship coefficients
	data("HapMap_ASW_MXL_kinMat")
	# find related set
	myrels <- findRelSet(HapMap_ASW_MXL_kinMat, thresh = 0.025)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
