\name{expRelsCategory}
\alias{expRelsCategory}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\title{
Assign relatedness from pedigree data: relations and categorical assignments
}
\description{
This function assigns relationships from pedigree data. Output includes the theoretical pairwise kinship coefficients.
 Categories of Dup (duplicate), PO (parent-offspring), FS (full sibling), Deg1, Deg2, Deg3, and U are also assigned.
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
expRelsCategory(pedigree)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identifiers for family, individual, individual's mother, individual's father and individual's sex (coded as "M" or "F") . 
     Identifiers can be integer, numeric or character but identifiers for mother and father for founders are assumed to be 0.  Error messages are returned for pedigree inconsistencies. See \code{\link{pedigreeCheck}}   }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{  Assigns specific relationships between individuals in a pedigree, including "U" = unrelated, "PO" = parent/offspring, "FS" = full siblings, "HS" = half siblings, "Av" = avuncular,  "FC" = first cousins, "GpGc" = grandparent-grandchild, "GGp" = great-grandparent-great-grandchild, "GAv" = grand-avuncular, "HSFC" = half-sib-first-cousin, "DFC" = double first cousin, among others.
 In addition, categories of relatedness (Deg0, PO, FS, Deg1, Deg2, Deg3, and U) are assigned based on thresholds given in Table 1 of Manichaikul, et al (2010).
 Here Deg0 is assigned (rather than Dup) to pairs with kinship > the Dup threshold (0.354) - these can occur with close inbreeding.
 Here Deg1 is assigned if pair is not PO or FS but kinship is greated than threshold (0.177) which can happen with inbreeding but also, e.g., half-sibs first-cousing.
Specific relationships are not determined calculated for inbred families but categories are assigned.
Output includes the theoretical pairwise kinship coefficients.
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{

A list with the following components:

\item{inbred.fam}{A vector of id's of families with inbreeding. }
 
\item{relprs.all}{A dataframe with columns "Individ1", "Individ2", "relation", "kinship", "family", and "exp.rel" containing the id's of the pair of individuals, the specific relationship between the individuals (see Details section), kinship coefficient, and assigned category of the relationship (see Details).
   For pairs of individuals in inbred families, "relation" = "inbred family". }

}


\author{Cecelia Laurie

}


\references{
 Ani Manichaikul, Josyf C. Mychaleckyj, Stephen S. Rich, Kathy Daly, Michele Sale, and Wei-Min Chen.
Robust relationship inference in genome-wide association studies. Bioinformatics, 26(22):2867-2873,
2010.
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\seealso{ \code{\link{pedigreeCheck}}, \code{\link{pedigreePairwiseRelatedness}}, \code{\link{pedigreeMaxUnrelated}}

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
family <- c(1,1,1,1,2,2,2,2,2,2,2)
individ <- c(1,2,3,4,5,6,7,8,9,10,11)
mother <- c(0,0,1,1,0,0,5,5,0,0,10)
father <- c(0,0,2,2,0,0,6,9,0,0,7)
sex <- c("F","M","F","F","F","M","M","M","M","F","F")
pedigree <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
expRelsCategory(pedigree)

# inbred family
family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp2 <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
expRelsCategory(samp2)

# deeper pedigree
ped<-data.frame("family"=rep(1,13),"individ"=1:13,"father"=rep(0,13),"mother"=rep(0,13),"sex"=rep(NA,13),stringsAsFactors=FALSE)
ped$sex<-c("F","M","F","F","M","F","M","F","M","F","F","M","F")
ped$father<-c(0,0,2,2,0,5,5,0,0,7,9,0,12)
ped$mother<-c(0,0,1,1,0,4,4,0,0,8,6,0,11)
expRelsCategory(ped)
}

\keyword{manip }


