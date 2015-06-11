\name{AssocResultsByChr}
\docType{class}

\alias{AssocResultsByChr-class}
\alias{AssocResultsByChr}

\alias{getAssocResults}
\alias{getAssocResults,AssocResultsByChr-method}

\alias{getValidChromosomes,AssocResultsByChr-method}

\alias{lookUpSnps,AssocResultsByChr,missing-method}
\alias{lookUpSnps,AnnotationByChr,AssocResultsByChr-method}

\title{Class AssocResultsByChr}

\description{
  The AssocResultsByChr class is a wrapper to work with association test results files stored by chromosome.
  The class relies on the underlying directory structure to be set up in specific way.
  All data sets must have one association results file per chromosome in one directory, with a naming convention BASE_chr1.RData, BASE_chr2.RData (where BASE can be specified for each directory).
}


\section{Constructor}{
  \describe{
    \item{}{
      \code{AssocResultsByChr(directory, base, chromSep, suffix)}:
      
      \code{directory} is the path to the directory.
      \code{base} is (optionally) the base name. The code will automatically attempt to detect base from the files, assuming the naming conventions are met.
      \code{chromSep} separates the base from the chromosome string. Default is "_chr".
      \code{suffix} comes after the chromosome string. Default is ".RData".

      The \code{AssocResultsByChr} constructor creates and returns a
      AssocResultsByChr instance that works with the files contained in the directory.
      
    }
  }
}

\section{Accessors}{
  
  In the code snippets below, \code{object} is an AssocResultsByChr object. The association results files are assumed to be named \code{object@base_chr-<chromosome>.RData}.

  \describe{
    \item{}{
      \code{getAssocResults(object, chromosome, returnColumns=NULL)}: Returns the data frame associated with \code{object} for a given chromosome.  If \code{chromosome=NULL} (default), returns all results in a single data frame. The \code{returnColumns} argument can be used to select columns to return; any columns that appear in the association data frames that intersect with \code{returnColumns} are returned. If \code{returnColumns=NULL}, all columns are returned by default.
    }

    \item{}{
      \code{getValidChromosomes(object)}: Returns the set of chromosomes with RData files in this directory.
    }

    \item{}{
      \code{lookUpSnps(object, snps=NULL, column="snpID", chromosome=getValidChromosomes(object1), returnColumns=NULL, ...)}: Returns a data frame of results associated with \code{object}. If \code{snps=NULL}, results are concatenated into a single data frame from all chromosomes specified by \code{chromosome}. If \code{snps} is not \code{NULL}, then only rows where the entries in \code{column} match \code{snps} are returned. \code{snps}, \code{column} and \code{chromosome} must be named arguments (a requirement of method dispatch).  \code{...} are passed to the \code{getSnpAnnotation} method for \code{AnnotationByChr}.
    }

    \item{}{
      \code{lookUpSnps(object1, object2, snps, column="snpID", chromosome=getValidChromosomes(object2), returnColumns=NULL, ...)}: \code{object1} is a \code{\link{AnnotationByChr}} object, and \code{object2} is a \code{AssocResultsByChr} object.  
      Returns a single data frame with columns from both \code{lookUpSnps} for \code{AnnotationByChr} and \code{AssocResultsByChr}.  \code{...} are passed to the \code{getSnpAnnotation} method for \code{AnnotationByChr}.
    }
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{AnnotationByChr}}
}

\examples{
  \dontrun{
    assocSet <- AssocResultsByChr("assoc_directory")  
    assoc <- getAssocResults(assocSet, chromosome=1)
    snps <- lookUpSnps(assocSet, snps=c(101,202,303), column="snpID")

    annotSet <- AnnotationByChr("annot_directory")  
    snps <- lookUpSnps(annotSet, assocSet, snps=c("rs101", "rs202", "rs303"), column="rsID")
  }
}

\keyword{methods}
\keyword{classes}

