\name{AnnotationByChr}
\docType{class}

\alias{AnnotationByChr-class}
\alias{AnnotationByChr}

\alias{getSnpAnnotation,AnnotationByChr-method}

\alias{getValidChromosomes,AnnotationByChr-method}

\alias{lookUpSnps}
\alias{lookUpSnps,AnnotationByChr,missing-method}

\title{Class AnnotationByChr}

\description{
  The AnnotationByChr class is a wrapper to work with SNP annotation files for stored by chromosome.
  The class relies on the underlying directory structure to be set up in specific way.
  All data sets must have one SNP annotation file per chromosome in one directory, with a naming convention BASE_chr1.RData, BASE_chr2.RData (where BASE can be specified for each directory).
}


\section{Constructor}{
  \describe{
    \item{}{
      \code{AnnotationByChr(directory, base=NULL)}:
      
      \code{directory} is the path to the directory.
      \code{base} is (optionally) the base name. The code will automatically attempt to detect base from the files, assuming the naming conventions are met.
        

      The \code{AnnotationByChr} constructor creates and returns a
      AnnotationByChr instance that works with the annotation files contained in the directory.
      
    }
  }
}

\section{Accessors}{
  
  In the code snippets below, \code{object} is an AnnotationByChr object. The SnpAnnotationDataFrame files are assumed to be named \code{object@base_chr-<chromosome>.RData}.

  \describe{
    \item{}{
      \code{getSnpAnnotation(object, chromosome)}: Returns the SnpAnnotationDataFrame associated with \code{object} for a given chromosome(s).
    }

    \item{}{
      \code{getValidChromosomes(object)}: Returns the set of chromosomes with RData files in this directory.
    }

    \item{}{
      \code{lookUpSnps(object, snps=NULL, column="snpID", chromosome=getValidChromosomes(object), extraCols=NULL, ...)}: Returns a SnpAnnotationDataFrame associated with \code{object}. If \code{snps=NULL}, all rows are returned from the SNP annotations for chromosomes specified with \code{chromosome}. If \code{snps=NULL}, only rows where \code{snps} matches the entries in \code{column} are returned. Specify \code{chromosome} to restrict the search to a subset of chromosomes.  \code{snps}, \code{column} and \code{chromosome} must be named arguments (a requirement of method dispatch). \code{...} are passed to \code{\link{getSnpAnnotation}}.
    }
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{GenotypeDataByChr}}
}

\examples{
  \dontrun{
    annot <- AnnotationByChr("annot_directory")  
    snps <- lookUpSnps(annot, snps=c("rs101", "rs202", "rs303"), column="rsID")
  }
}

\keyword{methods}
\keyword{classes}

