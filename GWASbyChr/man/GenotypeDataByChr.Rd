\name{GenotypeDataByChr}
\docType{class}

\alias{GenotypeDataByChr-class}
\alias{GenotypeDataByChr}

\alias{getGenoData}
\alias{getGenoData,GenotypeDataByChr-method}

\alias{getScanAnnotation,GenotypeDataByChr-method}

\alias{getSnpAnnotation,GenotypeDataByChr-method}

\alias{getValidChromosomes}
\alias{getValidChromosomes,GenotypeDataByChr-method}

\alias{getScanID,GenotypeDataByChr-method}

\alias{hasSnpID}
\alias{hasSnpID,GenotypeDataByChr-method}

\alias{getGenotypeFromSnpID}
\alias{getGenotypeFromSnpID,GenotypeDataByChr-method}

\title{Class GenotypeDataByChr}

\description{
  The GenotypeDataByChr class is a wrapper to work with imputed data gds files stored by chromosome.
  The class relies on the underlying directory structure to be set up in specific way.
  All imputed data sets must have one GDS file per chromosome in one directory, with a naming convention BASE_chr1.gds, BASE_chr2.gds (where BASE can be specified for each directory).
  If present, the SnpAnnotationDataFrame and ScanAnnotationDataFrame objects in the directory must be named as \code{BASE_chr1_snpAnnot.RData}, \code{BASE_chr2_scanAnnot.RData}, etc.
}


\section{Constructor}{
  \describe{
    \item{}{
      \code{GenotypeDataByChr(directory, base=NULL)}:
      
      \code{directory} is the path to the directory.
      \code{base} is (optionally) the base name. The code will automatically attempt to detect base from the files, assuming the naming conventions are met.
        

      The \code{GenotypeDataByChr} constructor creates and returns a
      GenotypeDataByChr instance that works with the imputed dosage files contained in the directory.
      
    }
  }
}

\section{Accessors}{
  
  In the code snippets below, \code{object} is an GenotypeDataByChr object. The specified directory should contain all GDS files, all SnpAnnotationDataFrames, and a ScanAnnotationDataFrame that matches all GDS files. The GDS files containing the genotypes are assumed to be named \code{object@base_chr-<chromosome>.gds}. The SnpAnnotationDataFrame associated with each chromosome GDS file are assumed to be named \code{object@base_chr-<chromosome>_snpAnnot.RData}. The ScanAnnotationDataFrame associated with the GDS file is assumed to be named \code{object@base_scanAnnot.RData}.

  \describe{
    \item{}{
      \code{getGenoData(object, chromosome, snpAnnot=NULL, scanAnnot=NULL)}: Returns a GenotypeData object pointing to the specified chromosome, with optional SNP and scan annotation objects attached. If snpAnnot is a SnpAnnotationDataFrame or scanAnnot is a ScanAnnotationDataFrame, that annotation is passed in the call to GenotypeData. If snpAnnot or scanAnnot is \code{TRUE}, the SnpAnnotationDataFrame or ScanAnnotationDataFrame is read from the directory (assuming they are named \code{object@base_chr-<chromosome>_snpAnnot.RData} or \code{object@base_scanAnnot.RData} and attached.
    }

    \item{}{
      \code{getScanAnnotation(object)}: Returns the scan annotation associated with \code{object}.
    }
    
    \item{}{
      \code{getSnpAnnotation(object, chromosome)}: Returns the SnpAnnotationDataFrame associated with \code{object} for a given chromosome.
    }

    \item{}{
      \code{getValidChromosomes(object)}: Returns the set of chromosomes with GDS files in this directory.
    }

    \item{}{
      \code{getScanID(object)}: Returns the scanIDs for these GDS files.
    }


    \item{}{
      \code{hasSnpID(object, snpID)}: Returns whether \code{snpID} is present in any of the GDS files in this directory.
    }


    \item{}{
      \code{getGenotypeFromSnpID(object, snpID, drop=FALSE, order=c("selection", "file"), transpose=FALSE)}: Returns genotypes for the snps in \code{snpID}.
    }

  }
}

\author{Adrienne Stilp}

\seealso{
  \code{\link{GenotypeData}}
}

\examples{
  \dontrun{
    dataset <- GenotypeDataByChr("gds_directory")
  
    genoData_chr1 <- getGenoData(dataset, 1)
    close(genoData_chr1)
  }
}

\keyword{methods}
\keyword{classes}

