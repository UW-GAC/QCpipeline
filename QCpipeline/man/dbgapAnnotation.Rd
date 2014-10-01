\name{dbgapAnnotation}
\alias{dbgapScanAnnotation}
\alias{dbgapSnpAnnotation}
\title{Write annotation files for dbGaP}
\description{\code{dbgapScanAnnotation} and \code{dbgapSnpAnnotation}
  create text files appropriate for posting on dbGaP.}
\usage{
dbgapScanAnnotation(scanAnnot, dir=".",
   consentVar=NULL, subjVar="subj.plink", dupVar="dup.post", omitVar="no.post",
   annotationCol="annotation", analysisCol="analysis")

dbgapSnpAnnotation(snpAnnot, dir=".",
   annotationCol="annotation", analysisCol="analysis")
}
\arguments{
  \item{scanAnnot}{A \code{\link{ScanAnnotationDataFrame}}.}
  \item{snpAnnot}{A \code{\link{SnpAnnotationDataFrame}}.}
  \item{dir}{A character string with the directory for file output.}
  \item{consentVar}{The variable in \code{scanAnnot} containing consent
    levels. If \code{NULL}, the output files will not be split by consent level.}
  \item{subjVar}{The logical variable in \code{scanAnnot} indicating unique
    subjects to post.}
  \item{dupVar}{The logical variable in \code{scanAnnot} indicating
    duplicate scans to post.}
  \item{omitVar}{The logical variable in \code{scanAnnot} indicating
    scans to be omitted from posting.}
  \item{annotationCol}{The logical column in the metadata indicating
    which variables should be included in the annotation files.}
  \item{analysisCol}{The logical column in the metadata indicating
    which variables should be included in the analysis files.}
}
\details{
  \code{dbgapScanAnnotation} writes the following files to \code{dir}:
  \itemize{
    \item{Sample_annotation.csv}
    \item{Sample_annotation_consent_*.csv}
    \item{Sample_annotation_duplicates.csv}
    \item{Sample_annotation_duplicates_consent_*.csv}
    \item{Sample_annotation_DD.txt}
    \item{Sample_analysis.csv}
    \item{Sample_analysis_duplicates.csv}
    \item{Sample_analysis_DD.txt}
  }

  \code{dbgapSnpAnnotation} writes the following files to \code{dir}:
  \itemize{
    \item{SNP_annotation.csv}
    \item{SNP_annotation_DD.txt}
    \item{SNP_analysis.csv}
    \item{SNP_analysis_DD.txt}
  }

  Which variables should be written to the annotation and analysis files
  are indicated in the metadata columns \code{annotationCol} and
  \code{analysisCol}.

  The data dictionary files are populated from the metadata.  The "type"
  column is automatically generated from the classes of the variables in
  \code{scanAnnot} and \code{snpAnnot}.
}
\author{Stephanie Gogarten}
%\seealso{}
%\examples{}
