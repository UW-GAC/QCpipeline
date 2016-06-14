## Convert AnnotatedDataFrame with sample.id column to ScanAnnotationDataFrame

ADFtoScanADF <- function(adf) {
    stopifnot(is(adf, "AnnotatedDataFrame"))
    stopifnot("sample.id" %in% varLabels(adf))
    varLabels(adf)[varLabels(adf) == "sample.id"] <- "scanID"
    ScanAnnotationDataFrame(pData(adf), varMetadata(adf))
}
