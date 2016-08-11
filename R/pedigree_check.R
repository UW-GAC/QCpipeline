library(GWASTools)
library(QCpipeline)
library(gdsfmt)
library(GENESIS)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
sessionInfo()


## read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

## check config and set defaults
required <- c("annot_scan_file",
              "exp_rel_file",
              "ibd_file")
optional <- c("annot_scan_subjectCol"="subjectID",
              "ibd_method"="pcrelate",
              "out_ibd_file"="ibd_obsrel.RData",
              "out_twins_file"="ibd_twins.RData",
              "scan_include_file"=NA,
              "unexpected_threshold"="deg2")
config <- setConfigDefaults(config, required, optional=names(optional), default=optional)
print(config)

annot <- getobj(config["annot_scan_file"])
## to-do: allow ScanAnnotationDataFrame also
if (is(annot, "ScanAnnotationDataFrame")) stop("annot_scan_file must be an AnnotatedDataFrame with sample.id")
annot <- pData(annot) %>%
    select_("sample.id", config["annot_scan_subjectCol"]) %>%
    rename_(Individ=config["annot_scan_subjectCol"])

if (!is.na(config["scan_include_file"])) {
    ids <- getobj(config["scan_include_file"])
    annot <- filter(annot, sample.id %in% ids)
} else {
    ids <- NULL
}

## read ibd file
## currently only PC-Relate format, maybe allow others also?
stopifnot(config["ibd_method"] == "pcrelate")
pcr <- openfn.gds(config["ibd_file"])
ibd <- pcrelateReadKinship(pcr, scan.include=ids)
closefn.gds(pcr)
nrow(ibd)

## expected relatives
rel <- getobj(config["exp_rel_file"])
rel <- rel %>%
    mutate(pair=pasteSorted(Individ1, Individ2)) %>%
    select(pair, family, relation, exp.rel, MZtwinID, twin.fam)

ibd <- select(ibd, ID1, ID2, kin, k0) %>%
    left_join(annot, by=c(ID1="sample.id")) %>%
    rename(Individ1=Individ) %>%
    left_join(annot, by=c(ID2="sample.id")) %>%
    rename(Individ2=Individ) %>%
    mutate(pair=pasteSorted(Individ1, Individ2)) %>%
    left_join(rel, by="pair") %>%
    select(-pair) %>%
    mutate(exp.rel=ifelse(is.na(exp.rel), "U", exp.rel))

## observed relatives
ibd <- ibd %>%
    mutate(obs.rel=ibdAssignRelatednessPCRelate(k0, kin)) %>%
    filter(!(exp.rel == "U" & obs.rel == "U"))
nrow(ibd)
save(ibd, file=config["out_ibd_file"])

## separate families with MZ twins
twin.fam.subj <- filter(ibd, twin.fam) %>%
    select(ID1, ID2) %>%
    unlist(use.names=FALSE) %>%
    unique()
ibd.twins <- filter(ibd, ID1 %in% twin.fam.subj | ID2 %in% twin.fam.subj)
save(ibd.twins, file=config["out_twins_file"])

ibd <- filter(ibd, !twin.fam)
table(ibd$exp.rel, ibd$obs.rel)
