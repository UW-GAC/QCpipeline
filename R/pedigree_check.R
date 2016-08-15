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
              "out_plot"="ibd.pdf",
              "scan_include_file"=NA)
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
    select(pair, family, relation, exp.rel, MZtwinID)

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

table(ibd$exp.rel, ibd$obs.rel)

## plot unexpected relatives
ibd <- mutate(ibd, unexp=ifelse(exp.rel == obs.rel, "expected", "unexpected"))

rels <- c("Dup", "PO", "FS", "Deg1", "Deg2", "Deg3", "Q", "U")
cols <- c(brewer.pal(length(rels)-1, "Dark2")[c(1, 2, 3, 6, 5, 4, 7)], "black")
cmap <- setNames(cols, rels)

theme_set(theme_bw() + theme(legend.position=c(1, 1), legend.justification=c(1,1), legend.background = element_rect(colour = "black")))

p <- ggplot(ibd, aes(k0, kin, color=exp.rel)) + facet_wrap(~unexp) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype='dashed', color="grey") +
    geom_point() +
    scale_color_manual(values=cmap, breaks=names(cmap)) +
    ylab("kinship estimate")

ggsave(config["out_plot"], plot=p, width=12, height=6)

