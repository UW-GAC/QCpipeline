##########
# Plate layout plots
# Usage: R --args config.file < plate_layout.R
##########

library(GWASTools)
library(QCpipeline)
sessionInfo()

# read configuration
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("missing configuration file")
config <- readConfig(args[1])

# check config and set defaults
required <- c("annot_scan_file", "ibd_obsrel_file")
optional <- c("annot_scan_subjectCol",
              "annot_scan_sexCol", "annot_scan_annotSexCol",
              "annot_scan_annotSexMale", "annot_scan_annotSexFemale", "annot_scan_annotSexUnknown", 
              "annot_scan_plateCol", "annot_scan_wellCol",
              "ibd_unobs_dup_file",
              "min_num_problems",
              "out_annot_file",
              "out_plate_plot",
              "out_plate_file",
              "scan_hilite_file",
              "scan_contaminated_file")
default <- c("subjectID",
  "sex", "annot.sex",
             "M", "F", NA,
             "Sample.Plate", "Sample.Well",
             NA,
             1,
             "out_scan_annot.RData",
             "plate_layout.pdf",
             "plate_layout.RData",
             NA,
             NA)
config <- setConfigDefaults(config, required, optional, default)
print(config)

(scanAnnot <- getobj(config["annot_scan_file"]))
scanID <- getScanID(scanAnnot)


## sex mismatches
# make sure annotated sex column is present, if requested
# annot.sex is not always kept if there are no sex mismatches.
stopifnot(hasVariable(scanAnnot, config["annot_scan_annotSexCol"]))
sex <- getVariable(scanAnnot, config["annot_scan_sexCol"])
annot.sex <- getVariable(scanAnnot, config["annot_scan_annotSexCol"])
# set annot.sex to match
annot.sex[annot.sex %in% config["annot_scan_annotSexMale"]] <- "M"
annot.sex[annot.sex %in% config["annot_scan_annotSexFemale"]] <- "F"
annot.sex[annot.sex %in% config["annot_scan_annotSexUnknown"]] <- NA
i_sex_mismatch <- !is.na(annot.sex) & !is.na(sex) & (sex != annot.sex)
table(i_sex_mismatch, exclude=NULL)



## unexpected but observed duplicate pairs - from IBD results
obsrels <- getobj(config["ibd_obsrel_file"])
dim(obsrels)
stopifnot("exp.rel" %in% names(obsrels))

# select out the unexpected duplicates and give a number to each group
dups <- obsrels[obsrels$obs.rel == "Dup", ]
names(dups)[names(dups) == "ID1"] <- "sample1"
names(dups)[names(dups) == "ID2"] <- "sample2"
names(dups)[names(dups) == "kinship"] <- "KC"
if (nrow(dups) > 1){
  dup.sets <- defineFamilies(dups, KC.threshold=0.45)
  names(dup.sets)[names(dup.sets) == "family"] <- "obsdup.id"
} else {
  # there are no dups
  dup.sets <- dups
  dup.sets$obsdup.id <- numeric()
}  
dup.sets$unexp.dup <- dup.sets$exp.rel != "Dup"
unexp.dup.sets <- dup.sets[dup.sets$exp.rel != "Dup", ]
i_unexp_dups <- scanID %in% c(unexp.dup.sets$sample1, unexp.dup.sets$sample2)
table(i_unexp_dups, exclude=NULL)

## add in check that each scan in in only one dup.set (family)


## expected but unobserved duplicate pairs
if (!is.na(config["ibd_unobs_dup_file"]) & file.exists(config["ibd_unobs_dup_file"])){
  unobs.dups <- getobj(config["ibd_unobs_dup_file"])
  i_unobs_dups <- scanID %in% unlist(unobs.dups)
} else {
  message("No unexpected duplicate file given or file not found.")
  unobs.dups <- NULL
  unobs.dup.id <- NULL
  i_unobs_dups <- rep(FALSE, nrow(scanAnnot))
}
table(i_unobs_dups, exclude=NULL)

## contaminated samples
if (!is.na(config["scan_contaminated_file"]) & file.exists(config["scan_contaminated_file"])){
  contam <- getobj(config["scan_contaminated_file"])
  i_contam <- scanID %in% contam
} else {
  message("No contaminated samples file given or file not found.")
  i_contam <- rep(FALSE, length(scanID))
}
table(i_contam, exclude=NULL)

## user-hilighted samles
if (!is.na(config["scan_hilite_file"]) & file.exists(config["scan_hilite_file"])){
  hilite <- getobj(config["scan_hilite_file"])
  i_hilite <- scanID %in% hilite
} else {
  message("No user-highlighted samples file given or file not found.")
  i_hilite <- rep(FALSE, length(scanID))
}
table(i_hilite, exclude=NULL)

# samples with any problems
i_problems <- i_sex_mismatch | i_unexp_dups | i_unobs_dups | i_contam | i_hilite
table(i_problems, exclude=NULL)


samp <- getVariable(scanAnnot,
                           varname=c("scanID", config["annot_scan_subjectCol"],
                                     config["annot_scan_plateCol"],
                                     config["annot_scan_wellCol"],
                                     config["annot_scan_annotSexCol"],
                                     config["annot_scan_sexCol"]))
stopifnot(allequal(samp$scanID, sort(samp$scanID)))
names(samp) <- c("scanID", "subjectID", "plate", "well", "annot.sex", "sex")

# add id prob
samp$id.prob <- i_problems

## add in obsdup.id and unexp.dup
x1 <- match(samp$scanID, dup.sets$sample1)
x2 <- match(samp$scanID, dup.sets$sample2)
#stopifnot(!any(!is.na(x1) & !is.na(x2)))
x <- x1
x[is.na(x1)] <- x2[is.na(x1)]
samp$obsdup.id <- as.integer(dup.sets$obsdup.id[x])
samp$unexp.dup <- dup.sets$unexp.dup[x]
table(samp$unexp.dup, exclude=NULL)


## add in unobsdup.id
samp$unobsdup.id <- NA
if (!is.null(unobs.dups)){
  for (i in 1:length(unobs.dups)){
    samp$unobsdup.id[samp$scanID %in% unobs.dups[[i]]] <- i
  }
}
table(samp$unobsdup.id, exclude=NULL)

## add in contaminated
samp$contaminated <- i_contam

## add in user highlighted
samp$hilite <- i_hilite

# annot.sex was updated earlier -- just change NA to U
samp$annot.sex[samp$annot.sex %in% config["annot_scan_annotSexMale"]] <- "M"
samp$annot.sex[samp$annot.sex %in% config["annot_scan_annotSexFemale"]] <- "F"
samp$annot.sex[samp$annot.sex %in% config["annot_scan_annotSexUnknown"]] <- "U"
table(samp$annot.sex, exclude=NULL)

# id problems per plate
table(samp$plate, samp$id.prob)

if (sum(samp$id.prob) > 1){
  # now run the plate map code
  outf <- config["out_plate_plot"]
  plate.file <- plateLayoutPlot(samp, outf, nprob=config["min_num_problems"])
  save(plate.file, file=config["out_plate_file"])
} else {
  message("No sample identity issues found.")
}

sampAnnot <- ScanAnnotationDataFrame(samp)
save(sampAnnot, file=config["out_annot_file"])

