plateLayoutPlot <- function(samples,out.file,nprob){
  # samples must be a data.frame with columns labeled as "scanID","plate","well","annot.sex","sex","obsdup.id","unexp.dup","unobsdup.id"
  # where "annot.sex" is annotated (expected) sex and "sex" is observed sex; both must have values of "F" (female), "M" (male) or "U" (unknown)
  # "obsdup.id" is "integer identifier for observed duplicate set (i.e. all members of a set are dups of one another)"
  # "unobsdup.id" is "integer identifier for expected duplicate set in which the members were not observed to be duplicates"
  # "unexp.dup" is "logical indicator for dup sets that were not expected to be duplicates"
  # "nprob" is the minimum number of problematic samples per plate required for inclusion in the plots
  
  # check rows of data frame
  if(!class(samples)=="data.frame") stop("samples must be a data.frame")
  if(!(all(is.element(c("scanID","plate","well","annot.sex","sex","obsdup.id","unexp.dup","unobsdup.id"),names(samples))))) stop("samples must have column names scanID, plate, well, annot.sex, sex, obsdup.id, unexp.dup, unobsdup.id, [contaminated, hilite]")
  
  # defaults for optional rows
  if (!("hilite" %in% names(samples))) samples$hilite <- FALSE
  if (!("contaminated" %in% names(samples))) samples$contaminated <- FALSE
  
  # input values
  if(!(all(is.element(samples$sex,c("M","F","U")) &
           is.element(samples$annot.sex,c("F","M","U"))))) stop("sex and annot.sex values must be M, F or U")
  if(!is.logical(samples$unexp.dup)) stop("unexp.dup must be logical")
  if(!is.integer(samples$obsdup.id)) stop("obsdup.id must be integer values")
  if (!is.logical(samples$contaminated)) stop("contaminated must be logical")
  if (!is.logical(samples$hilite)) stop("hilite must be logical")
  
  nums <- c(paste("0", 1:9, sep=""), 10:12)
  wells <- NULL
  for(i in nums) wells <- c(wells, paste(LETTERS[1:8], i, sep=""))
  # check wells
  if(!(all(is.element(samples$well, wells)))) stop("wells must be A-H combined with 1-12")
  
  # parse out sample plate row and col
  samples$samp.col <- as.numeric(substr(samples$well, 2, 3))
  samples$samp.row <- substr(samples$well, 1, 1)
  # make numeric row (to represent A-H
  tmp <- data.frame(samp.row=sort(unique(samples$samp.row)), samp.row.num=8:1, stringsAsFactors=F)  # reverse so that A=8 and H=1
  samples <- merge(samples, tmp, all.x=T); dim(samples)
  # make data for plotting positions of all wells
  wells <- samples[!duplicated(samples$well), c("well", "samp.row", "samp.col", "samp.row.num")]
  
  # identify plates containing problematic samples
  if (!("id.prob" %in% names(samples))){
    samples$id.prob <- (samples$annot.sex!="U" & samples$annot.sex!=samples$sex) | !is.na(samples$unobsdup.id) | (!is.na(samples$unexp.dup) & samples$unexp.dup) | samples$contaminated | samples$hilite
  }
  
  # determine plates with problems
  pl <- table(samples$plate, samples$id.prob)[,"TRUE"]
  pl <- names(pl)[pl >= nprob]
  n <- length(pl)
  cex.well <- 4.5 # cex for wells on plate
  
  # plot plate layout with positions of problematic samples
  pdf(out.file)
  for(i in 1:n){
    plate <- pl[i]
    
    mtxt1 <- "fill color=genetic sex (pink=F,blue=M,yellow=U,white=empty), red circle=sex mismatch\nnumbers: black in well=exp & obs dup, red in well=unexp & obs dup, above well=exp & unobs dup\ncrosses: red=contaminated sample, black=hilite sample"
    
    chk <- samples[is.element(samples$plate, plate), ]
    
    plot(wells$samp.col, wells$samp.row.num, xlim=c(0.5,12.5), ylim=c(0.5,8.5), col="gray",
         cex=cex.well, xlab=mtxt1, ylab="", main=plate, yaxt="n", xaxt="n",
         cex.sub=0.7, cex.lab=0.7)
    # row and column numbers/letters
    axis(side=2, at=1:8, labels=c("H","G","F","E","D","C","B","A"), las=2)
    axis(side=1, at=1:12, labels=1:12)
    
    ## gender mismatch
    sel <- chk$annot.sex != chk$sex
    points(chk$samp.col[sel], chk$samp.row.num[sel], col="red", cex=cex.well+0.5, pch=16)
    # fill according to annotated sex
    sex.col <- ifelse(chk$sex == "M", "lightblue",
                      ifelse(chk$sex == "F", "pink",
                             ifelse(chk$sex == "U", "yellow", NA)))
    points(chk$samp.col, chk$samp.row.num, col=sex.col, cex=cex.well, pch=16)
    
    ## contaminated samples
    sel <- chk$contaminated
    points(chk$samp.col[sel], chk$samp.row.num[sel], cex=cex.well, col="red", pch=4, lwd=1)
    
    ## hilite samples
    sel <- chk$hilite
    points(chk$samp.col[sel], chk$samp.row.num[sel], cex=cex.well, col="black", pch=4, lwd=1)
    
    
    ## expected and observed dups
    sel <- !is.na(chk$obsdup.id) & !chk$unexp.dup
    lbl <- chk$obsdup.id[sel]
    if(length(lbl) > 0) text(chk$samp.col[sel], chk$samp.row.num[sel], col="black", labels=lbl, cex=1)
    # label unexpected dups
    sel <- !is.na(chk$obsdup.id) & chk$unexp.dup
    lbl <- chk$obsdup.id[sel]
    if(length(lbl) > 0) text(chk$samp.col[sel], chk$samp.row.num[sel], col="red", labels=lbl, cex=1)	
    # label expected but unobserved
    sel <- !is.na(chk$unobsdup.id) 
    lbl <- chk$unobsdup.id[sel]
    if(length(lbl) > 0) text(chk$samp.col[sel], chk$samp.row.num[sel]+0.5, col="purple", labels=lbl, cex=1)	
    
    
  }
  dev.off()
  print(paste(n, "plates displayed"))
  data.frame(plate=pl, page=1:length(pl), stringsAsFactors=F)
}
