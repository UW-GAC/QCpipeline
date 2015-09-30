# R function to run command line BLAT tool on a series of FASTA DNA sequences
# S Nelson, UWGCC, 8/9/2013
# migrated to QCannot package 9/21/2015

runBlat <- function(seq.input, make.fasta=TRUE,
                    resource.path="/projects/geneva/gcc-fs2/SNP_annotation/UCSC_downloads/blat_resources/",
                     output.fn="blat_query",  ooc=TRUE,
                     tileSize=11, stepSize=tileSize,
                     repMatch=1024, minScore=30, minIdentity=90,
                     find.best=TRUE, out.type="psl", return.result=FALSE)
{
  
  # Currently script only works when running R on non-Windows platform
  sys.windows <-  Sys.info()["sysname"]=="Windows"
  if(sys.windows) stop("I've detected you are on a Windows system -- due to dependency on command-line BLAT tool, this function can only be run on non-Windows platforms")

  # requires Biostrings and bitops. importing both with QCannot 
  
  ## check function arguments
  # check logicals
  stopifnot(is.logical(make.fasta))
  stopifnot(is.logical(ooc))
  stopifnot(is.logical(find.best))
  stopifnot(is.logical(return.result))

  # check BLAT program specs
  stopifnot(is.numeric(tileSize))
  if(tileSize < 8 | tileSize > 12) warning("tileSize is usually between 8 and 12")
  
  stopifnot(is.numeric(stepSize))
  if(stepSize < 8 | stepSize > 12) warning("stepSize is usually between 8 and 12")

  stopifnot(is.numeric(repMatch))
  if(repMatch < 256 | repMatch > 4096) warning("Typically repMatch is 256 for tileSize 12, 1024 for tile size 11, 4096 for tile size 10. Default is 1024")

  stopifnot(is.numeric(minScore))

  stopifnot(is.numeric(minIdentity))
  if(minIdentity < 0 | minIdentity > 100) warning("minIdentity should range between 0 and 100; 90 is the default for nucleotide searches")

  # if pointing to pre-exisiting FASTA file, check that file exists
  if(!make.fasta){
    if(file.exists(seq.input)) stop("Specified input FASTA file does not exist: ", seq.input)  }
  
   # if asking to find best hit(s), check that output type is psl
  if (find.best & out.type!="psl")
   {warning("Cannot apply 'near best in genome filter' (find.best=TRUE) unless BLAT output type is psl (out.type='psl')")}

   # construct file paths to BLAT resource files and check they exist
  db.path <- file.path(resource.path, "hg19.2bit")
  if(!file.exists(db.path)) stop("resource file ", db.path," cannot be found")
  ooc.path <- file.path(resource.path, "11.ooc")
  if(ooc & !file.exists(ooc.path)) stop("resource file ", ooc.path," cannot be found")
  
  ### end checks

  # check whether input needs to be converted to FASTA format
  if (make.fasta) {
    message("Converting user input into FASTA formatted file\n")
     
  # convert to FASTA sequence
    seq <- seq.input[,2]
    names(seq) <- seq.input[,1]
    #dna=DNAStringSet(seq)
    ## replacing with more generic "BStringSet" to allow for "[A/G]" type SNP notation
    dna <- BStringSet(seq)
    fasta.fn <- paste(output.fn,".fa",sep="")
    message("Writing FASTA formatted input file to ",fasta.fn,"\n")
    writeXStringSet(dna, fasta.fn)
  }

  # if "seq.input" input is already a FASTA file saved on the filesystem
  if (!make.fasta){
    message("User input bexing interpreted as file path to pre-existing FASTA file\n")
    fasta.fn <- seq.input
  }

  ### run BLAT query tool
  
  # create system command to submit the BLAT job
  ooc.arg <- ""
  if(ooc) ooc.arg <- paste0(" -ooc=", ooc.path)
  
  filo <- paste(output.fn, out.type, sep=".")
  
  cmnd <- paste0("blat ", ooc.arg, " -minScore=", minScore, " -minIdentity=", minIdentity, " -tileSize=", tileSize, " -stepSize=", stepSize, " -repMatch=", repMatch, " -out=", out.type, " ", db.path, " ", fasta.fn, " ", filo)
  
  message(paste("Submitting following command line BLAT query:\n", cmnd, "\n"))
  system(cmnd)

  ### if also finding "near best in genome" results, run pslReps program
  if(find.best & out.type=="psl"){

    # build second command line argument
    message("Refining query results to get 'near best in genome' matches\n")

    topalign <- paste0(output.fn, "_topalignments.psl")
    topalign_psr <- paste0(output.fn, "_topalignments.psr")

    cmnd2 <- paste0("pslReps -minCover=0.15 -minAli=0.96 -nearTop=0.001 ", filo, " ", topalign, " ", topalign_psr)

    message(paste("Submitting following command line to run pslReps:\n", cmnd2, "\n"))
    system(cmnd2)

    # delete psr file
    file.remove(topalign_psr)
  }

  ### process results file (if psl)
  if (return.result & out.type=="psl") {
    rslt <- read.table(file=paste(output.fn,out.type,sep="."),skip=5)
    psl.names <-  c("match","mismatch","repmatch","Ns","Q_gap_count","Q_gap_bases","T_gap_count","T_gap_bases", "strand","Q_name","Q_size","Q_start","Q_end","T_name","T_size","T_start","T_end","block_count", "blocksizes","qStarts","tStarts") 
    names(rslt) <- psl.names

    #### calculate score
    ## see http://genome.ucsc.edu/FAQ/FAQblat.html#blat4 (C code)
    cat("\n\tCalculating score...\n")
    psl <- rslt
    psl$score <- psl$match + bitShiftR(psl$repmatch, 1) - psl$mismatch - psl$Q_gap_count - psl$T_gap_count

    #### calcluate a percent identity
    ## see https://groups.google.com/a/soe.ucsc.edu/forum/#!searchin/genome/BLAT$20psl$20qnuminsert/genome/-N6XU5Qh_Lo/f6123GfxsRIJ, title "BLAT vs Browser off by one base"
    message("\tCalculating percent identity")
# , per formula from Jim Kent posted to UCSC Genome Browser discussion list (note may differ from Web result)...\n")
    psl$pct.identity <- round(100*((psl$match + psl$repmatch)/
                                   (psl$mismatch + psl$match + psl$repmatch + psl$Q_gap_count + psl$T_gap_count)),1)

      rslt <- psl

    # if found best alignment, read in top results file and add flag
    if (find.best) {
      top.rslt <-  read.table(file=topalign, skip=5)
      names(top.rslt) <- psl.names

      rslt$id <- 1:nrow(rslt)
      
      # merge top hits with original list to identify
      comb <- merge(rslt, top.rslt, by=c("T_name","T_start","T_end","block_count"),
                    all.x=FALSE, all.y=TRUE, sort=FALSE)
      rslt$top.hit <- is.element(rslt$id, comb$id)
    } # close if loop in reading in top alignment

    return(rslt)
  } # close if loop on psl output type

  # if not psl, direct user to results file on Unix system
  if (!return.result | out.type!="psl") {
    cat("see BLAT",out.type,"output at:",paste(output.fn,out.type,sep="."), "\n")
  }

  # if not psl and user requested results returned as R object, echo warning
  if (return.result | out.type!="psl") {
    warning("will only return R object of BLAT results if using .psl output type\n")
  }    

} # close function
