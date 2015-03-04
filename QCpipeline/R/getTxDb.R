# Function to create a 'TxDb' transcript annotation object
# E.g., for passing to 'defineExomeVars' function

# SN UW GAC February 2015

# Arguments
# 1) genome build - "hg18" (build 36) or "hg19" (build 37)
# 2) database_table - name of database table to get:
# "refGene" (default) for  RefSeq "refGene" table, under track "Ref Seq Genes"
# "knownGene" for "knownGene" table, under track "UCSC genes"
# ...other tables listed by 'supportedUCSCtables'

# Output
# 'txdb' - a TxDb object ('GenomicFeatures' package)

# Notes - UCSC 'knownGene' table exists as the library 'TxDb.Hsapiens.UCSC.hg19.knownGene'
# For other database tables, including refGene, requires making the TxDb object

################# Start function definition

getTxDb <- function(build, database_table) {

  # check required packages
  require(GenomicFeatures)

  # check packages that may be required based on "database_table" argument
  if(database_table=="knownGene") {
      if(build=="hg19") {require(TxDb.Hsapiens.UCSC.hg19.knownGene)}
      if(build=="hg18") {require(TxDb.Hsapiens.UCSC.hg18.knownGene)}    
  }

  # check specified build
    if(!is.element(build, c("hg18","hg19")))
       {stop("Only implemented for human genome builds hg18 or hg19, not ", build)}  

  # check if specified database table is supported
  if(!is.element(database_table,
                 rownames(supportedUCSCtables()))){
    stop("Database table ", database_table," not found in list of 'supportedUCSCtables()'")
  }

  if(database_table=="knownGene") {
      message("Getting UCSC database table 'knownGene'")
      if(build=="hg19") {txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene}
      if(build=="hg18") {txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene}    
  }

  else {
    message("Fetching database table ", database_table)
    txdb <- makeTranscriptDbFromUCSC(genome = build, tablename = database_table)
}

  return(txdb)
}
