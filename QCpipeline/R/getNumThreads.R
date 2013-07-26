# function to return the number of threads

# function to check if x is a valid positive integer
checkPositiveInt <- function(x){

  if (is.null(x)) return(FALSE)
  if (!is.numeric(x)) return(FALSE)
  # covers NA, NaN, Inf
  if (!is.finite(x)) return(FALSE)
  if (x <= 0) return(FALSE)
  if (abs(x %% 1) > .Machine$double.eps^0.5) return(FALSE)

  # made it past all the checks.
  return(TRUE)
}

# for testing
checkPositiveIntTest <- function() {
  stopifnot(!checkPositiveInt(NA))
  stopifnot(!checkPositiveInt(NULL))
  stopifnot(!checkPositiveInt(checkPositiveInt))
  stopifnot(!checkPositiveInt(0))
  stopifnot(!checkPositiveInt(-1))
  stopifnot(!checkPositiveInt(-1.0))
  stopifnot(!checkPositiveInt(-0.5))
  stopifnot(!checkPositiveInt(-Inf))
  stopifnot(!checkPositiveInt(Inf))
  stopifnot(!checkPositiveInt(NaN))
  stopifnot(!checkPositiveInt("-1"))  
  stopifnot(!checkPositiveInt("1"))
  stopifnot(!checkPositiveInt(3.5))
  
  stopifnot(checkPositiveInt(3))
  
}

# function to return the number of threads to for in multithreading
getNumThreads <- function(x=NA){
  # valid options are 
  # - n-m (where n, m are numbers; in which case it will check to see how many slots have been assigned by SGE)
  # - n (in which case it will just set the number of threads to n)
  #
  # default is 1
  
  num_threads <- 1
  
  # this will be a string
  if (is.character(x)) {
    
    x_int <- strtoi(x)
    
    if(checkPositiveInt(x_int)) {
      num_threads <- x_int
    } else {
      
      x_strspl <- strsplit(x, "-")
      if (length(x_strspl) == 2){
        if (checkPositiveInt(x_strspl[1]) & checkPositiveInt(x_strpl[2])){
          checkSlots <- strtoi(Sys.getenv("NSLOTS"))
          
          if (checkPositiveInt(checkSlots)) {
            num_threads <- checkSlots
          }
        }
      }
    }
  } else if (is.numeric(x)){
    if (checkPositiveInt(x)){
      num_threads <- x
    }
  }
  
  return(num_threads)
}

