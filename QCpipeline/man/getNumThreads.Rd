\name{getNumThreads}
\alias{getNumThreads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find the number of threads to use for multithreaded functions
}
\description{
Returns the number of cores to use based on user input and (if requested) environment variables set by SGE.
}
\usage{
getNumThreads(x = NA)
}

\arguments{
  \item{x}{ By default, returns 1; assumes no multithreading. If \code{x} or \code{strtoi(x)} is an integer and \code{x > 0}, it returns that integer. If \code{x = "n-m"}, where \code{n} and \code{m} are integers, it will detect the number of cores assigned to your job by SGE and return that number. }
}
\details{
This function was designed to choose the number of threads to run based on user input.

If \code{x = 1} (or any other integer), it will check to see if the number is a whole number and \code{> 0}. If so, it will return that number.

If \code{x = "n"}, getNumThreads will attempt to convert it to an integer and return that integer.

If \code{x = "n-m"}, where \code{"n"} and \code{"m"} are integers, it will check to see if \code{NSLOTS} has been set as an environment variable (by the SGE cluster scheduler).

In general, \code{"x"} must pass a number of checks to return a value other than 1.
}
\value{
An integer describing the number of threads to use, based on input to the function. Default is \code{1}.
}
\author{Adrienne Stilp}
\note{}
\seealso{}
\examples{
# returns the number of cores allocated to your job on the SGE cluster
getNumThreads("2-4")

# to specify a fixed number of threads regardless of the NSLOTS variable
# both will return 4
getNumThreads(4)
getNumThreads("4")

# All these return 1
getNumThreads()
getNumThreads(x=NULL)
getNumThreads(-1)
getNumThreads(2.4)
}
\keyword{}
