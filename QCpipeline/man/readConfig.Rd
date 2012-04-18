\name{readConfig}
\alias{readConfig}
\alias{setConfigDefaults}
\title{Read a configuration file}
\description{Read a configuration file}
\usage{
readConfig(file, ...)

setConfigDefaults(config, required, optional, default)
}
\arguments{
  \item{file}{file where column 1 is parameter name
    and column 2 is value.}
  \item{...}{additional arguments to \code{read.table}}
  \item{config}{named character vector}
  \item{required}{character vector of required parameter names}
  \item{optional}{character vector of optional parameter names}
  \item{default}{character vector of default values for parameters in
    \code{optional}}
}
\details{}
\value{
  \code{readConfig} returns a named character vector of parameter
  values.

  \code{setConfigDefaults} takes a named character vector returned by
  \code{readConfig} and adds additional parameters in
  \code{optional} with values in \code{default}.
  An error will result if a parameter in \code{required} is
  missing from \code{config}.
}
\author{Stephanie Gogarten}
\seealso{}
\examples{
f <- tempfile()
write.table(cbind(letters[1:10], 1:10), file=f, quote=FALSE,
            row.names=FALSE, col.names=FALSE)
config <- readConfig(f)

required <- letters[1:5]
optional <- letters[11:15]
default <- 11:15
config <- setConfigDefaults(config, required, optional, default)

unlink(f)
}
\keyword{}
