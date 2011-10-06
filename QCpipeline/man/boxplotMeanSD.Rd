\name{boxplotMeanSD}
\alias{boxplotMeanSD}
\title{Boxplot with mean and SD}
\description{Boxplot with mean and SD}
\usage{
boxplotMeanSD(x, y, data=NULL, xlab=NULL, ylab=NULL, nSD=1, ...)
}
\arguments{
  \item{x}{vector or character string denoting column in data}
  \item{y}{vector or character string denoting column in data}
  \item{data}{data.frame}
  \item{xlab}{title for x axis (defaults to x if data is not \code{NULL}}
  \item{ylab}{title for y axis (defaults to y if data is not \code{NULL}}
  \item{nSD}{number of standard deviations to plot}
  \item{...}{additional plotting arguments}
}
\details{}
\author{Jess Shen}
\seealso{}
\examples{
age <- sample(25:55, 100, replace=TRUE)
sex <- sample(c("M", "F"), 100, replace=TRUE)
boxplotMeanSD(sex, age)

data <- data.frame(age, sex)
boxplotMeanSD("sex", "age", data)
}
\keyword{}
