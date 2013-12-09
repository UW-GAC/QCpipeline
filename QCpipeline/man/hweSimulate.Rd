\name{hweSimulate}
\alias{hweSimulate}
\title{Simulations for HWE}
\description{Simulates a set of HWE results for a set of samples and SNPs, under the null hypothesis that the genotypes are in HWE.}
\usage{
hweSimulateGenotypeMatrix(nsnp, nsamp, aFreq)
hweSimulate(genotypes, p_in)
}
\arguments{
  \item{nsnp}{Number of SNPs to simulate}
  \item{nsamp}{Number of samples to simulate}
  \item{aFreq}{vector of allele A frequencies for the \code{nsnp} SNPs}
}
\details{}
\value{
\code{hweSimulateGenotypeMtarix} returns a random matrix of genotypes for \code{nsamp} samples and \code{nsnp} SNPs with input allele frequencies \code{aFreq}, assuming that the samples' genotypes are in HWE for each SNP. There are no missing calls.
\itemize{
\code{hweSimulate} returns a data frame with columns:
\item{ \code{aFreq.exp}: Expected (input) allele frequency}
\item{ \code{aFreq.sim}: Simulated allele frequency, using genotype/allele counts}
\item{ \code{maf.sim}: Simulated MAF, using genotype/allele counts}
\item{ \code{maf.exp}: Input MAF}
\item{ \code{f.sim}: Simulated inbreeding coefficient, using genotype/allele counts to determine expected number of heterozygotes}
\item{ \code{f.exp}: Simulated inbreeding coefficient, using input allele frequency to determine expected number of heterozygotes}
\item{ \code{pval.sim}: Simulated HWE p value, using genotype/allele counts to determine expected number of heterozygotes}
}
}
\author{Adrienne Stilp}
\seealso{}
\examples{
nsnp <- 100
nsamp <- 50
aFreq <- runif(nsnp)
genotypes <- hweSimulateGenotypeMatrix(nsnp, nsamp, aFreq)
sim <- hweSimulate(genotypes, aFreq)

plot(density(sim$f.sim, na.rm=TRUE), col="red")
lines(density(sim$f.exp), col="black")

}
\keyword{}
