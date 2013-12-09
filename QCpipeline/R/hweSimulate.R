# counts genotypes in the same way as the HWE GWASTools function
.countGenotypes <- function(genotypes) {
  nAA <- apply(genotypes, 1, function(x) sum(x == 2, na.rm=TRUE))
  nAa <- apply(genotypes, 1, function(x) sum(x == 1, na.rm=TRUE))
  naa <- apply(genotypes, 1, function(x) sum(x == 0, na.rm=TRUE))
  return(data.frame(nAA, nAa, naa))
}

# counts alleles in the same way as the HWE GWASTools function
.countAlleles <- function(gc) {
  nA <- 2 * gc$nAA + gc$nAa
  na <- 2 * gc$naa + gc$nAa
  return(data.frame(nA, na))
}

# simulates a genotype matrix with a given snp and sample size, as well as a given allele frequency for each snp
hweSimulateGenotypeMatrix <- function(nsnp, nsamp, aFreq) {
  # genotypes. hm.
  tmp.1 <- scale(matrix(runif(nsnp*nsamp), nrow=nsamp, ncol=nsnp), center=aFreq)
  tmp.2 <- scale(matrix(runif(nsnp*nsamp), nrow=nsamp, ncol=nsnp), center=aFreq)
  
  alle.1 <- ifelse(tmp.1 < 0, 1,
                   ifelse(tmp.1 > 0, 0, sample(c(0,1),1)))
  alle.2 <- ifelse(tmp.2 < 0, 1,
                   ifelse(tmp.2 > 0, 0, sample(c(0,1),1)))
  
  geno <- t(alle.1 + alle.2)
  
  geno
  
}

# returns simulated HWE results
hweSimulate <- function(genotypes, p_in){
  
  genotypeCounts <- .countGenotypes(genotypes)
  alleleCounts <- .countAlleles(genotypeCounts)
  aFreq <- alleleCounts[,1]/(as.matrix(alleleCounts) %*% cbind(c(1, 1)))
  
  obs.het <- genotypeCounts[,2]
  geno.tot <- apply(genotypeCounts, 1, sum)
  exp.het <- 2 * aFreq * (1 - aFreq) * geno.tot
  f_obs <- 1 - (obs.het/exp.het)
  
  f_exp <- 1 - (obs.het / (2 * p_in * (1 - p_in) * geno.tot))
  
  hwePs <- GWASExactHW::HWExact(genotypeCounts)
  
  return(data.frame(f.sim=f_obs, f.exp=f_exp, aFreq.sim=aFreq, aFreq.exp=p_in,
                    maf.exp=ifelse(p_in < 0.5, p_in, 1-p_in),
                    maf.sim=ifelse(aFreq < 0.5, aFreq, 1-aFreq),
                    pval.sim=hwePs))
  
}