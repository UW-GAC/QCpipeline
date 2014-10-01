test_calculateLambda <- function(){
  
  # df=1
  stat <- rchisq(1e6, df=1)
  stat[sample(10e6, 100)] <- NA
  pval <- pchisq(stat, df=1, lower.tail=F)  
  
  lambda <- median(stat, na.rm=T) / qchisq(0.5, df=1)
  checkEquals(lambda, calculateLambda(stat, df=1))

  # df=2
  stat <- rchisq(1e6, df=2)
  stat[sample(10e6, 100)] <- NA
  pval <- pchisq(stat, df=2, lower.tail=F)  
  
  lambda <- median(stat, na.rm=T) / qchisq(0.5, df=2)
  checkEquals(lambda, calculateLambda(stat, df=2))

  n <- sample(3:10, 1)
  stat <- rchisq(1e6, df=n)
  stat[sample(10e6, 100)] <- NA
  pval <- pchisq(stat, df=n, lower.tail=F)  
  
  lambda <- median(stat, na.rm=T) / qchisq(0.5, df=n)
  checkEquals(lambda, calculateLambda(stat, df=n))
  
}