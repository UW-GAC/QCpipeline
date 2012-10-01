
defineFamilies <- function(ibd, KC.threshold=0.1, start.fam.id=1) {
  # keep only pairs for which KC > KC.threshold
  dat <- ibd[ibd$KC > KC.threshold,]

  ids <- sort(unique(c(dat$sample1,dat$sample2))); 
  n <- length(ids)
  clust <- list()

  # make a cluster for each sample as all those with which it is related (including itself)
  for(i in 1:n){
    id <- ids[i]
    tmp <- dat[is.element(dat$sample1, id)  | is.element(dat$sample2, id), c("sample1", "sample2")]
    clust[[i]] <- sort(union(tmp$sample1,tmp$sample2))
  }
		
  # check each cluster against all the others and merge any that have overlaps
  cl <- list()
  for(i in 1:n){
    cl[[i]] <- clust[[i]]
    for(j in 1:n){
      if(i==j) next
      # merge two clusters if they overlap
      if(any(is.element(cl[[i]], clust[[j]]))) cl[[i]] <- sort(union(cl[[i]], clust[[j]]))
    }
  }

  # remove duplicate clusters
  clust <- cl
  cl <- list()
  cl[[1]] <- clust[[1]]
  k <- 1
  for(i in 2:n){
    tmp <- clust[[i]]
    chk <- rep(NA, k)
    for(j in 1:k) { # for each unique cluster stored in cl
      if(length(tmp) <= length(cl[[j]])) { chk[j] <- all(tmp %in% cl[[j]]) }  else { chk[j] <- FALSE }  
    }
    if(!any(chk)) { k <- k + 1; cl[[k]] <- tmp }
  }
  
  # check
  chk <- unlist(cl)
  stopifnot(length(chk) == unique(length(chk)))
  tmp <- unlist(lapply(cl, length))
  stopifnot(all(is.element(c(dat$sample1,dat$sample2),chk)))

  n <- length(cl)
  dat$family <- NA
  for(i in 1:n) dat[is.element(dat$sample1,cl[[i]]) & is.element(dat$sample2,cl[[i]]), "family"] <- i+start.fam.id-1
  
  return(dat)
}
