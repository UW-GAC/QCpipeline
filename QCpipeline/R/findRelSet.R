findRelSet <-
function(kinMat, thresh, unrel.set=NULL){
	
	# set diagonal to 0
	diag(kinMat) <- 0
	
	# check if matrix is symmetric
	if(!all(kinMat == t(kinMat))){
		stop("kinMat must be a Symmetric Matrix!")
	}
	
	# number of divergent pairs
	ndiv <- apply(kinMat,1,function(x){ sum(x < -thresh) })
	
	# "total kinship" values
	kinMat[kinMat < thresh] <- 0
	kinsum <- rowSums(kinMat)
	# indicator for relatives
	indKIN <- kinMat != 0
	
	rels <- NULL
	if(!is.null(unrel.set)){
		for(i in 1:length(unrel.set)){
			# identify relatives of individual set to the unrelated set
			rel.idx <- which(indKIN[unrel.set[i],])
			# don't count relatives that are also set to the unrelated set
			rel.idx <- rel.idx[!is.element(rel.idx,unrel.set)]
			rels <- append(rels, rel.idx)
			indKIN[rel.idx,] <- FALSE; indKIN[,rel.idx] <- FALSE
		}
		indKIN[unrel.set,unrel.set] <- FALSE
	}

    # number of relatives	
	numrel <- rowSums(indKIN)
	while(max(numrel) > 0){
		idx <- which(numrel == max(numrel))
		if(length(idx) > 1){
			ndiv2 <- ndiv[idx]
			didx <- which(ndiv2 == min(ndiv2))
			if(length(didx) == 1){
				idx <- idx[didx]
			}else{
				kinsum2 <- kinsum[idx[didx]]
				kidx <- which(kinsum2 == min(kinsum2))[1]
				idx <- idx[didx[kidx]]
			}
		}
		rels <- append(rels, idx)
		indKIN[idx,] <- FALSE; indKIN[,idx] <- FALSE;
		numrel <- rowSums(indKIN)
	}
	
	if(!is.null(rels)){
		rels <- as.numeric(sort(rels))
	}
	return(rels)
}
