#####
# Identify relatives from ibd coefficient estimates
#####


ibdAssignRelatednessPCRelate <- function(
  k0,  		# vector of k0 estimates
  kc,  			# vector of kinship coefficient estimates
  cut.kc.dup=1/(2^(3/2)),  # kinship coefficient threshold for duplicates
  cut.kc.fs=1/(2^(5/2)), # kc threshold for deg 1 relatives
  cut.kc.deg2=1/(2^(7/2)), # kc threshold for deg2 relatives
  cut.kc.deg3=1/(2^(9/2)), # kc threshold for deg3 relatives
  cut.k0=0.1 #
){
  
  ## returns a vector of assignments to "PO", "FS", "HS", "FC","U" (unrelated) and "Q" for everything else, for each pair of (k0,k1)
  
  ## default thresholds for assigning relationships use kinship coefficients in table 1 of Manichaikul (2010) - KING paper
  
  # if(length(k0)!=length(kc)) stop("k0 and kc must be parallel vectors of the same length")  
  
  
  dup.sel <- kc > cut.kc.dup
  po.sel <- kc <= cut.kc.dup & kc > cut.kc.fs & k0 <= cut.k0
  fs.sel <- kc <= cut.kc.dup & kc > cut.kc.fs & k0 > cut.k0
  d2.sel <- kc <= cut.kc.fs & kc > cut.kc.deg2
  d3.sel <- kc <= cut.kc.deg2 & kc > cut.kc.deg3
  un.sel <- kc <= cut.kc.deg3
  
  # check for overlap - should be none
  rels <- po.sel & dup.sel & fs.sel & d2.sel & d3.sel
  if (any(rels)) stop("one or more pairs assigned to more than one relationship")
  
  rels <- po.sel | dup.sel | fs.sel | d2.sel | d3.sel
  un.sel <- !rels 
  
  asnmt <- rep("Q", length(kc))
  asnmt[dup.sel] <- "Dup"
  asnmt[po.sel] <- "PO"
  asnmt[fs.sel] <- "FS"
  asnmt[d2.sel] <- "Deg2"
  asnmt[d3.sel] <- "Deg3"
  asnmt[un.sel] <- "U"
  
  return(asnmt)
  
}
