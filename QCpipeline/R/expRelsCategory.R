expRelsCategory<-function(pedigree){
 # input pedigree
 # assign relationship categories to inbred as well as other relations
 # use KING thresholds
 # this is for ease of comparison with expected KC

  prwise<-pedigreePairwiseRelatedness(pedigree)
  cut.dup <- 1/(2^(3/2))
  cut.deg1 <- 1/(2^(5/2))
  cut.deg2 <- 1/(2^(7/2))
  cut.deg3 <- 1/(2^(9/2))

  kc<-NULL
# assign relationship categories to inbred
  inb.fam<-prwise$inbred.fam 
  if(!is.null(inb.fam)){ 
     kc<-prwise$inbred.KC 
     unrel<-kc$kinship<=cut.deg3
     deg3<-kc$kinship>cut.deg3 & kc$kinship<=cut.deg2
     deg2<-kc$kinship>cut.deg2 & kc$kinship<=cut.deg1
     fspo<-kc$kinship>cut.deg1 & kc$kinship<=cut.dup  # analyzed further below
     dup<-kc$kinship>cut.dup

     kc$exp.rel<-NA
     kc$exp.rel[unrel]<-"U"
     kc$exp.rel[deg3]<-"Deg3"
     kc$exp.rel[deg2]<-"Deg2"
     kc$exp.rel[dup]<-"Deg0"

     chk<-kc[fspo,]
     fams<-unique(chk$family)
     for(f in fams){
       tmp<-chk[chk$family==f,]
       pp<-pedigree[pedigree$family==f,]
       for(j in 1:nrow(tmp)){
          id1<-tmp$Individ1[j]
          id2<-tmp$Individ2[j]
          p1<-pp[pp$individ==id1,]
          p2<-pp[pp$individ==id2,]
          if(id1 %in% c(p2$mother,p2$father) | id2 %in% c(p1$mother,p1$father)){
            sel<-kc$Individ1 %in% c(id1,id2) & kc$Individ2 %in% c(id1,id2)
            kc$exp.rel[sel]<-"PO"
          } else {
          if(p1$mother==p2$mother & p1$father==p2$father){
             sel<-kc$Individ1 %in% c(id1,id2) & kc$Individ2 %in% c(id1,id2)
             kc$exp.rel[sel]<-"FS"
          }
          }
       }
     }
  kc$exp.rel[is.na(kc$exp.rel)]<-"Deg1" # larger KC than Deg2 but not FS or PO
  kc$relation<-"inbred family" 
  }
 
# assign relationship categories to relative pairs
  relprs<-prwise$relativeprs 
  if(!is.null(relprs)){
     relprs$exp.rel<-NA
     unrel<-relprs$kinship<=cut.deg3
     deg3<-relprs$kinship>cut.deg3 & relprs$kinship<=cut.deg2
     deg2<-relprs$kinship>cut.deg2 & relprs$kinship<=cut.deg1
     fspo<-relprs$kinship>cut.deg1 & relprs$kinship<=cut.dup  
     dup<-relprs$kinship>cut.dup
     relprs$exp.rel<-NA
     relprs$exp.rel[unrel]<-"U"
     relprs$exp.rel[deg3]<-"Deg3"
     relprs$exp.rel[deg2]<-"Deg2"
     relprs$exp.rel[dup]<-"Dup"
     sel<-fspo & relprs$relation %in% c("PO","FS")
     relprs$exp.rel[sel]<-relprs$relation[sel]
     relprs$exp.rel[is.na(relprs$exp.rel)]<-"Deg1"  
   }

  relprs.all<-rbind(relprs,kc)
 
  out<-list(inb.fam,relprs.all)
  names(out)<-c("inbred.fam","relprs.all")
  return(out)
}  
  
