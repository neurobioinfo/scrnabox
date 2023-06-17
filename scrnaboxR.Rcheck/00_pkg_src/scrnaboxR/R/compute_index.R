# Compute the indices 
compute_index <-function(X,Y) {
con1<-table(X,Y)
library(aricode)
library(mclust)
avcorrec<-function(trueclu,obsclus)  1-classError(trueclu, obsclus)$errorRate
CramerVF<-function(con1){
    sqrt(mean(chisq.test(con1)$statistic)/(min((dim(con1)[1]-1),(dim(con1)[2]-1))*sum(con1)))
    # sqrt(chisq.test(con1)$statistic/(chisq.test(con1)$statistic+n))
}

return(list(
chisq_test=chisq.test(con1), 
contigency_table=con1,
prop_table=prop.table(con1), 
prop_table_row=prop.table(con1, margin = 1), 
prop_table_col=prop.table(con1, margin = 2), 
cramerv=CramerVF(con1), 
accuracy=avcorrec(X,Y),
adjusted_mutual_information=AMI(X,Y),
normalized_mutual_information=NMI(X,Y),
rand_index=RI(X,Y),
adjusted_rand_index=ARI(X,Y)
))
}
