# Generate the distance between the gene of contrast 
dist_contrast<-function(con) {
 dist_const<-function(x,y){
  x2 <- list(first0<-x$genes,second0<-y$genes)
  v.table <- venn(x2)
  al<-v.table[1:4]
  return(1-al[4]/sum(al))
  }
  len0<-length(con)
  dist0<-matrix(NA, ncol=len0, nrow = len0)
  for(i in 1:len0){
    for(j in i:len0){
        if (i == j) next
      dist0[i,j]<- dist_const(con[[i]],con[[j]]) 
    }
  }
  dst<-data.matrix(as.dist(t(dist0)))
}

