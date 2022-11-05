separ0<-function(file0,dest_path,xob,cut_pvalue,cut_log_low,cut_log_up,n_row){
  xob2<-xob
  colnames(xob2)[1]<-"genes"
  # colnames(xob2)<-c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj')
  # rownames(xob2)<-xob$X
  fil2<-xob2[xob2$adj.P.Val< cut_pvalue,]
  filn<-fil2[fil2$logFC< cut_log_low,]
  filp<-fil2[fil2$logFC > cut_log_up,]
  filn2<-filn[1:n_row,]
  filp2<-filp[1:n_row,]
  file0<-gsub('.{4}$', '', file0)
  write.csv(filn2,file=paste(paste(dest_path,file0,sep='/'),'down.csv',sep=''),row.names = FALSE)
  write.csv(filp2,file=paste(paste(dest_path,file0,sep='/'),'up.csv',sep=''),row.names = FALSE)
}

split_up_down<-function(file_path,dest_path,cut_pvalue,cut_log_low=0,cut_log_up=0, n_row){
  files <-list.files(file_path,pattern="*.csv")
  lf0<-length(files)
  if (missing(cut_pvalue)) {
    cut_pvalue=1.0
  } 
  for( i in 1:lf0){
      xob<-read.csv(paste(file_path,files[i],sep="/"))
      if (missing(n_row)) {
        n_row=nrow(xob)
      } 
      separ0(file0=files[i],dest_path,xob,cut_pvalue,cut_log_low,cut_log_up,n_row)
    }
}

