extract_file<-function(pwd){
files <-list.files(pwd)
con<-list()
lf0<-length(files)
for( i in 1:lf0){
    con[[i]]<-read.csv(paste(pwd,files[i],sep="/"),row.names = NULL)
}
aa<-list(files,con)
}

