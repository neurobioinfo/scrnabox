msg <- function(x, startup = FALSE) {
  # the following function are from tidyverse
  if (startup) {
    if (!isTRUE(getOption("scrnaboxR.quiet"))) {
      rlang::inform(x, class = "packageStartupMessage")
    }
  } else {
    rlang::inform(x)
  }
}

extract_file<-function(pwd){
files <-list.files(pwd)
con<-list()
lf0<-length(files)
for( i in 1:lf0){
    con[[i]]<-read.csv(paste(pwd,files[i],sep="/"),row.names = NULL)
}
aa<-list(files,con)
}



# library("xlsx")
csv_to_xlsx<-function(files_path,dest_path,name='sheet'){
  files <-list.files(files_path,pattern="*.csv")
  xlsx::write.xlsx(files, file=("final.xlsx"),sheetName=paste0('list'), row.names=TRUE,append=TRUE)
  lf0<-length(files)
  for( i in 1:lf0){
      xob<-read.csv(paste(files_path,files[i],sep="/"))
      write.xlsx(xob, file=paste("final.xlsx"),sheetName=paste0('sheet',i), row.names=TRUE,append=TRUE)
    }
}


