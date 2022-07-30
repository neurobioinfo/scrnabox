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


