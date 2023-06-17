


dittplot<-function(gwas_genes,file_path,des_path){
    extract_file2<-function(pwd){
    files <-list.files(pwd)
    con<-list()
    lf0<-length(files)
    for( i in 1:lf0){
        xob<-read.csv(paste(pwd,files[i],sep="/"),row.names = NULL)
        need<-c('genes','logFC')
        con[[i]]<-xob[,which(colnames(xob)%in%need)]
    }
    file0<-gsub('.{4}$', '', files)
    aa<-list(file0,con)
        }
    gwas<-data.frame(genes=gwas_genes)
    # pwd='/Users/sam/Desktop/final_scrna_dge/test'
    con0<-list()
    con0[[1]]<- gwas
    file0<-extract_file2(file_path)
    con0[2:(length(file0[[2]])+1)]<-file0[[2]]
    df_list <- con0
    #merge all data frames together
    ax<- df_list %>% reduce(left_join, by='genes')
    colnames(ax) <- c('genes',file0[[1]])
    ax0<-ax[,-1]
    aaa<-which(rowSums(is.na(ax0))<dim(ax0)[2])
    ax00<-ax[aaa,]
    data_mod <- reshape2::melt(ax00, id.var = c('genes'))
    d <- data_mod
    d$genes <- data_mod[,1]
    d$contrast <- as.factor(data_mod[,2])
    ggplot(d, aes(contrast,genes ,  fill = value)) + geom_point(shape = 21) + theme_light() + guides(x =  guide_axis(angle = 90))
    ggsave( paste(des_path,"/fig_ditt_1.pdf",sep=""))
    data_mod2<-data_mod
    data_mod2$value2<- NA
    data_mod2$value2[data_mod2$value<0]<-'Blue'
    data_mod2$value2[data_mod2$value>0]<-'red'
    data_mod2<-data_mod2[,c(1,2,4)]
    d <- data_mod2
    d$genes <- data_mod2[,1]
    d$contrast <- as.factor(data_mod2[,2])
    # plot
    ggplot(d, aes(contrast,genes ,  fill = value2)) +geom_point(shape = 21) + theme_light() + guides(x =  guide_axis(angle = 90))
    ggsave( paste(des_path,"/fig_ditt_2.pdf",sep=""))
}

