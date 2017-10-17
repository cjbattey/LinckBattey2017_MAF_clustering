#faststructure results analysis
setwd("~/Dropbox/structure_simulations/")
source("./R_scripts/ggthemes.R");source("./R_scripts/structurePlot.R")
library(ggplot2);library(plyr);library(magrittr);library(data.table);library(ggjoy);library(gridExtra)

####################################
######### simulated data ###########
####################################
#index of files by maf
files <- list.files("fstr_out",full.names = T)
files <- grep("meanQ",files,value = T)
maf_index <- list()
#maf_levels <- c(0.0,0.0125,0.025,0.0375,0.05,0.1,0.25)
maf_levels <- c(1,2,3,4,5,8,20)
b <- 1
for(i in maf_levels){
  maf_index[[b]] <- grep(paste0("maf",i,"."),files,fixed=T,value=T)
  b <- b+1
}

#read in q-matrices
q_matrices <- list()
a <- 1
for(i in maf_index){
  for(j in 1:100){
    fstr <- fread(i[j],data.table = F)
    fstr$pop <- c(rep(1,10),rep(2,10),rep(3,10))
    fstr$id <- 1:30
    
    #swap column order to minimize label switching
    tmp <- ddply(fstr,.(pop),function(e) colMeans(e[,1:3]))
    popclust <- c(which.max(tmp[1,2:4]),which.max(tmp[2,2:4]),which.max(tmp[3,2:4]))
    if(length(unique(popclust))!=3){
      popclust[duplicated(popclust)] <- c(1,2,3)[c(1,2,3) %in% c(popclust) == F]
    }
    fstr <- fstr[,c(popclust,4,5)]
    colnames(fstr) <- c("q1","q2","q3","pop","id")
    
    q_matrices[[a]] <- fstr
    a <- a+1
  }
}
names(q_matrices) <- unlist(maf_index)

#get distance bw populations in q-matrix space
mean_q_list <- lapply(q_matrices,function(e) {
  ddply(e,"pop",summarize,q1=mean(q1),q2=mean(q2),q3=mean(q3))[-1]
})
q_dist_list <- sapply(mean_q_list,function(e){
  mean(dist(e))
})

#df for plotting 
fstr_out <- data.frame(maf=unlist(lapply(maf_levels,function(e) rep(e,100))),popclustdist=q_dist_list)
#write.table(fstr_out,"fstr_out_summary.csv",sep=",",row.names = F)

#joyplot
joyplot <- ggplot(data=fstr_out,aes(x=popclustdist,y=factor(maf)))+
  theme_minimal()+
  #ggtitle("Faststructure Population Discrimination\nby Minimum Minor Allele Count")+
  xlab("Population Discrimination")+ylab("Minimum Minor Allele Count")+
  scale_y_discrete(expand=c(0.01,0))+xlim(0,1.5)+
  geom_joy()

#structure plots
q <- do.call(rbind,q_matrices)
q$maf <- unlist(lapply(rownames(q),function(e) strsplit(e,"maf") %>% unlist() %>% .[2] %>% strsplit(.,"\\.") %>% unlist() %>% .[1]))
q$run <- unlist(lapply(rownames(q), function(e) strsplit(e,"sim") %>% unlist() %>% .[2] %>% strsplit(.,"\\.") %>% unlist() %>% .[1]))
q$maf <- factor(q$maf,levels(factor(q$maf))[c(3,7,6,5,4,2,1)])
meltq <- melt(q,id.vars=c("id","pop","maf","run"))
structure_plot <- ggplot(data=subset(meltq,run=="03"),aes(x=id,y=value,fill=variable))+facet_grid(maf~.)+
  theme_minimal()+theme(axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        axis.text.x=element_text(angle=45,hjust=1,size=7.5))+
  ylab("")+
  scale_fill_manual(values = grey.colors(3))+
  geom_bar(stat="identity")

layout <- matrix(c(1,1,1,1,2,2,2,
                   1,1,1,1,2,2,2,
                   1,1,1,1,2,2,2),nrow=3,byrow=T)

grid.arrange(joyplot,structure_plot,layout_matrix=layout)

####################################################
########## empirical (R. satrapa) data #############
####################################################
setwd("~/Dropbox/structure_simulations/satrapa")
files <- list.files("fstr_out",full.names = T)
files <- grep("meanQ",files,value = T)
maf_index <- list()
maf_levels <- c(0.0125,0.025,0.0375,0.05,0.1,0.25)
#maf_levels <- c(1,2,3,4,5,8,20)
b <- 2
for(i in maf_levels){
  maf_index[[b]] <- grep(paste0("maf",i,"_"),files,fixed=T,value=T)
  b <- b+1
}
maf_index[[1]] <- files[1:10]

#read in q-matrices
names <- readLines("satrapa_sampleID_rowOrder.txt")
q_matrices <- list()
a <- 1
for(i in maf_index){
  for(j in 1:10){
    fstr <- fread(i[j],data.table = F)
    fstr$pop <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,2,1,1,1)
    fstr$id <- names
    fstr <- arrange(fstr,desc(pop))
    #swap column order to minimize label switching
    tmp <- ddply(fstr,.(pop),function(e) colMeans(e[,1:3]))
    popclust <- c(which.max(tmp[1,2:4]),which.max(tmp[2,2:4]),which.max(tmp[3,2:4]))
    if(length(unique(popclust))!=3){
      popclust[duplicated(popclust)] <- c(1,2,3)[c(1,2,3) %in% popclust == F]
    }
    fstr <- fstr[,c(popclust,4,5)]
    colnames(fstr) <- c("q1","q2","q3","pop","id")
    
    q_matrices[[a]] <- fstr
    a <- a+1
  }
}
names(q_matrices) <- unlist(maf_index)

#get Euclidean distance bw populations in q-matrix space
mean_q_list <- lapply(q_matrices,function(e) {
  ddply(e,"pop",summarize,q1=mean(q1),q2=mean(q2),q3=mean(q3))[-1]
})
q_dist_list <- sapply(mean_q_list,function(e){
  mean(dist(e))
})

#df for plotting 
fstr_out <- data.frame(maf=c(rep(0,10),unlist(lapply(maf_levels,function(e) rep(e,10)))),popclustdist=q_dist_list)
#write.table(fstr_out,"fstr_out_satrapa_summary.csv",sep=",",row.names = F)

#joyplot
joyplot <- ggplot(data=fstr_out,aes(x=popclustdist,y=factor(maf)))+
  theme_minimal()+
  #ggtitle("Faststructure Population Discrimination\nby Minimum Minor Allele Count")+
  xlab("Population Discrimination")+ylab("Minimum Minor Allele Count")+
  scale_y_discrete(expand=c(0.01,0))+
  geom_joy()

#structure plots (need to fix label switching in read-in loop for matching colors to populations)
q <- do.call(rbind,q_matrices)
q$maf <- c(rep(0,10*33),unlist(lapply(maf_levels,function(e) rep(e,10*33))))
q$maf <- factor(q$maf,levels(factor(q$maf))[c(7,6,5,4,3,2,1)])
q$run <- unlist(lapply(rownames(q), function(e) strsplit(e,"_") %>% unlist() %>% .[length(.)] %>% strsplit("\\.") %>% unlist() %>% .[1])) #I'm sorry
q$id <- factor(q$id,levels(factor(q$id))[c(18:21,1:17,26:30,22:25,31:33)]) #order factor levels to group populations (why is this so difficult?)
meltq <- melt(q,id.vars=c("id","pop","maf","run"))
structure_plot <- ggplot(data=subset(meltq,run==2),aes(x=id,y=value,fill=variable))+facet_grid(maf~.)+
  theme_minimal()+theme(axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        axis.text.x=element_text(angle=45,hjust=1,size=6))+
  ylab("")+xlab("")+
  scale_fill_manual(values = grey.colors(3))+
  geom_bar(stat="identity")

layout <- matrix(c(1,1,1,1,2,2,2,
                   1,1,1,1,2,2,2,
                   1,1,1,1,2,2,2),nrow=3,byrow=T)

grid.arrange(joyplot,structure_plot,layout_matrix=layout)
