##
setwd("~/Dropbox/structure_simulations/satrapa/")
source("~/Dropbox/structure_simulations/R_scripts/ggthemes.R")
library(tidyr);library(plyr);library(adegenet);library(data.table);library(foreach);library(ggjoy)
library(gridExtra)
#library(doMC)
#registerDoMC(cores=8)

#satrapa pops c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,2,1,1,1)


#mode function for checking k-means cluster assignments (via https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#function to run 100 k-means assignments and dapc x-validations on a designated infile (single-threaded)
cluster_multivar <- function(infile,pop.info=F,pop){
  tmp <- read.table(infile)
  nloci <- ncol(tmp)-(1+pop.info)
  nind <- nrow(tmp)/2
  if(pop.info==T){
    str <- read.structure(infile,n.ind=nind,n.loc=nloci,col.lab=1,col.pop=2,row.marknames=0,onerowperind=F,ask=F)
  } else {
    str <- read.structure(infile,n.ind=nind,n.loc=nloci,col.lab=1,col.pop=0,row.marknames=0,onerowperind=F,ask=F)
    str@pop <- factor(pop)
  }
  #k-means cluster, find % assigned to correct cluster after accounting for label switching
  kmeans_accuracy <- vector(length=100)
  for(i in 1:100){
    clust <- find.clusters(str,n.pca=length(str@pop),n.clust=3)$grp
    
    popclust <- c(Mode(clust[str@pop==1]),Mode(clust[str@pop==2]),Mode(clust[str@pop==3])) 
    if(length(unique(popclust)) < 3){ #if any pop lacks a unique modal cluster assignment, correct popclust
      missingclust <- c(1:3)[!c(1:3) %in% popclust]
      popclustct <- c(sum(as.integer(clust[str@pop==1]==popclust[1])),
              sum(as.integer(clust[str@pop==2]==popclust[2])),
              sum(as.integer(clust[str@pop==3]==popclust[3])))
      popclust[which.min(popclustct)] <- missingclust
    }
    kmeans_accuracy[i] <- sum(sum(as.integer(clust[str@pop==1]==popclust[1])),
                    sum(as.integer(clust[str@pop==2]==popclust[2])),
                    sum(as.integer(clust[str@pop==3]==popclust[3])))/length(str@pop)
  }
  
  #dapc using original (pre-defined) clusters, using n PC's w/highest success in first 10 reps. 
  str_noMD <- scaleGen(str,NA.method="mean",center=T,scale=F)
  xval <- xvalDapc(str_noMD,grp=str@pop,n.da=2,n.pca.max=length(str@pop),training.set=0.5,center=F,scale=F,xval.plot=F,n.rep=10)
  xval_npc <- xval$`Number of PCs Achieving Highest Mean Success` %>% as.integer()
  xval <- xvalDapc(str_noMD,grp=str@pop,n.da=2,n.pca=xval_npc,center=F,scale=F,xval.plot=F,n.rep=100,training.set=0.5)
  xval <- xval$`Cross-Validation Results`[,2]
  
  out <- data.frame(kmeans=kmeans_accuracy,xval=xval,file=infile)
  out
}

#get list of files to cluster (edit as needed before running)
files <- list.files("str_in",full.names = T)
files <- grep("maf",files,value=T)
files <- files[1:7]

out <- foreach(i=files,.combine = rbind) %dopar% cluster_multivar(i,pop.info = F,pop=c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,2,1,1,1))

#summarize across runs/maf levels (simulated data)
out <- read.csv("tmp.csv")
out$maf <- strsplit(as.character(out$file),'maf') %>% lapply(function(e) as.integer(unlist(strsplit(unlist(e)[2],"\\."))[1])) %>% unlist()

#maf column for satrapa data
out$maf <- c(rep(0,100),rep(0.0125,100),rep(0.025,100),rep(0.0375,100),rep(0.05,100),rep(0.1,100),rep(0.25,100))

### plotting (EBL edits)
out <- read.csv("~/Dropbox/structure_simulations/clust_multivar_out.csv")
files <- as.character(out$file)
maf <- lapply(files, function(i){
  gsub(pattern = "(.*maf)(.*)(.str)", replacement = "\\2", x = i)
})
out$maf <- as.numeric(unlist(maf))

plot1 <- ggplot(data=out,aes(x=kmeans,y=factor(maf)))+
  theme_minimal()+
  xlim(0.25,1.1)+
  ggtitle("K-means cluster\n assignment accuracy")+
  ylab("Minimum Minor Allele Count")+xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy(bandwidth=0.025)

plot2 <- ggplot(data=out,aes(x=xval,y=factor(maf)))+
  theme_minimal()+
  xlim(0.25,1.1)+
  ggtitle("DAPC x-validation\n accuracy")+
  ylab("Minimum Minor Allele Count")+xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy(bandwidth=0.03)

out_reg <- read.csv("./satrapa/satrapa_cluster_multivar_out.csv")
out_reg$maf <- as.factor(out_reg$maf)
levels(out_reg$maf) <- c("0","1","2","3","4","8","20")

plot3 <- ggplot(data=out_reg,aes(x=kmeans,y=factor(maf)))+
  theme_minimal()+
  xlim(0.0,1.1)+ 
  #tweaked scale across all to avoid clipping
  #ggtitle("K-means cluster\n assignment accuracy")+
  ylab("Minimum Minor Allele Count")+xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy(bandwidth=0.025)

plot4 <- ggplot(data=out_reg,aes(x=xval,y=factor(maf)))+
  theme_minimal()+
  xlim(0.25,1.1) +
  #ggtitle("DAPC x-validation\n accuracy")+
  ylab("Minimum Minor Allele Count")+xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy(bandwidth=0.03)

### PCA plots

# list files
files <- list.files("./satrapa/str_in/", pattern = ".str", full.names = TRUE)
files <- files[1:7]

# count loci
loci <- lapply(files, function(x){
  ncol(read.table(x))-1 #subtract sample ID column
})
names(loci) <- files

#ID clusters
maf0 <- read.structure("./satrapa/str_in/satrapa_unlinked.str", n.ind=33, n.loc=3898, onerowperind = F,col.lab=1,col.pop=0, row.marknames=0,NA.char="-9",ask=F)
k3.maf0 <- find.clusters(maf0, n.pca=95,n.clust = 3,choose.n.clust = F)

# read each file, run & plot PCA
plots <- list()
maf <- c(0,1,2,3,4,8,20)
for(i in 1:length(loci)){
  seq <- read.structure(names(loci[i]), n.ind=33, n.loc=loci[[i]], onerowperind = F,col.lab=1,col.pop=0,row.marknames=0,NA.char="-9",ask=F)
  text <- paste("MAC=",maf[i],sep="")
  seq.scaled <- scaleGen(seq,NA.method="mean",scale=F)
  pca <- prcomp(seq.scaled,center=F,scale=F)
  pc <- data.frame(pca$x[,1:3])
  p <- ggplot(pc,aes(x=PC1,y=PC2,col=k3.maf0$grp))+xlab("PC1")+ylab("PC2")+
    scale_color_grey() +
    theme(panel.border = element_rect(color="black",fill=NA),
          legend.position="none",
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    geom_point(shape=k3.maf0$grp, size=1.5,stroke=1) +
    scale_shape_manual(values = c(15,16,17)) +
    ggtitle(text)
  print(p)
  plots <- c(plots, list(p))
}


# list files
simfiles <- list.files("~/Dropbox/structure_simulations/str_in/k3", pattern = "sim01.*", full.names = TRUE)

# count loci
simloci <- lapply(simfiles, function(x){
  ncol(read.table(x))-2 #subtract sample ID column + pop ID
})



names(simloci) <- simfiles

#ID clusters
sim.maf0 <- read.structure("/Users/ethanlinck/Dropbox/structure_simulations/str_in/k3/sim01.str" , n.ind=30, n.loc=4530, onerowperind = F,col.lab=1,col.pop=1, row.marknames=0,NA.char="-9",ask=F)
sim.k3.0 <- find.clusters(sim.maf0, n.pca=95,n.clust = 3,choose.n.clust = F)

# read each file, run & plot PCA
simplots <- list()
simmaf <- c(0,1,2,20,3,4,5,8)
for(i in 1:length(simloci)){
  sim.seq <- read.structure(names(simloci[i]), n.ind=30, n.loc=simloci[[i]], onerowperind = F,col.lab=1,col.pop=1,row.marknames=0,NA.char="-9",ask=F)
  sim.text <- paste("MAC=",simmaf[i],sep="")
  sim.seq.scaled <- scaleGen(sim.seq,NA.method="mean",scale=F)
  sim.pca <- prcomp(sim.seq.scaled,center=F,scale=F)
  sim.pc <- data.frame(sim.pca$x[,1:3])
  x <- ggplot(sim.pc,aes(x=PC1,y=PC2,col=sim.k3.0$grp))+xlab("PC1")+ylab("PC2")+
    scale_color_grey() +
    theme(panel.border = element_rect(color="black",fill=NA),
          legend.position="none",
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    geom_point(size=1.5,shape=sim.k3.0$grp,stroke=1,alpha=0.9) +
    ggtitle(sim.text)
  print(x)
  simplots <- c(simplots, list(x))
}

#side-by-side plots, 3 PCAs
layout <- cbind(c(1,1,1),c(1,1,1),c(2,2,2),c(2,2,2),c(3,4,5))
r2 <- cbind(c(6,6,6),c(6,6,6),c(7,7,7),c(7,7,7),c(8,9,10))
layout <- rbind(layout,r2)
grid.arrange(plot1,plot2,plots[[7]],plots[[5]],plots[[1]],plot3,plot4,simplots[[4]],simplots[[6]],simplots[[1]],layout_matrix=layout)


