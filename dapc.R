### discriminant analysis of principal components

library(rgdal);library(adegenet);library(plyr);library(ggplot2);
library(viridis);library(grid);library(magrittr);library(raster);
library(foreach);library(stringr);library(data.table);library(broom);
library(gridExtra);library(dplyr);library("ggthemes");library(ggjoy);
library(devtools);library(cleangeo);library(mapdata);library(maps)

setwd("~/Dropbox/structure_simulations/satrapa/")

localities <- read.csv("satrapa_locations.csv", stringsAsFactors = FALSE)

# list files
files <- list.files("./str_in", pattern = ".str", full.names = TRUE)
files <- files[1:7]

# count loci
loci <- lapply(files, function(x){
  ncol(read.table(x))-1 #subtract sample ID column
})

#read in files
range <- shapefile("Regulus_satrapa_22712594.shp")
state <- map_data("state")
map <- map_data("world")
range.df <- fortify(range)
names(loci) <- files
maf <- c(0,1,2,3,4,8,20)
for(i in 1:length(loci)){
  assign(paste("maf",maf[[i]],sep = ""), 
         read.structure(names(loci[i]), n.ind=33, n.loc=loci[[i]], 
                        onerowperind = F,col.lab=1,col.pop=0,
                        row.marknames=0,NA.char="-9",ask=F))
}

# scale data 
seq.scaled.0 <- scaleGen(maf0,NA.method="mean",scale=F)
seq.scaled.0125 <- scaleGen(maf0.0125,NA.method="mean",scale=F)
seq.scaled.025 <- scaleGen(maf0.025,NA.method="mean",scale=F)
seq.scaled.0375 <- scaleGen(maf0.0375,NA.method="mean",scale=F)
seq.scaled.05 <- scaleGen(maf0.05,NA.method="mean",scale=F)
seq.scaled.1 <- scaleGen(maf0.1,NA.method="mean",scale=F)
seq.scaled.25 <- scaleGen(maf0.25,NA.method="mean",scale=F)

# find clusters
k3.maf0 <- find.clusters(maf0, n.pca=95,n.clust = 3,choose.n.clust = F)
k3.maf0.012 <- find.clusters(maf0.0125, n.pca=95,n.clust = 3,choose.n.clust = F)
k3.maf0.025 <- find.clusters(maf0.025, n.pca=95,n.clust = 3,choose.n.clust = F)
k3.maf0.037 <- find.clusters(maf0.0375, n.pca=95,n.clust = 3,choose.n.clust = F)
k3.maf0.05 <- find.clusters(maf0.05, n.pca=95,n.clust = 3,choose.n.clust = F)
k3.maf0.1  <- find.clusters(maf0.1, n.pca=95,n.clust = 3,choose.n.clust = F)
k3.maf0.25 <- find.clusters(maf0.25, n.pca=95,n.clust = 3,choose.n.clust = F)

satr.west <- seppop(maf0,pop=rownames(maf0@tab) %in% subset(localities,region=="west")$sample.ID)[[2]] #[[2]] selects condition==T
satr.east <- seppop(maf0,pop=rownames(maf0@tab) %in% subset(localities,region=="east")$sample.ID)[[2]] 
satr.mex <- seppop(maf0,pop=rownames(maf0@tab) %in% subset(localities,region=="mexico")$sample.ID)[[2]]

#run PCA, merge output with specimen locations
pca <- prcomp(seq.scaled.0,center=F,scale=F)
screeplot(pca)
pc <- data.frame(pca$x[,1:3])
pc$sample.ID <- make.names(rownames(pc))
localities$sample.ID <- gsub(" ","",localities$sample.ID)
pc <- merge(pc,localities,by="sample.ID")

#plot PCA with sample IDs
ggplot(data=pc,aes(x=PC1,y=PC2,col=k3.maf0$grp))+geom_text(aes(label=sample.ID))

#PC1 by longitude
ggplot(data=pc,aes(x=long,y=PC1))+theme_minimal()+geom_text(aes(label=sample.ID))

#linear regression of PC1 against longitude
model <- lm(pc$PC1~pc$long)
summary(model)
ggplot(data=pc,aes(y=PC1,x=long))+theme_minimal()+geom_point()+stat_smooth(method="lm")

# check for correlation between missing data and PC
# a <- fread("./str_in/satrapa_unlinked.str",header=TRUE) #unlinked SNP's from same pyrad step 7 run as seq file. 
# nchar <- as.numeric(names(a)[2])
# colnames(a) <- c("sample","seq")
# b <- lapply(a$seq,FUN=function(e) str_count(e,"N"))
# c <- as.numeric(b)/nchar
# md <- data.frame(sample=a$sample,md=c) %>% cbind(pc=pc$PC1)
# lm(abs(md$pc)~md$md) %>% summary()
# plot(y=abs(md$pc),x=md$md)

#dapc
dapc <- dapc(maf0,pop=k3.maf0$grp,n.pca=1,n.da=1)
tmp <- data.frame(sample.ID=rownames(dapc$tab),dapc.coord=dapc$ind.coord,grp=dapc$grp)
pc <- merge(pc,tmp,by="sample.ID")

#cross validation across allele frequencies
xval0 <- xvalDapc(seq.scaled.0, k3.maf0$grp, n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval0125 <- xvalDapc(seq.scaled.0125, k3.maf0.012$grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval025 <- xvalDapc(seq.scaled.025, k3.maf0.025$grp, n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval0375 <- xvalDapc(seq.scaled.0375, k3.maf0.037$grp, n.pca.max = 300, training.set = 0.9,
                    result = "groupMean", center = TRUE, scale = FALSE,
                    n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval05 <- xvalDapc(seq.scaled.05, k3.maf0.05$grp, n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval1 <- xvalDapc(seq.scaled.1, k3.maf0.1$grp, n.pca.max = 300, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval25 <- xvalDapc(seq.scaled.25, k3.maf0.25$grp, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)

# nasty data manipulation because I need to spend more time in the tidyverse
maf0 <- xval0$`Cross-Validation Results`[,2]
maf0125 <- xval0125$`Cross-Validation Results`[,2]
maf025 <- xval025$`Cross-Validation Results`[,2]
maf037 <- xval0375$`Cross-Validation Results`[,2]
maf05 <- xval05$`Cross-Validation Results`[,2]
maf1 <- xval1$`Cross-Validation Results`[,2]
maf25 <- xval25$`Cross-Validation Results`[,2]
col <- rep("0",length(maf0))
maf0 <- cbind.data.frame(maf0,col)
colnames(maf0) <- c("success","MAF")
col1 <- rep("1",length(maf0125))
maf0125 <- cbind.data.frame(maf0125,col1)
colnames(maf0125) <- c("success","MAF")
col2 <- rep("2",length(maf025))
maf025 <- cbind.data.frame(maf025,col2)
colnames(maf025) <- c("success","MAF")
col3 <- rep("3",length(maf037))
maf037 <- cbind.data.frame(maf037,col3)
colnames(maf037) <- c("success","MAF")
col4 <- rep("4",length(maf05))
maf05 <- cbind.data.frame(maf05,col4)
colnames(maf05) <- c("success","MAF")
col5 <- rep("8",length(maf1))
maf1 <- cbind.data.frame(maf1,col5)
colnames(maf1) <- c("success","MAF")
col6 <- rep("20",length(maf25))
maf25 <- cbind.data.frame(maf25,col6)
colnames(maf25) <- c("success","MAF")
df <- rbind(maf0,maf0125,maf025,maf037,maf05,maf1,maf25)

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
          axis.title=element_text(size = 10)) +
    geom_point(size=1.5,shape=k3.maf0$grp,stroke=1,alpha=0.9) +
    annotate(geom = "text", y = max(pc$PC2), x = max(pc$PC1),label=text,hjust=1)
  print(p)
  plots <- c(plots, list(p))
}


# plot
boxplot <- ggplot(df,aes(x=MAF,y=success))+theme_minimal() +
  xlab("Minimum Minor Allele Count")+ylab("Population Discrimination")+
  coord_flip()+
  geom_boxplot()

joyplot <- ggplot(data=df,aes(x=success,y=MAF))+theme_minimal()+
  #facet_grid(~variable,scales="free")+
  ylab("Minimum Minor Allele Count")+xlab("Population Discrimination")+
  scale_y_discrete(expand = c(0.01, 0))+
  scale_x_discrete(expand = c(0.01, 0),limits=c(0.0,0.2,0.4,0.6,0.8,1))+
  geom_joy()

#side-by-side plots, 3 PCAs
layout <- cbind(c(1,1,1),c(1,1,1),c(2,3,4))
grid.arrange(joyplot,plots[[7]],plots[[4]],plots[[1]],layout_matrix=layout)

#layout <- cbind(c(1,1,1,1,1,1,1,1,1),c(1,1,1,1,1,1,1,1,1),c(1,1,1,1,1,1,1,1,1),c(1,1,1,1,1,1,1,1,1),c(NA,2,3,4,5,6,7,8,NA))
#grid.arrange(joyplot,plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],layout_matrix=layout)

# run lm 
model <- lm(df$success~df$long)
summary(model)
ggplot(data=df,aes(y=success,x=MAF))+theme_minimal()+geom_point()+stat_smooth(method="lm")

# plot localities
range<-readOGR("Regulus_satrapa_22712594.shp")
range.crop <- crop(range,c(-135,-50,13,55))
range.df <- tidy(range.crop)
state <- map_data("state")
map <- map_data("world")
range.df <- subset(range.df,id != 1)
range.df$season <- gsub("0","Year Round",range.df$id) %>% gsub("2","Summer",.) %>% gsub("3","Winter",.)
localities$region <- as.factor(localities$region)
levels(localities$region) <- c("East","Mexico","West")
rmap <- ggplot()+coord_map()+theme_bw()+ylim(13,55)+xlim(-135,-50)+
  scale_size(breaks=c(1,2,3,4),name="Samples",range=c(1.5,6))+
  scale_color_grey() +
  scale_shape_discrete(solid = T, name = "Cluster")+
  scale_fill_manual(name="Season",values = c("grey55","grey85","grey40"))+
  theme(panel.grid=element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
        axis.text = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  geom_polygon(data=range.df,
               aes(x=long,y=lat,group=group,fill=range.df$season),col=NA)+
  geom_path(data=state,aes(x=long,y=lat,group=group),col="grey",lwd=0.25)+
  geom_path(data=map,aes(x=long,y=lat,group=group))+
  geom_point(data=localities, aes(x=long, y=lat, size=count, shape=region),stroke=1,alpha=0.9)

#plot w/ range map
layout <- cbind(c(1,1,1),c(1,1,1),c(2,3,4))
layout <- rbind(layout, c(5,5,5),c(5,5,5),c(5,5,5))
grid.arrange(joyplot,plots[[7]],plots[[4]],plots[[1]],rmap,layout_matrix=layout)

### simulations!

# list files
simfiles <- list.files("~/Dropbox/structure_simulations/str_in", pattern = "sim01.*", full.names = TRUE)
simfiles <- simfiles[-c(7)]

# count loci
simloci <- lapply(simfiles, function(x){
  ncol(read.table(x))-2 #subtract sample ID column + pop ID
})

# match empirical
simloci <- simloci[-c(7)]


names(simloci) <- simfiles
maf <- c(0,1,2,20,3,4,8)
for(i in 1:length(simloci)){
  assign(paste("sim_maf",maf[[i]],sep = ""), 
         read.structure(names(simloci[i]), n.ind=30, n.loc=simloci[[i]], 
                        onerowperind = F,col.lab=1,col.pop=1,
                        row.marknames=0,NA.char="-9",ask=F))
}

# scale data 
seq.scaled.sim0 <- scaleGen(sim_maf0,NA.method="mean",scale=F)
seq.scaled.sim1 <- scaleGen(sim_maf1,NA.method="mean",scale=F)
seq.scaled.sim2 <- scaleGen(sim_maf2,NA.method="mean",scale=F)
seq.scaled.sim20 <-scaleGen(sim_maf20,NA.method="mean",scale=F)
seq.scaled.sim3 <- scaleGen(sim_maf3,NA.method="mean",scale=F)
seq.scaled.sim4 <- scaleGen(sim_maf4,NA.method="mean",scale=F)
seq.scaled.sim8 <- scaleGen(sim_maf8,NA.method="mean",scale=F)

# find clusters
sim.k3.0  <- find.clusters(sim_maf0,n.pca=95,n.clust = 3,choose.n.clust = F)
sim.k3.1  <- find.clusters(sim_maf1,n.pca=95,n.clust = 3,choose.n.clust = F)
sim.k3.2  <- find.clusters(sim_maf2,n.pca=95,n.clust = 3,choose.n.clust = F)
sim.k3.20 <- find.clusters(sim_maf20, n.pca=95,n.clust = 3,choose.n.clust = F)
sim.k3.3  <- find.clusters(sim_maf3,n.pca=95,n.clust = 3,choose.n.clust = F)
sim.k3.4  <- find.clusters(sim_maf4,n.pca=95,n.clust = 3,choose.n.clust = F)
sim.k3.8  <- find.clusters(sim_maf8,n.pca=95,n.clust = 3,choose.n.clust = F)


pop1 <- c("1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8","1_9","1_10")
pop2 <- c("2_1","2_2","2_3","2_4","2_5","2_6","2_7","2_8","2_9","2_10")
pop3 <- c("3_1","3_2","3_3","3_4","3_5","3_6","3_7","3_8","3_9","3_10")
sim_1 <- seppop(sim_maf0,pop=rownames(sim_maf0@tab) %in% pop1)[[2]] #[[2]] selects condition==T
sim_2 <- seppop(sim_maf0,pop=rownames(sim_maf0@tab) %in% pop2)[[2]] 
sim_3 <- seppop(sim_maf0,pop=rownames(sim_maf0@tab) %in% pop3)[[2]] 

#dapc
dapc <- dapc(sim_maf0,pop=sim.k3.0$grp,n.pca=1,n.da=1)
tmp <- data.frame(sample.ID=rownames(dapc$tab),dapc.coord=dapc$ind.coord,grp=dapc$grp)


#cross validation across allele frequencies
simxval0 <- xvalDapc(seq.scaled.sim0, sim.k3.0$grp, n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
simxval1 <- xvalDapc(seq.scaled.sim1, sim.k3.1$grp, n.pca.max = 300, training.set = 0.9,
                  result = "groupMean", center = TRUE, scale = FALSE,
                  n.pca = NULL, n.rep = 30, xval.plot = TRUE)
simxval2 <- xvalDapc(seq.scaled.sim2, sim.k3.2$grp, n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
simxval20 <- xvalDapc(seq.scaled.sim20, sim.k3.20$grp, n.pca.max = 300, training.set = 0.9,
                    result = "groupMean", center = TRUE, scale = FALSE,
                    n.pca = NULL, n.rep = 30, xval.plot = TRUE)
simxval3 <- xvalDapc(seq.scaled.sim3, sim.k3.3$grp, n.pca.max = 300, training.set = 0.9,
                     result = "groupMean", center = TRUE, scale = FALSE,
                     n.pca = NULL, n.rep = 30, xval.plot = TRUE)
simxval4 <- xvalDapc(seq.scaled.sim4, sim.k3.4$grp, n.pca.max = 300, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = NULL, n.rep = 30, xval.plot = TRUE)
simxval8 <- xvalDapc(seq.scaled.sim8, sim.k3.8$grp, n.pca.max = 300, training.set = 0.9,
                   result = "groupMean", center = TRUE, scale = FALSE,
                   n.pca = NULL, n.rep = 30, xval.plot = TRUE)

# nasty data manipulation because I need to spend more time in the tidyverse
simr0 <- simxval0$`Cross-Validation Results`[,2]
simr1 <- simxval1$`Cross-Validation Results`[,2]
simr2 <- simxval2$`Cross-Validation Results`[,2]
simr20 <- simxval20$`Cross-Validation Results`[,2]
simr3 <- simxval3$`Cross-Validation Results`[,2]
simr4 <- simxval4$`Cross-Validation Results`[,2]
simr8 <- simxval8$`Cross-Validation Results`[,2]

simcol <- rep("0",length(simr0))
simr0 <- cbind.data.frame(simr0,simcol)
colnames(simr0) <- c("success","MAF")

simcol1 <- rep("1",length(simr1))
simr1 <- cbind.data.frame(simr1,simcol1)
colnames(simr1) <- c("success","MAF")

simcol2 <- rep("2",length(simr2))
simr2 <- cbind.data.frame(simr2,simcol2)
colnames(simr2) <- c("success","MAF")

simcol20 <- rep("20",length(simr20))
simr20 <- cbind.data.frame(simr20,simcol20)
colnames(simr20) <- c("success","MAF")

simcol3 <- rep("3",length(simr3))
simr3 <- cbind.data.frame(simr3,simcol3)
colnames(simr3) <- c("success","MAF")

simcol4 <- rep("4",length(simr4))
simr4 <- cbind.data.frame(simr4,simcol4)
colnames(simr4) <- c("success","MAF")

simcol8 <- rep("8",length(simr8))
simr8 <- cbind.data.frame(simr8,simcol8)
colnames(simr8) <- c("success","MAF")

simdf <- rbind(simr0,simr1,simr2,simr3,simr4,simr8,simr20)

# read each file, run & plot PCA
simplots <- list()
simmaf <- c(0,1,2,20,3,4,8)
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
          axis.title=element_text(size = 10)) +
    geom_point(size=1.5,shape=sim.k3.0$grp,stroke=1,alpha=0.9) +
    annotate(geom = "text", y = max(pc$PC2), x = max(pc$PC1),label=sim.text,hjust=1)
  print(x)
  simplots <- c(simplots, list(x))
}

simjoyplot <- ggplot(data=simdf,aes(x=success,y=MAF))+theme_minimal()+
  #facet_grid(~variable,scales="free")+
  ylab("Minimum Minor Allele Count")+xlab("Population Discrimination")+
  scale_y_discrete(expand = c(0.01, 0))+
  scale_x_discrete(expand = c(0.01, 0),limits=c(0.0,0.2,0.4,0.6,0.8,1))+
  geom_joy()


#side-by-side plots, 3 PCAs
layout <- cbind(c(1,1,1),c(1,1,1),c(2,3,4))
r2 <- cbind(c(5,5,5),c(5,5,5),c(6,7,8))
layout <- rbind(layout,r2)
grid.arrange(joyplot,plots[[7]],plots[[5]],plots[[1]],simjoyplot,simplots[[4]],simplots[[6]],simplots[[1]],layout_matrix=layout)




