#PCA MAF threshold comparisons

#install.packages
#load libraries
library(adegenet);library(plyr);library(ggplot2);library(viridis);library(grid);library(magrittr);library(raster)
library(foreach);library(stringr);library(data.table);library(gridExtra)
#library(doMC)
setwd("~/Dropbox/structure_simulations/satrapa/")

#colors
satrapa.palette <- c("#FFAE00","#3E9C3E","#000000")

#drop NA column
files <- list.files("./str_in/", full.names = TRUE)
filelist <- lapply(files, read.table)
#filelist <- lapply(filelist, function(i){
#  i[,-c(2)]
#})

# append original file names to list
names(filelist) <- files

# rewwrite files 
# sapply(names(filelist), function(x){
#  write.table(filelist[[x]], file=paste(x), col.names = F, row.names = F, quote = F)
#})

# calculate missing data per row for multiple files
# md <- list()
# md <- lapply(files, function(x){
#  a <- read.table(x)
#  b <- ncol(a)-1
#  d <- apply(a,1,function(y){
#    a <- a[2:nrow(a)]
#    c <- length(which(y==-9))
#    md[x] <- c/b
#  })
#})

# calculate missing data from full alignment
a <- read.table("./str_in/satrapa_unlinked.str")
b <- ncol(a)-1
md <- vector()
md <- apply(a,1,function(y){
  a <- a[2:nrow(a)]
  c <- length(which(y==-9))
  md <- c/b
})
taxa_md <- cbind(as.data.frame(a$V1), md)
taxa_md <- unique(taxa_md)
taxa_md <- taxa_md[order(taxa_md$md, decreasing = TRUE),]
colnames(taxa_md) <- c("sample", "missing_data")

#plot md by taxa
ggplot(taxa_md, aes(x = reorder(sample, -missing_data), y = missing_data)) +
  xlab("sample ID")+
  ylab("proportion missing data")+
  geom_bar(stat = "identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

# remove 7 individuals
#r <- as.factor(c("sat_AZ_56","sat_WA_16","sat_WA_15","sat_CO_14","sat_AZ_64","sat_MX_71,""sat_WV_4"))
#filelist <- lapply(filelist, function(x){
#    x[!(x$V1 %in% r),]
#})

# rewrite files 
#sapply(names(filelist), function(x){
#  write.table(filelist[[x]], file=paste(x), col.names = F, row.names = F, quote = F)
#})

# count loci
files <- list.files("./str_in", pattern = "*.str", full.names = TRUE)
files <- files[1:7]
loci <- lapply(files, function(x){
  ncol(read.table(x))-1
})
names(loci) <- files

#subset
loci <- loci[c(1:7)]

# read each file, run & plot PCA
plots <- list()
maf <- c(0,1/80,2/80,3/80,4/80,.1,.25)
for(i in 1:length(loci)){
  seq <- read.structure(names(loci[i]), n.ind=33, n.loc=loci[[i]], onerowperind = F,col.lab=1,col.pop=0,row.marknames=0,NA.char="-9",ask=F)
  text <- paste("MAF=",maf[i],sep="")
  seq.scaled <- scaleGen(seq,NA.method="mean",scale=F)
  clust.opt <- find.clusters(seq,n.pca=95,n.clust=3,choose.n.clust = F)
  clust <- cbind(sampleID=rownames(seq@tab),clust.opt=unname(clust.opt$grp)) %>% data.frame()
  pca <- prcomp(seq.scaled,center=F,scale=F)
  pc <- data.frame(pca$x[,1:3])
  p <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.opt$grp))+xlab("PC1")+ylab("PC2")+
    scale_color_manual(values = satrapa.palette) +
    theme(panel.border = element_rect(color="black",fill=NA),
          legend.position="bottom",
          legend.title = element_blank(),
          panel.background = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_text(size=10))+
    geom_point(size=1.5) +
    annotate(geom = "text", y = max(pc$PC2), x = max(pc$PC1),label=text,hjust=1)
  print(p)
  plots <- c(plots, list(p))
}
print(clust.opt$grp)
levels(clust.opt$grp) <- c("Eastern","Western","Mexico")


# calculate maf
a <- read.table("./str_in/satrapa_unlinked.str")
a <- a[,-1]
maf <- list()
maf <- apply(a, 2, function(x){
  df <- x
  df[df==-9] <- NA
  df <- as.data.frame(table(df))
  df <- as.vector(df$Freq)
  if (length(df)==1){
    ma <- 0
  } else {
    ma <- min(df)
  }
  ma/nrow(a)
})

# summary stats
min(maf)
max(maf)
mean(maf)
unique(maf)

# clean data
maf <- as.vector(maf)
maf <- as.data.frame(table(maf))
maf <- maf[-c(1),]
maf <- cbind(maf, c(1:nrow(maf)))
colnames(maf) <- c("maf","count","minor_allele_frequency")

# plot
ggplot(maf, aes(x = minor_allele_frequency, y = count)) +
  xlab("Minor allele frequncy")+
  ylab("Count")+
  geom_bar(stat = "identity") +
  theme(panel.border = element_rect(color="black",fill=NA),
        panel.background = element_blank(),
        axis.title=element_text(size=10))+
  scale_x_continuous(breaks = c(1:32,1), limits = c(0,32))

### plot simulated results

# count loci
simfiles <- list.files("~/Dropbox/structure_simulations/str_in", pattern = "sim01", full.names = TRUE)
simloci <- lapply(simfiles, function(x){
  ncol(read.table(x))-2
})
names(simloci) <- simfiles

#plots <- list()
#maf <- c(1/80,2/80,3/80,4/80,.1,.25,0)
#for(i in 1:length(simloci)){
#  seq <- read.structure(names(simloci[i]), n.ind=33 ,n.loc=simloci[[i]], onerowperind = F,col.lab=1,col.pop=0,row.marknames=0,NA.char="-9",ask=F)
#  text <- paste("MAF=",maf[i],sep="")
#  seq.scaled <- scaleGen(seq,NA.method="mean",scale=F)
#  clust.opt <- find.clusters(seq,n.pca=95,choose.n.clust = F)
#  clust <- cbind(sampleID=rownames(seq@tab),clust.opt=unname(clust.opt$grp)) %>% data.frame()
#  pca <- prcomp(seq.scaled,center=F,scale=F)
#  pc <- data.frame(pca$x[,1:3])
#  p <- ggplot(pc,aes(x=PC1,y=PC2,col=clust.opt$grp))+xlab("PC1")+ylab("PC2")+
#    scale_color_manual(values = satrapa.palette) +
#    theme(panel.border = element_rect(color="black",fill=NA),
#          legend.position="bottom",
#          legend.title = element_blank(),
#          panel.background = element_blank(),
#          axis.text=element_blank(),
#          axis.ticks=element_blank(),
#          axis.title=element_text(size=10))+
#    geom_point(size=1.5) +
#    annotate(geom = "text", y = max(pc$PC2), x = max(pc$PC1),label=text,hjust=1)
#  print(p)
#  plots <- c(plots, list(p))
#}

# read in true pops
pops <- vector(mode = "list", length = length(simloci))
for(i in 1:length(simloci)){
  a <- read.table(names(simloci[i]))
  pops[[i]] <- as.vector(a[,2])
}
names(pops) <- names(simloci)

# all pop assignments true
infpops <- vector(mode = "list", length = length(simloci))
for(i in 1:length(simloci)){
  seq <- read.structure(names(simloci[i]), n.ind=40 ,n.loc=simloci[[i]], onerowperind = F,col.lab=1,col.pop=0,row.marknames=0,NA.char="-9",ask=F)
  seq.scaled <- scaleGen(seq,NA.method="mean",scale=F)
  clust.opt <- find.clusters(seq,n.pca=95,choose.n.clust = F)
  print(clust.opt$grp)
  b <- as.numeric(clust.opt$grp)
  infpops[[i]] <- b
}


# calculate maf
b <- read.table("~/Dropbox/structure_simulations/str_in/sim01.str")
b <- b[,-c(1:2)]
mafsim <- list()
mafsim <- apply(b, 2, function(y){
  df <- y
  df[df==-9] <- NA
  df <- as.data.frame(table(df))
  df <- as.vector(df$Freq)
  if (length(df)==1){
    ma <- 0
  } else {
    ma <- min(df)
  }
  ma/nrow(b)
})

# summary stats
min(mafsim)
max(mafsim)
mean(mafsim)
unique(mafsim)

# clean data
mafsim <- as.vector(mafsim)
mafsim <- as.data.frame(table(mafsim))
mafsim <- mafsim[-c(1),]
mafsim <- cbind(mafsim, c(1:nrow(mafsim)))
colnames(mafsim) <- c("maf","count","minor_allele_frequency")

# plot
ggplot(mafsim, aes(x = minor_allele_frequency, y = count)) +
  xlab("Minor allele frequncy")+
  ylab("Count")+
  geom_bar(stat = "identity") +
  theme(panel.border = element_rect(color="black",fill=NA),
        panel.background = element_blank(),
        axis.title=element_text(size=10))+
  scale_x_continuous(breaks = c(1:40,1), limits = c(0,40))

# plot all on grid
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(plots[[7]],plots[[4]],plots[[1]])
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

grid_arrange_shared_legend(plots[[7]],,plots[[4]],plots[[1]],
                          ncol=1,nrow=3,position = c("bottom"))

### fst

