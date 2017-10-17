### summary stats for results section

library(ggplot2);library(reshape)
source("~/Dropbox/structure_simulations/R_scripts/ggthemes.R")

# read in file
infile <- readLines("~/Downloads/satrapa_stats.txt",warn=FALSE)
a <- grep("Final Sample stats summary",infile)
stats <- read.table("~/Downloads/satrapa_stats.txt",skip=a-1)

# stats 
mean_filtered_reads <- mean(stats$reads_passed_filter)
mean_clusters <- mean(stats$clusters_total)
avg_depth <- mean_filtered_reads/mean_clusters
mean_assembly <- mean(stats$loci_in_assembly)

# calculate mean # loci per filtered sim file
simfiles <- list.files("~/Dropbox/structure_simulations/str_in/k3", full.names = TRUE)
simloci <- lapply(simfiles, function(x){
  ncol(read.table(x))-2
})
names(simloci) <- simfiles
maf <- list("maf1","maf2","maf3","maf4","maf5","maf8","maf20")
mean <- lapply(maf, function(i){
  a <- grepl(i, names(simloci))
  b <- simloci[a]
  mean <- mean(as.numeric(b))
})
names(mean) <- maf

# for empirical data
files <- list.files("~/Dropbox/structure_simulations/satrapa/str_in/", pattern = "*.str", full.names = TRUE)
files <- files[1:7]
loci <- lapply(files, function(x){
  ncol(read.table(x))-1
})
names(loci) <- files

#subset
loci <- loci[c(1:7)]

#plot SFS
# simulated data
b <- read.table("~/Dropbox/structure_simulations/str_in/k3/sim01.str")
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
mafsim[1,2]/ncol(b) #% singletons: 0.4459161
mafsim[2,2]/ncol(b) #% doubletons: 0.1264901
mafsim[3,2]/ncol(b) #% tripletons: 0.06843267
mafsim[4,2]/ncol(b) #% 0.0410596
mafsim[5,2]/ncol(b) #% 0.02626932

mafsim <- cbind(mafsim, c(1:nrow(mafsim)))
colnames(mafsim) <- c("maf","count","minor_allele_frequency")

# empirical data
c <- read.table("~/Dropbox/structure_simulations/satrapa/str_in/satrapa_unlinked.str")
c <- c[,-c(1:2)]
mafreg <- list()
mafreg <- apply(c, 2, function(y){
  df <- y
  df[df==-9] <- NA
  df <- as.data.frame(table(df))
  df <- as.vector(df$Freq)
  if (length(df)==1){
    ma <- 0
  } else {
    ma <- min(df)
  }
  ma/nrow(c)
})

# summary stats
min(mafreg)
max(mafreg)
mean(mafreg)
unique(mafreg)

# clean data
mafreg <- as.vector(mafreg)
mafreg <- as.data.frame(table(mafreg))
mafreg[1,2]/ncol(c) #% singletons = 0.3666923
mafreg[2,2]/ncol(c) #% doubletons = 0.1696177
mafreg[3,2]/ncol(c) #% tripletons = 0.0805748
mafreg[4,2]/ncol(c) #% 0.05542725
mafreg[5,2]/ncol(c) #% 0.03900436


mafreg <- mafreg[-c(1),]
mafreg <- cbind(mafreg, c(1:nrow(mafreg)))
colnames(mafreg) <- c("maf","count","minor_allele_frequency")

MAF <- c(1:20)
df <- as.data.frame(cbind(MAF, mafsim$count[1:20],mafreg$count[1:20]))
colnames(df) <- c("minor_allele_frequency","Simulated","Empirical")
df.m <- melt(df, id.vars='minor_allele_frequency')

#plot paired barcharts
ggplot(df.m, aes(minor_allele_frequency, value))+
  geom_bar(aes(fill=variable), position="dodge",stat="identity")+
  xlab("Minor allele count")+
  ylab("Count")+
  scale_fill_grey()+
  theme(panel.border = element_rect(color="black",fill=NA),
        panel.background = element_blank(),
        axis.title=element_text(size=10))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.text=element_text(size=10)) +
  theme(legend.position = c(0.9, 0.9))+
  scale_x_continuous(breaks = c(1:20,1), limits = c(0,21))

