#satrapa analysis pipeline
library(magrittr);library(data.table);library(foreach);library(doMC);
library(ggplot2);library(plyr);library(ggridges);library(reshape)
setwd("~/Dropbox/structure_simulations/satrapa/")
registerDoMC(cores=8)
source("../R_scripts/maf/maf_functions.R")

#read in empirical structure files and filter by maf
filter_by_mac(infile="str_in/satrapa_unlinked.str",pop.info = F,mac=c(1,2,3,4,5,8,11,15))

#run structure (copy to wopr & run in a screen)
setwd("/media/burke/bigMac/Dropbox/structure_simulations/satrapa/")
library(foreach);library(doMC);library(data.table);library(magrittr)
registerDoMC(cores=30)

reps <- 10 #number of independent analyses per input file
structure_path <- "/media/burke/bigMac/Dropbox/structure_kernel_src/structure"   #"/Applications/structure/structure"
mainparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/satrapa/str_params/mainparams.txt"
extraparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/satrapa/str_params/extraparams.txt"
params_dir <- "/media/burke/bigMac/Dropbox/structure_simulations/satrapa/str_params/"
str_in <- "/media/burke/bigMac/Dropbox/structure_simulations/satrapa/str_in/"
str_out <- "/media/burke/bigMac/Dropbox/structure_simulations/satrapa/str_out/"

files <- list.files(str_in) %>% grep("mac",.,value=T)
n_loci <- c()
for(i in files){
  tmp <- fread(paste0(str_in,i)) #use fread() for significant speed increase if str files don't have extra whitespace columns (otherwise read.table())
  n_loci <- append(n_loci,ncol(tmp)-1) #-1 if no pop info, -2 if pop info present.
}
names(n_loci) <- files

structure_commands <- c()
og_params <- readLines(mainparams_path)
ex_params <- readLines(extraparams_path)
for(i in files){
  for(j in 1:reps){
    og_params[grep("INFILE",og_params)] <- paste0("#define INFILE  ",str_in,i) #edit infile line 
    og_params[grep("NUMLOCI",og_params)] <- paste0("#define NUMLOCI    ",n_loci[i]) # edit n loci
    og_params[grep("OUTFILE",og_params)] <- paste0("#define OUTFILE ",str_out,i,"_",
                                                   formatC(j,digits=1,flag="0",format = "d")) #edit outfile line
    ex_params[grep("SEED",ex_params)] <- paste0("#define SEED     ",sample(1:10000,1))
    writeLines(og_params,paste0(params_dir,"mainparams_",i,"_",
                                formatC(j,digits=1,flag="0",format = "d"),".txt"))
    writeLines(ex_params,paste0(params_dir,"extraparams_",i,"_",
                                formatC(j,digits=1,flag="0",format = "d"),".txt"))
    structure_commands <- append(structure_commands,paste0(structure_path,
                                                           " -m ",paste0(params_dir,"mainparams_",i,"_",
                                                                         formatC(j,digits=1,flag="0",format = "d"),".txt"),
                                                           " -e ",paste0(params_dir,"extraparams_",i,"_",
                                                                         formatC(j,digits=1,flag="0",format = "d"),".txt")))
  }
}

foreach(i=structure_commands) %dopar% system(i)

#run multivariate clustering
files <- list.files("str_in",full.names = T) %>% grep("mac",.,value = T)
clust <- foreach(i=files,.combine = rbind) %dopar% cluster_multivar(i,pop.info=F,nreps=10,
                                                                    pop=c(2,2,2,2,2,2,2,2,2,
                                                                         2,2,2,2,2,2,2,2,3,
                                                                         3,3,3,1,1,1,1,2,2,
                                                                         2,2,2,1,1,1))

#summarize structure output
files <- list.files("str_out",full.names=T) %>% grep("mac",.,value=T)
q_matrices <- list()
sum_stats <- list()
j <- 1
for(i in files){
  tmp <- readLines(i,warn=F)
  mac <- i %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit("\\.") %>% unlist() %>% .[1] %>% as.numeric()
  run <- i %>% strsplit("_") %>% unlist() %>% .[5] %>% as.numeric()
  alpha <- tmp[grep("Mean value of alpha",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
  lnL_mean <- tmp[grep("Mean value of ln likelihood",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
  lnL <- tmp[grep("Estimated Ln Prob of Data",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
  
  q <- tmp[(grep("Inferred ancestry of individuals:",tmp)+2):(grep("Inferred ancestry of individuals:",tmp)+34)] %>% 
    lapply(function(e){ strsplit(e," * ") %>% unlist() %>% .[-c(1,2,4,5,9,10,11)]}) %>% 
    do.call(rbind.data.frame,.)
  colnames(q) <- c("id","1","2","3")
  q$pop <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,2,1,1,1)
  q[,2:4] <- q[,2:4] %>% apply(2,function(e) as.numeric(as.character(e)))
  
  #fraction individuals assigned to unique and mutually exclusive populations based on majority ancestry cluster
  ind_assignments <- apply(q[,2:4],1,function(z) as.numeric(names(z[which(z==max(z))])[1]))
  popclust <- c(Mode(ind_assignments[c(1:17,26:30)]),Mode(ind_assignments[c(18:21,31:33)]),Mode(ind_assignments[22:25]))
  missingclust <- c(1:3)[!c(1:3) %in% popclust]
  if(length(unique(popclust))==1){
    accuracy <- 0
  } else if (length(unique(popclust))==2){
    pop1 <- sum(ind_assignments[c(1:17,26:30)]==popclust[1])
    pop2 <- sum(ind_assignments[c(18:21,31:33)]==popclust[2])
    pop3 <- sum(ind_assignments[22:25]==missingclust)
    accuracy <-(pop1+pop2+pop3)/33
  } else if (length(unique(popclust))==3){
    pop1 <- sum(ind_assignments[c(1:17,26:30)]==popclust[1])
    pop2 <- sum(ind_assignments[c(18:21,31:33)]==popclust[2])
    pop3 <- sum(ind_assignments[22:25]==popclust[3])
    accuracy <-(pop1+pop2+pop3)/33
  }
  
  #swap column names to minimize label switching in structure plots
  clustnames <- c("q1","q2","q3")
  newclustnames <- c(rep(NA,nlevels(factor(q$pop))))
  for(i in 1:3){
    d <- subset(q,pop==i)
    e <- colMeans(d[2:4])
    f <- as.numeric(names(e[which(e==max(e))]))
    if(is.na(newclustnames[f])){
      newclustnames[f] <- clustnames[i]
    }
  }
  newclustnames[which(is.na(newclustnames))] <- clustnames[which(clustnames %in% newclustnames==F)]
  colnames(q) <- c("id",newclustnames,"pop")
  
  q$mac <- mac
  q$run <- run
  q$lnL <- lnL
  q$lnL_mean <- lnL_mean
  
  q_matrices[[j]] <- q
  sum_stats[[j]] <- c(mac,run,alpha,lnL,lnL_mean,accuracy)
  j <- j+1
}
sum_stats <- do.call(rbind.data.frame,sum_stats)
names(sum_stats) <- c("mac","run","alpha","lnL","lnL_mean","accuracy")
sum_stats$log_alpha <- log(sum_stats$alpha)

mean_q_list <- lapply(q_matrices,function(e) {
  ddply(e,"pop",summarize,q1=mean(q1),q2=mean(q2),q3=mean(q3))[-1]
})
q_dist_list <- sapply(mean_q_list,function(e){
  mean(dist(e))
})
sum_stats$`Population Discrimination` <- q_dist_list/1.414214


sum_stats_wide <- sum_stats
sum_stats2 <- melt(sum_stats,id.vars = c("mac","run"))

#ridge plot
ggplot(data=subset(sum_stats2,variable %in% c("alpha","Population Discrimination")),
       aes(y=factor(mac),x=value))+
  theme_minimal()+theme(strip.background = element_blank())+
  facet_wrap(~variable,scales="free")+
  scale_y_discrete(expand=c(0.01,0))+
  ylab("Minimum Minor Allele Count")+xlab("")+
  geom_density_ridges(alpha=0.5,rel_min_height=0)

#structure plots (highest lnL per mac for each simulation)
strplotdata <- do.call(rbind,q_matrices)
strplotdata$q_dist <- unlist(lapply(q_dist_list,function(e) unlist(rep(e,33))))
best_runs <- ddply(strplotdata,.(mac,sim),summarize,max_lnL=max(lnL),max_q_dist=max(q_dist),mean_lnL=max(lnL_mean))
#strplotdata <- subset(strplotdata,q_dist %in% best_runs$max_q_dist)
strplotdata <- subset(strplotdata,lnL %in% best_runs$max_lnL)
meltq <- melt(strplotdata[-c(8,9,10)],id.vars=c("id","pop","mac","run"))
ggplot(data=meltq,aes(x=id,y=value,fill=variable))+
  facet_grid(mac~.)+
  theme_minimal()+theme(axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        strip.background = element_blank(),
                        axis.text.x=element_blank(),
                        rect = element_blank())+
  ylab("mac")+xlab("")+
  scale_fill_manual(values = grey.colors(3)[c(2,1,3)])+
  geom_bar(stat="identity",width=.9)

