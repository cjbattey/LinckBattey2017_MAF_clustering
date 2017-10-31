#3-population simulation and structure inference pipeline
library(magrittr);library(data.table);library(foreach);library(doMC);
library(ggplot2);library(plyr);library(ggridges);library(reshape)
setwd("~/Dropbox/structure_simulations/")
registerDoMC(cores=8)
source("./R_scripts/maf/maf_functions.R")

####################################
######## simulate data #############

#run fsc simulation using parameters in ./fsc_params/sim.par
system("cd ~/Dropbox/structure_simulations/;
       /Applications/fsc_mac64/fsc25221 -i ./fsc_params/sim.par -g -n 10") #edit sim.par (or new par) as needed.

#convert fastsimcoal (arlequin) to structure format
files <- list.files("sim",full.names = T)
files <- files[grep(".arp",files)]
foreach(i=files) %dopar% arp2structure(infile=i,md=NULL,npops=3,samples_per_pop=10,out_directory="str_in/")

#Filter sites by minimum minor allele count
files <- list.files("./str_in",full.names = T)
foreach(i=files) %dopar% filter_by_mac(infile = i,mac=c(2,3,4,5,8,11,15),pop.info=T)


#########################################
### run structure in parallel on wopr ###

setwd("/media/burke/bigMac/Dropbox/structure_simulations/")
library(foreach);library(doMC);library(data.table)
registerDoMC(cores=30)

#set options and file paths here (use full paths)
reps <- 10 #number of independent analyses per input file
structure_path <- "/media/burke/bigMac/Dropbox/structure_kernel_src/structure"   #"/Applications/structure/structure"
mainparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/str_params/mainparams.txt"
extraparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/str_params/extraparams.txt"
params_dir <- "/media/burke/bigMac/Dropbox/structure_simulations/str_params/"
str_in <- "/media/burke/bigMac/Dropbox/structure_simulations/str_in/"
str_out <- "/media/burke/bigMac/Dropbox/structure_simulations/str_out/"

#scan files in str_in for n loci
files <- list.files(str_in)
n_loci <- c()
for(i in files){
  tmp <- fread(paste0(str_in,i)) #use fread() for significant speed increase if str files don't have extra whitespace columns (otherwise read.table())
  n_loci <- append(n_loci,ncol(tmp)-2) #-1 if no pop info, -2 if pop info present.
}
names(n_loci) <- files

#write new params files and paste together structure commands
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

#run in parallel
foreach(i=structure_commands) %dopar% system(i)


###########################################
######## multivariate clustering ##########

files <- list.files("str_in",full.names = T)
clust <- foreach(i=files,.combine = rbind) %dopar% cluster_multivar(i,pop.info=T,nreps=10)
write.table(clust,"multivariate_clusters.csv",sep=",",row.names = F,col.names = T,quote = F)

clust <- read.csv("multivariate_clusters.csv")
clust$sim <- sapply(clust$file,function(e) {as.character(e) %>% #moved this to the clust_multivar() function. 
                                              strsplit("sim") %>% 
                                                unlist() %>% 
                                                  .[2] %>% 
                                                    strsplit("\\.") %>% 
                                                      unlist() %>%
                                                        .[1] %>%
                                                          as.numeric()}) #sorry
clust$mac <- sapply(clust$file,function(e) {if(!grepl("mac",e)){1} else{as.character(e) %>% 
                                                                          strsplit("mac") %>% 
                                                                            unlist() %>% 
                                                                              .[2] %>% 
                                                                                strsplit("\\.") %>% 
                                                                                  unlist() %>%
                                                                                    .[1] %>%
                                                                                      as.numeric()}})
melt_clust <- melt(clust[-3],id.vars = c("sim","mac"))
ggplot(data=melt_clust,aes(x=value,y=factor(mac)))+
  theme_minimal()+
  theme(strip.background = element_blank())+
  facet_wrap(~variable)+
  xlim(0,1)+
  xlab("Assignment accuracy")+ylab("Minimum minor allele count")+
  scale_y_discrete(expand = c(0.01,0))+
  geom_density_ridges()
                                    
                                                  

#########################################
###### summarize structure output #######

files <- list.files("./str_out",full.names=T)
q_matrices <- list()
sum_stats <- list()
j <- 1
for(i in files){
  tmp <- readLines(i,warn=F)
  if(!grepl("mac",i)){
    mac <- 1
    run <- i %>% strsplit("_") %>% unlist() %>% .[3] %>% as.numeric()
  } else{
    mac <- i %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit("\\.") %>% unlist() %>% .[1] %>% as.numeric()
    run <- i %>% strsplit("_") %>% unlist() %>% .[4] %>% as.numeric()
  }
  sim <- i %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("\\.") %>% unlist() %>% .[1] %>% as.numeric()
  alpha <- tmp[grep("Mean value of alpha",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
  lnL_mean <- tmp[grep("Mean value of ln likelihood",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
  lnL <- tmp[grep("Estimated Ln Prob of Data",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
  
  q <- tmp[(grep("Inferred ancestry of individuals:",tmp)+2):(grep("Inferred ancestry of individuals:",tmp)+31)] %>% 
        lapply(function(e){ strsplit(e," * ") %>% unlist() %>% .[-c(1,2,4,6,10,11,12)]}) %>% 
          do.call(rbind.data.frame,.)
  colnames(q) <- c("id","pop","1","2","3")
  q[,3:5] <- q[,3:5] %>% apply(2,function(e) as.numeric(as.character(e)))

  #fraction individuals assigned to unique and mutually exclusive populations based on majority ancestry cluster
  ind_assignments <- apply(q[,3:5],1,function(z) as.numeric(names(z[which(z==max(z))])[1]))
  popclust <- c(Mode(ind_assignments[1:10]),Mode(ind_assignments[11:20]),Mode(ind_assignments[21:30]))
  missingclust <- c(1:3)[!c(1:3) %in% popclust]
  if(length(unique(popclust))==1){
    accuracy <- 10/30
  } else if (length(unique(popclust))==2){
    pop1 <- sum(ind_assignments[1:10]==popclust[1])
    pop2 <- sum(ind_assignments[11:20]==popclust[2])
    pop3 <- sum(ind_assignments[21:30]==missingclust)
    accuracy <-(pop1+pop2+pop3)/30
  } else if (length(unique(popclust))==3){
    pop1 <- sum(ind_assignments[1:10]==popclust[1])
    pop2 <- sum(ind_assignments[11:20]==popclust[2])
    pop3 <- sum(ind_assignments[21:30]==popclust[3])
    accuracy <-(pop1+pop2+pop3)/30
  }
  
  #swap column names to minimize label switching in structure plots
  clustnames <- c("q1","q2","q3")
  newclustnames <- c(rep(NA,nlevels(factor(q$pop))))
  for(i in 1:3){
    d <- subset(q,pop==i)
    e <- colMeans(d[3:5])
    f <- as.numeric(names(e[which(e==max(e))]))
    if(is.na(newclustnames[f])){
      newclustnames[f] <- clustnames[i]
    }
  }
  newclustnames[which(is.na(newclustnames))] <- clustnames[which(clustnames %in% newclustnames==F)]
  colnames(q) <- c("id","pop",newclustnames)
  
  q$mac <- mac
  q$sim <- sim
  q$run <- run
  q$lnL <- lnL
  q$lnL_mean <- lnL_mean
  
  q_matrices[[j]] <- q
  sum_stats[[j]] <- c(sim,mac,run,alpha,lnL,lnL_mean,accuracy)
  j <- j+1
}
sum_stats <- do.call(rbind.data.frame,sum_stats)
names(sum_stats) <- c("sim","mac","run","alpha","lnL","lnL_mean","accuracy")
sum_stats$log_alpha <- log(sum_stats$alpha)

mean_q_list <- lapply(q_matrices,function(e) {
  ddply(e,"pop",summarize,q1=mean(q1),q2=mean(q2),q3=mean(q3))[-1]
})
q_dist_list <- sapply(mean_q_list,function(e){
  mean(dist(e))
})
sum_stats$`Population Discrimination` <- q_dist_list/1.414214


sum_stats_wide <- sum_stats
sum_stats2 <- melt(sum_stats,id.vars = c("sim","mac","run"))

#ridge plot
ggplot(data=subset(sum_stats2,variable %in% c("log_alpha","accuracy")),
       aes(y=factor(mac),x=value))+
  theme_minimal()+theme(strip.background = element_blank())+
  facet_wrap(~variable,scales="free")+
  scale_y_discrete(expand=c(0.01,0))+
  ylab("Minimum Minor Allele Count")+xlab("")+
  geom_density_ridges(alpha=0.5,rel_min_height=0)

#structure plots (highest lnL per mac for each simulation)
strplotdata <- do.call(rbind,q_matrices)
strplotdata$q_dist <- unlist(lapply(q_dist_list,function(e) unlist(rep(e,30))))
best_runs <- ddply(strplotdata,.(mac,sim),summarize,max_lnL=max(lnL),max_q_dist=max(q_dist),mean_lnL=max(lnL_mean))
#strplotdata <- subset(strplotdata,q_dist %in% best_runs$max_q_dist)
strplotdata <- subset(strplotdata,lnL %in% best_runs$max_lnL)
meltq <- melt(strplotdata[-c(9,10,11)],id.vars=c("id","pop","sim","mac","run"))
ggplot(data=meltq,aes(x=id,y=value,fill=variable))+
  facet_grid(mac~sim)+
  theme_minimal()+theme(axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        strip.background = element_blank(),
                        axis.text.x=element_blank(),
                        rect = element_blank())+
  ylab("mac")+xlab("")+
  scale_fill_manual(values = grey.colors(3)[c(2,1,3)])+
  geom_bar(stat="identity",width=1)

#exploratory plots for likelihoods and pop dist
df <- sum_stats_wide 

df <- ddply(df,.(mac),function(e) subset(e,lnL>quantile(lnL,.2)))
ggplot(df,aes(x=lnL,y=`Population Discrimination`))+
  facet_wrap(~mac,scales="free")+
  theme_minimal()+theme(strip.background = element_blank())+
  geom_point()+
  geom_smooth(method="lm")
#so... the highest level of pop discrimination (and lowest alpha) is *not* from the highest or lowest lnL runs
