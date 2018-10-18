#3-population simulation and structure inference pipeline
library(magrittr);library(data.table);library(foreach);library(doMC);library(cowplot)
library(ggplot2);library(plyr);library(ggridges);library(reshape);library(gridExtra)
setwd("~/Dropbox/structure_simulations/")
registerDoMC(cores=10)
source("./R_scripts/maf/maf_functions.R")

####################################
######## simulate data #############

#run fsc simulation using parameters in ./fsc_params/sim.par
system("cd ~/Dropbox/structure_simulations/;
       /Applications/fsc_mac64/fsc25221 -i ./fsc_params/sim.par -g -n 10") #edit sim.par (or new par) as needed.

#convert fastsimcoal (arlequin) to structure format
files <- list.files("sim",full.names = T)
files <- files[grep(".arp",files)]
foreach(i=files) %dopar% arp2structure(infile=i,md=.25,npops=3,samples_per_pop=10,out_directory="str_in/5kloci/")

#Filter sites by minimum minor allele count
files <- list.files("./str_in/5kloci/",full.names = T)
foreach(i=files) %dopar% filter_by_mac(infile = i,mac=c(2,3,4,5,8,11,15),pop.info=T,max_md=.25)

# #randomly sample sites to 1000bp for each alignment
# files <- list.files("./str_in",full.names = T)
# foreach(i=files) %dopar% sample_sites(i,1000)

#########################################
### run structure in parallel on wopr ###
setwd("/media/burke/bigMac/Dropbox/structure_simulations/")
library(foreach);library(doMC);library(data.table)
registerDoMC(cores=20)

#set options and file paths here (use full paths)
reps <- 10 #number of independent analyses per input file
structure_path <- "/media/burke/bigMac/Dropbox/tools/structure_kernel_src/structure"   #"/Applications/structure/structure"
mainparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/str_params/mainparams.txt"
extraparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/str_params/extraparams.txt"
params_dir <- "/media/burke/bigMac/Dropbox/structure_simulations/str_params/"
str_in <- "/media/burke/bigMac/Dropbox/structure_simulations/str_in/5kloci/"
str_out <- "/media/burke/bigMac/Dropbox/structure_simulations/str_out/5kloci/"

#scan files in str_in for n loci
files <- list.files("str_in/5kloci/",full.names = T)
files <- files[grepl("subsample",files)]
n_loci <- c()
for(i in files){
  tmp <- fread(i) #use fread() for significant speed increase if str files don't have extra whitespace columns (otherwise read.table())
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

# ###########################################
# ######## faststructure ####################
# #strip taxa names and population column, add six empty columns to structure input to match faststructure input reqs (srsly...)
# files <- list.files("str_in/5kloci",full.names = T)
# str2faststr <- function(file){
#   str <- read.table(file)
#   str <- str[,-c(1:2)] #use this row if there's a population column
#   #str <- str[,-1] #use this if no population column
#   blank <- data.frame(matrix(nrow=nrow(str),ncol=6,data="faststructuremademeputthishere"))
#   str <- cbind(blank,str)
#   outname <- basename(file) %>% tools::file_path_sans_ext()
#   write.table(str,paste0("./fstr_in/",outname,".str"),row.names = F,col.names = F)
# }
# 
# foreach(i=files) %dopar% str2faststr(i)
# 
# #build list of commands to run faststructure in parallel
# files <- list.files("fstr_in",full.names = T)
# commands <- c()
# nreps <- 10
# for(i in files){
#   for(j in 1:nreps){
#     outname <- basename(i) %>% tools::file_path_sans_ext()
#     outname <- paste0(outname,"_",j)
#     command <- paste0("python /Applications/fastStructure/structure.py -K 3 --format=str --input=/users/cj/Dropbox/structure_simulations/",
#                       tools::file_path_sans_ext(i),
#                       " --output=/users/cj/Dropbox/structure_simulations/fstr_out/",
#                       outname,
#                       " --seed=",sample(1:1e6,1))
#     commands <- append(commands,command)
#   }
# }
# 
# #run in parallel
# foreach(i=commands) %dopar% system(i)

###########################################
######## multivariate clustering ##########

files <- list.files("str_in/40kloci",full.names = T)
files <- files[grepl("subsample",files)]
clust <- foreach(i=files,.combine = rbind) %dopar% cluster_multivar(i,pop.info=T,nreps=10)
write.table(clust,"sim_mv_clust_constlength.csv",sep=",",row.names = F,col.names = T,quote = F)

#ridge plot
clust <- clust[,-2]
melt_clust <- melt(clust,id.vars = c("sim","mac","run"))
ridgeplot <- ggplot(data=melt_clust,aes(x=value,y=factor(mac)))+
  theme_classic()+
  theme(strip.background = element_blank())+
  facet_wrap(~variable,scales="free_x")+
  ylab("Minimum minor allele count")+xlab("")+
  scale_y_discrete(expand = c(0.01,0))+
  geom_density_ridges()

#PC plots
tmp <- files[grepl("sim4",files)]
pcs <- get_pc_from_structure(tmp,pop.info=T)
pcs <- do.call(rbind.data.frame,pcs)
pcs$mac[is.na(pcs$mac)] <- 1
pcs$mac <- factor(pcs$mac,levels=c(15,11,8,5,4,3,2,1))
pcplots <- ggplot(data=pcs,aes(x=PC1,y=PC2,col=clust))+
  theme_classic()+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        axis.title.y.right=element_text())+
  scale_x_continuous(breaks=c(-5,0,5))+scale_y_continuous(breaks=c(-5,5))+
  scale_color_grey()+
  facet_grid(mac~.)+
  geom_point(size=0.5)+stat_ellipse()

pdf("fig/sim_pc_constlen.pdf",width=6.5,height=4)
ggdraw()+
  draw_plot(ridgeplot,0,0,.75,1)+
  draw_plot(pcplots,.75,0,.25,.875)
dev.off()                               
                                                  

#########################################
###### summarize structure output #######
files <- list.files("./str_out/5kloci/",full.names=T)# %>% grep("subsample",.,value=T)
q_matrices <- list()
sum_stats <- list()
j <- 1
for(i in files){
  tmp <- readLines(i,warn=F)
  if(grepl("subsample",i)){
    #version for subsampled data
    if(!grepl("mac",i)){
      mac <- 1
      run <- i %>% strsplit("_") %>% unlist() %>% .[4] %>% as.numeric()
      sim <- i %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    } else{
      mac <- i %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit("\\_") %>% unlist() %>% .[1] %>% as.numeric()
      run <- i %>% strsplit("_") %>% unlist() %>% .[5] %>% as.numeric()
      sim <- i %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    }
  } else {
    #version for nonsubsampled data (slightly changed file naming bc I'm dumb)
    if(!grepl("mac",i)){
      mac <- 1
      run <- i %>% strsplit("_") %>% unlist() %>% .[3] %>% as.numeric()
      sim <- i %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit(".str") %>% unlist() %>% .[1] %>% as.numeric()
    } else{
      mac <- i %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit(".str") %>% unlist() %>% .[1] %>% as.numeric()
      run <- i %>% strsplit("_") %>% unlist() %>% .[4] %>% as.numeric()
      sim <- i %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    }
  }
  alpha <- tmp[grep("Mean value of alpha",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
  lnL <- tmp[grep("Estimated Ln Prob of Data",tmp)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
  
  q <- tmp[(grep("Inferred ancestry of individuals:",tmp)+2):(grep("Inferred ancestry of individuals:",tmp)+31)] %>% 
        lapply(function(e){ strsplit(e," * ") %>% unlist() %>% .[-c(1,2,4,6,10,11,12)]}) %>% 
          do.call(rbind.data.frame,.)
  colnames(q) <- c("id","pop","1","2","3")
  q[,3:5] <- q[,3:5] %>% apply(2,function(e) as.numeric(as.character(e)))
  
  #swap column names to minimize label switching in structure plots
  clustnames <- c("q1","q2","q3")
  newclustnames <- c(rep(NA,nlevels(factor(q$pop))))
  for(k in 1:3){
    d <- subset(q,pop==k)
    e <- colMeans(d[3:5])
    f <- as.numeric(names(e[which(e==max(e))])[1])
    if(is.na(newclustnames[f])){
      newclustnames[f] <- clustnames[k]
    }
  }
  newclustnames[which(is.na(newclustnames))] <- clustnames[which(clustnames %in% newclustnames==F)]
  colnames(q) <- c("id","pop",newclustnames)
  
  #get pc_st (within-population distance/total distance)
  pcst <- pc_st(q[,c("q1","q2","q3")],q$pop)
  
  #get assignment accuracy (% individuals assigned to "correct" population)
  assignments <- apply(q[,c("q1","q2","q3")],1,function(e){
    match(max(e),e)
  })
  accuracy <- sum(assignments==c(rep(1,10),rep(2,10),rep(3,10)))/30
  
  #extra columns for convenient q matrix filtering (actually less efficient? tbd)
  q$pcst <- pcst
  q$mac <- mac
  q$sim <- sim
  q$run <- run
  q$lnL <- lnL
  #q$lnL_mean <- lnL_mean
  
  q_matrices[[j]] <- q
  sum_stats[[j]] <- c(sim,mac,run,alpha,lnL,pcst,accuracy)
  j <- j+1
}
sum_stats <- do.call(rbind.data.frame,sum_stats)
names(sum_stats) <- c("sim","mac","run","alpha","lnL","pcst","accuracy")
sum_stats$log_alpha <- log(sum_stats$alpha)
sum_stats_wide <- sum_stats
sum_stats2 <- melt(sum_stats,id.vars = c("sim","mac","run"))

#ridge plot
ridgeplot <- ggplot(data=subset(sum_stats2,variable %in% c("accuracy","pcst")),
       aes(y=factor(mac),x=value))+
  theme_classic()+theme(strip.background = element_blank())+
  facet_wrap(~variable,scales="free_x")+
  scale_y_discrete(expand=c(0.07,0))+
  ylab("Minimum Minor Allele Count")+xlab("")+
  geom_density_ridges(alpha=0.5,rel_min_height=0)

#structure plots (highest lnL per mac for each simulation)
strplotdata <- do.call(rbind,q_matrices)
best_runs <- ddply(subset(strplotdata,sim==5),.(mac),function(e){
  #z <- e$run[which(max(e$lnL)==e$lnL)][1]
  #x <- e$sim[which(max(e$lnL)==e$lnL)][1]
  z <- e$run[which(max(e$pcst)==e$pcst)][1]
  x <- e$sim[which(max(e$pcst)==e$pcst)][1]
  subset(e,e$run==z & sim==x)
})
meltq <- melt(best_runs[-c(6,10)],id.vars=c("id","pop","sim","mac","run"))
#meltq <- melt(strplotdata[-c(6,10)],id.vars=c("id","pop","sim","mac","run"))
meltq$mac <- factor(meltq$mac, levels=rev(levels(factor(meltq$mac))))
strplot <- ggplot(data=meltq,aes(x=id,y=value,fill=variable))+
              facet_grid(mac~.)+
              ggtitle("  ")+
              theme_classic()+theme(axis.text.y=element_blank(),
                                    axis.ticks=element_blank(),
                                    strip.background = element_blank(),
                                    axis.text.x=element_blank(),
                                    strip.text.y=element_blank(),
                                    rect = element_blank(),
                                    legend.position = "none")+
              ylab("")+xlab("")+
              scale_fill_manual(values = grey.colors(3)[c(2,1,3)])+
              geom_bar(stat="identity",width=.9,col="black",lwd=0.25)

pdf("~/Dropbox/structure_simulations/fig/sim_struct_varlen.pdf",width=6,height=3.5)
ggdraw()+
  draw_plot(ridgeplot,0,0,.7,1)+
  draw_plot(strplot,0.7,.075,.3,.925)
dev.off()


#full structure plot per sim
meltq <- melt(strplotdata[-c(6,10)],id.vars=c("id","pop","sim","mac","run"))
meltq$mac <- factor(meltq$mac, levels=rev(levels(factor(meltq$mac))))
ggplot(data=subset(meltq,sim==5),aes(x=id,y=value,fill=variable))+
  facet_grid(mac~run)+
  ggtitle("  ")+
  theme_classic()+theme(axis.text.y=element_blank(),
                        axis.ticks=element_blank(),
                        strip.background = element_blank(),
                        axis.text.x=element_blank(),
                        strip.text.y=element_blank(),
                        rect = element_blank(),
                        legend.position = "none")+
  ylab("")+xlab("")+
  scale_fill_manual(values = grey.colors(3)[c(2,1,3)])+
  geom_bar(stat="identity",width=.9,col="black",lwd=0.25)




