setwd("~/Dropbox/structure_simulations/")
library(plyr);library(data.table);library(foreach);library(doMC);library(ggplot2);library(magrittr);library(ggridges)
source("./R_scripts/maf/ggthemes.R")

mig_rates <- c(0,5e-10,5e-8,5e-6,5e-4,5e-2)
og_params <- readLines("fsc_params/sim_k3.par")

for(i in mig_rates){
  par <- og_params
  mig_matrix <- par[(grep("Migration matrix 0",par)+1):(grep("Migration matrix 0",par)+3)]
  mig_matrix <- lapply(mig_matrix,function(e) gsub("0.00005",i,e))
  par[(grep("Migration matrix 0",par)+1):(grep("Migration matrix 0",par)+3)] <- mig_matrix
  writeLines(unlist(par),paste0("./fsc_params/mig_",i,".par"))
  Sys.sleep(1) 
  system(paste0("cd ~/Dropbox/structure_simulations/;
                /Applications/fsc_mac64/fsc25221 -i ",paste0("./fsc_params/mig_",i,".par")," -n 10;
                mv ",paste0("mig_",i,"/*")," sample_migration/sim/ ;
                rm -r ",paste0("mig_",i))) #edit sim.par (or new par) as needed.
}

#convert to structure format
setwd("sample_migration/")

files <- list.files("sim",full.names = T)
files <- files[grep(".arp",files)]
j <- 1
for(i in files){
  arp <- readLines(i)
  tmp <- grep("SampleData",arp)
  arp <- arp[(tmp[1]+1):(tmp[3]+20)] #last value here should be n samples
  arp <- arp[-c(21:25,46:50)] #nsamples+1:nsamples+5
  arp <- gsub("*\t\t|*\t1\t","",arp)
  for(k in seq(2,60,2)){ #second value is 2*nsamples
    arp[k] <- paste0(unlist(strsplit(arp[(k-1)]," "))[1],arp[k])
  }
  names <- unlist(lapply(arp,function(e) unlist(strsplit(e," "))[1]))
  names <- data.frame(names,
                      pop=sapply(strsplit(names,"_"),function(e) as.numeric(e[1])))
  seq <- unlist(lapply(arp,function(e) unlist(strsplit(e," "))[2]))
  seq <- strsplit(seq,"")
  str <- do.call(rbind.data.frame, seq)
  
  #simple missing data simulation - drop 
  # str <- apply(str,2,function(e) {
  #   col <- as.numeric(as.character(e))
  #   drop <- sample(1:60,15)
  #   col[drop] <- -9
  #   col})
  # str <- data.frame(str)
  
  #join names and rename columns
  colnames(str) <- 1:ncol(str)
  str <- cbind(names,str)
  
  #write to file
  write.table(str,paste0("./str_in/",basename(i),".str"), 
              sep="\t",row.names = F,col.names = F,quote = F,append=F)
  j <- j+1
}

######################################
### generate maf-filtered datasets ###
######################################
library(magrittr);library(data.table);library(foreach);library(doMC)
setwd("~/Dropbox/structure_simulations/sample_migration/")
registerDoMC(cores=8)

#function to filter by minor allele count
filter_by_mac <- function(infile,mac,pop.info=T){
  for(maf in mac){
    tmp <- fread(infile,header=F) %>% data.frame()
    names <- tmp[1:(1+pop.info)]
    seq <- tmp[-(1:(1+pop.info))]
    
    # seq <- apply(seq,2,function(e){ #drop non-biallelic loci, make all sites 0/1 (skip if using simulation w/o missing data)
    #   e[e==-9] <- NA
    #   e <- factor(e)
    #   if(nlevels(e)==2){
    #     levels(e) <- c(0,1)
    #     e <- as.numeric(as.character(e))
    #     e[is.na(e)] <- -9
    #     e
    #   }
    # })
    # seq <- seq[!sapply(seq,is.null)] #drop null columns (non-biallelic sites)
    # seq <- do.call(cbind.data.frame,seq)
    
    mafs <- apply(seq,2,function(e) { #build index of minor allele frequencies/counts
      e[e==-9] <- NA
      allele_count <- sum(e,na.rm=T)
      if((allele_count/length(na.omit(e)))>0.5){ #if the derived allele is has frequency > 0.5, take 1-frequency
        allele_count <- length(na.omit(e))-allele_count
      }
      allele_count})
    
    #filter and write to file
    pi_seq <- seq[,mafs>=maf]
    pi_str <- cbind(names,pi_seq)
    write.table(pi_str,
                paste0(infile,
                       "_maf",maf,"_",
                       ".str"),
                sep="\t",row.names = F,col.names = F,quote = F,append = F)
  }
}

#get list of files
files <- list.files("./str_in",full.names = T)
files <- files[!grepl("maf",files)]
files <- grep(".str",files,value=T)

#run in parallel
foreach(i=files) %dopar% filter_by_mac(infile = i,mac=c(2),pop.info=T)

#####################################
######### run in structure ##########
#####################################
#multithreaded structure runs
setwd("/media/burke/bigMac/Dropbox/structure_simulations/")
#setwd("~/Dropbox/structure_simulations/sample_growth/")
library(foreach);library(doMC);library(data.table)
registerDoMC(cores=30)

#set options and file paths here (use full file paths. no ~ or .)
reps <- 1 #number of independent analyses per input file
structure_path <- "/media/burke/bigMac/Dropbox/structure_kernel_src/structure"   #"/Applications/structure/structure"
mainparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/sample_migration/str_params/mainparams.txt"
extraparams_path <- "/media/burke/bigMac/Dropbox/structure_simulations/sample_migration/str_params/extraparams.txt"
params_dir <- "/media/burke/bigMac/Dropbox/structure_simulations/sample_migration/str_params/"
str_in <- "/media/burke/bigMac/Dropbox/structure_simulations/sample_migration/str_in/"
str_out <- "/media/burke/bigMac/Dropbox/structure_simulations/sample_migration/str_out/"

#scan files in str_in for n loci
files <- list.files(str_in)
files <- grep(".str",files,value = T)
n_loci <- c()
for(i in files){
  tmp <- fread(paste0(str_in,i)) #can use fread() for significant speed increase if str files don't have extra whitespace columns
  n_loci <- append(n_loci,ncol(tmp)-2) #-1 if no pop info, -2 if pop info present.
}
names(n_loci) <- files
#tmp <- data.frame(n_loci,growth=unlist(lapply(names(n_loci),function(e) strsplit(e,"_") %>% unlist() %>% .[2]))) for debugging

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
    ex_params[grep("SEED",ex_params)] <- paste0("#define SEED     ",sample(1:1e6,1))
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

####################################
## post-hoc analysis and plotting ##
####################################
setwd("~/Dropbox/structure_simulations/sample_migration/")
files <- list.files("str_out",full.names = T)
q_matrices <- list()
sum_stats <- data.frame(matrix(ncol=3))
for(i in 1:length(files)){
  str <- readLines(files[i])
  mig <- strsplit(files[i],"_") %>% unlist() %>% .[3] %>% as.numeric()
  run <- strsplit(files[i],"_") %>% unlist() %>% .[5] %>% strsplit("\\.") %>% unlist() %>% .[1] %>% as.numeric()
  if(grepl("maf",files[i])){
    maf <- strsplit(files[i],"maf") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
  } else {
    maf <- 1
  }
  q <- str[(grep("Inferred ancestry",str)+2):(grep("Inferred ancestry",str)+31)]
  q <- lapply(q,function(e) strsplit(e," * ") %>% unlist() %>% .[7:9] %>% as.numeric())
  q <- do.call(rbind.data.frame,q)
  names(q) <- c("q1","q2","q3")
  q$pop <- c(rep(1,10),rep(2,10),rep(3,10))
  q$mig <- mig
  q$run <- run
  q$id <- 1:30
  
  #swap column order to minimize label switching
  tmp <- ddply(q,.(pop),function(e) colMeans(e[,1:3]))
  popclust <- c(which.max(tmp[1,2:4]),which.max(tmp[2,2:4]),which.max(tmp[3,2:4]))
  if(length(unique(popclust))!=3){
    popclust[duplicated(popclust)] <- c(1,2,3)[c(1,2,3) %in% c(popclust) == F]
  }
  q <- q[,c(popclust,4:ncol(q))]
  
  q_matrices[[i]] <- q
  sum_stats[i,] <- c(mig=mig,run=run,maf=maf)
}

mean_q_list <- lapply(q_matrices,function(e) {
  ddply(e,"pop",summarize,q1=mean(q1),q2=mean(q2),q3=mean(q3))[-1]
})
q_dist_list <- sapply(mean_q_list,function(e){
  mean(dist(e))
})
sum_stats$q_dist <- q_dist_list/max(q_dist_list)
names(sum_stats) <- c("mig","run","maf","q_dist")

strplotdata <- do.call(rbind,q_matrices)
meltq <- melt(strplotdata,id.vars=c("id","pop","mig","run"))
ggplot(data=meltq,aes(x=id,y=value,fill=variable))+facet_grid(run~mig)+
  theme_minimal()+theme(axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        strip.background = element_blank(),
                        strip.text.y = element_blank(),
                        axis.text.x=element_text(angle=45,hjust=1,size=7.5))+
  ylab("")+
  scale_fill_manual(values = grey.colors(3))+
  geom_bar(stat="identity")

ggplot(data=sum_stats,aes(y=factor(mig),x=q_dist,fill=factor(maf)))+
  theme_minimal()+theme(legend.position = "right")+
  xlim(-0.1,1.1)+
  scale_y_discrete(expand=c(0.025,0))+
  geom_density_ridges(alpha=0.5,stat="binline",aes(height=..density..),scale=.9)

