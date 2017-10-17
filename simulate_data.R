#fastsimcoal -> structure analysis pipeline
#setwd("/media/burke/bigMac/Dropbox/structure_simulations/")
setwd("~/Dropbox/structure_simulations/")

#run fsc simulation using parameters in ./fsc_params/sim.par
system("cd ~/Dropbox/structure_simulations/;
        /Applications/fsc_mac64/fsc25221 -i ./fsc_params/sim_k2.par -g -n 100") #edit sim.par (or new par) as needed.

#convert fastsimcoal (arlequin) to structure format, write to sim0X.str
#if 3 populations of 10 samples each
files <- list.files("sim/")
files <- files[grep(".arp",files)]
j <- 1
for(i in files){
  arp <- readLines(paste0("sim/",i))
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
  str <- apply(str,2,function(e) {
    col <- as.numeric(as.character(e))
    drop <- sample(1:60,15)
    col[drop] <- -9
    col})
  str <- data.frame(str)
  
  #join names and rename columns
  colnames(str) <- 1:ncol(str)
  str <- cbind(names,str)
  
  #write to file
  write.table(str,paste0("./str_in/sim",formatC(j,digits=1,flag="0",format = "d"),".str"), #formatC() adds leading 0's 
                  sep="\t",row.names = F,col.names = F,quote = F,append=F)
  j <- j+1
}


#if 2 populations of 15 each
files <- list.files("sim_k2/")
files <- files[grep(".arp",files)]
j <- 1
for(i in files){
  arp <- readLines(paste0("sim_k2/",i))
  tmp <- grep("SampleData",arp)
  arp <- arp[(tmp[1]+1):(tmp[2]+30)] #last value here should be nsamples*2
  arp <- arp[-c(31:35)] #nsamples+1:nsamples+5
  arp <- gsub("*\t\t|*\t1\t","",arp)
  for(k in seq(2,60,2)){ #second value is 2*nsamples
    arp[k] <- paste0(unlist(strsplit(arp[(k-1)]," "))[1],arp[k])
  }
  names <- unlist(lapply(arp,function(e) unlist(strsplit(e," "))[1]))
  names <- data.frame(names,
                      pop=sapply(strsplit(names,"_"),function(e) as.numeric(e[1])))
  seq <- unlist(lapply(arp,function(e) unlist(strsplit(e," "))[2]))
  seq <- strsplit(seq,"")
  #str <- do.call(rbind.data.frame, seq) #somehow the two lines below are faster
  str <- data.frame(matrix(nrow=60,ncol=length(seq[[1]])))
  for(i in 1:60){str[i,]<-seq[[i]]}
  
  # #missing data simulation
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
  write.table(str,paste0("./str_in/k2/sim_",formatC(j,digits=1,flag="0",format = "d"),".str"), #formatC() adds leading 0's 
              sep="\t",row.names = F,col.names = F,quote = F,append=F)
  j <- j+1
}
