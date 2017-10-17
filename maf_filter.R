### DEPRECATED. See maf_filter_parallel.R

#fastsimcoal -> structure analysis pipeline
#setwd("/media/burke/bigMac/Dropbox/structure_simulations/")
library(magrittr);library(data.table)
setwd("~/Dropbox/structure_simulations/")

pop.info=T #T/F for presence of pop info column

#filter for minor allele frequency and write to file 
files <- list.files("./str_in/")
files <- files[!grepl("maf",files)]
files <- grep(".str",files,value=T)
#for(maf in c(1/80,2/80,3/80,4/80,.1,.25)){
for(maf in c(1,2,3,4,5,8,20)){
  j <- 1
  for(i in files){
    tmp <- fread(paste0("./str_in/",i),header=F) %>% data.frame()
    names <- tmp[1:(1+pop.info)]
    seq <- tmp[-(1:(1+pop.info))]

    seq <- apply(seq,2,function(e){ #drop non-biallelic loci, make all sites 0/1
      e[e==-9] <- NA
      e <- factor(e)
      if(nlevels(e)==2){
        levels(e) <- c(0,1)
        e <- as.numeric(as.character(e))
        e[is.na(e)] <- -9
        e
      }
    })
    seq <- seq[!sapply(seq,is.null)] 
    seq <- do.call(cbind.data.frame,seq)
    
    mafs <- apply(seq,2,function(e) { #build index of minor allele frequencies/counts
      e[e==-9] <- NA
      af <- sum(e,na.rm=T)
      if((af/length(!is.na(e)))>0.5){
        af <- length(!is.na(e))-af
      }
      af})
    
    #filter and write to file
    pi_seq <- seq[,mafs>=maf]
    pi_str <- cbind(names,pi_seq)
    write.table(pi_str,
                paste0("./str_in/",i,
                       "_maf",maf,
                       ".str"),
                sep="\t",row.names = F,col.names = F,quote = F,append = F)
    j <- j+1
  }
}
