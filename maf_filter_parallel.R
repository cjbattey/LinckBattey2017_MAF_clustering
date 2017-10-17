#setwd("/media/burke/bigMac/Dropbox/structure_simulations/")
library(magrittr);library(data.table);library(foreach);library(doMC)
setwd("~/Dropbox/structure_simulations/")
registerDoMC(cores=8)

#function to filter by minor allele count
filter_by_mac <- function(infile,mac=c(1,2,3,4,5,8,20),pop.info=T){
  for(maf in mac){
  tmp <- fread(infile,header=F) %>% data.frame()
  names <- tmp[1:(1+pop.info)]
  seq <- tmp[-(1:(1+pop.info))]
  
  seq <- apply(seq,2,function(e){ #drop non-biallelic loci, make all sites 0/1 (skip if using simulation w/o missing data)
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
                     "_maf",maf,
                     ".str"),
              sep="\t",row.names = F,col.names = F,quote = F,append = F)
  }
}

#get list of files
files <- list.files("./str_in/k2",full.names = T)
files <- files[!grepl("maf",files)]
files <- grep(".str",files,value=T)

#run in parallel
foreach(i=files) %dopar% filter_by_mac(infile = i,mac=c(1,2,3,4,5,10,15),pop.info=T)
