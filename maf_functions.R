#functions for filtering, processing, and clustering structure-formatted genotype alignments files



#convert arp format from fastsimcoal2 output to structure format. requires same n inds per pop.
#md is the average fraction of individuals missing data per site.
arp2structure <- function(infile,out_directory="str_in/",md=NULL,npops=3,samples_per_pop=10){
  require(truncnorm);require(tools)
  arp <- readLines(infile)
  tmp <- list()
  for(i in grep("SampleData",arp)){
    tmp <- append(tmp,arp[(i+1):(i+samples_per_pop*2)])
  }
  arp <- tmp
  arp <- gsub("*\t\t|*\t1\t","",arp)
  for(k in seq(2,2*npops*samples_per_pop,2)){ 
    arp[k] <- paste0(unlist(strsplit(arp[(k-1)]," "))[1],arp[k])
  }
  names <- unlist(lapply(arp,function(e) unlist(strsplit(e," "))[1]))
  names <- data.frame(names,
                      pop=sapply(strsplit(names,"_"),function(e) as.numeric(e[1])))
  seq <- unlist(lapply(arp,function(e) unlist(strsplit(e," "))[2]))
  seq <- strsplit(seq,"")
  str <- do.call(rbind.data.frame, seq)
  
  #missing data simulation (assumes the proportion of missing individuals is normally distributed with mean=md & sd=md/2)
  if(!is.null(md)){
    str <- apply(str,2,function(e) {
      col <- as.numeric(as.character(e))
      n_inds_dropped <- round(rtruncnorm(1,mean=round(npops*samples_per_pop*md),sd=(npops*samples_per_pop*md)/2,a=0,b=npops*samples_per_pop))
      drop <- sample(1:(npops*samples_per_pop),n_inds_dropped)
      drop <- seq(1,60,2)[drop]
      drop <- lapply(drop,function(e) c(e,e+1)) %>% unlist()
      col[drop] <- -9
      col})
    str <- data.frame(str)
  }
  
  #join names and rename columns
  colnames(str) <- 1:ncol(str)
  str <- cbind(names,str)
  
  #write to file
  write.table(str,paste0(out_directory,gsub("_1_","",file_path_sans_ext(basename(infile))),".str"),
              sep="\t",row.names = F,col.names = F,quote = F,append=F)
}



#filter structure file by minor allele count 
#accepts structure files with -9 as the NA character
#two rows per individual, and columns: ID, population, genotypes...
#writes new mac-filtered files in the location of the infile appended with mac1/mac2/etc
filter_by_mac <- function(infile,mac=c(2,3,5,8,11,15),pop.info=T,na.char=-9,max_md){
  library(data.table);library(magrittr);require(pbapply)
  for(mac_selected in mac){
    tmp <- fread(infile,header=F) %>% data.frame()
    names <- tmp[1:(1+pop.info)]
    seq <- tmp[-(1:(1+pop.info))]
    
    seq <- pbapply(seq,2,function(e){ #drop non-biallelic loci, make all sites 0/1 (can skip if using simulation w/o missing data)
      e[e==na.char] <- NA
      md <- 1-(length(na.omit(e))/length(e))
      e <- factor(e)
      if(nlevels(e)==2 & md<=max_md){
        levels(e) <- c(0,1)
        e <- as.numeric(as.character(e))
        e[is.na(e)] <- -9
        e
      }
    })
    seq <- seq[!sapply(seq,is.null)]
    seq <- data.frame(seq)

    mac_index <- pbapply(seq,2,function(e) { #build index of minor allele counts
      e[e==-9] <- NA
      allele_count <- sum(e,na.rm=T)
      if((allele_count/length(na.omit(e)))>0.5){ #if the derived allele is has frequency > 0.5, take 1-frequency
        allele_count <- length(na.omit(e))-allele_count
      }
      allele_count})
    
    #filter and write to file
    pi_seq <- seq[,mac_index>=mac_selected]
    pi_str <- cbind(names,pi_seq)
    write.table(pi_str,
                paste0(tools::file_path_sans_ext(infile),
                       "_mac",mac_selected,
                       ".str"),
                sep="\t",row.names = F,col.names = F,quote = F,append = F)
  }
}

#randomly sample n sites out of the alignment and write to file
sample_sites <- function(infile,n,pop.info=T){
  tmp <- fread(infile,header=F,data.table=F)
  names <- tmp[1:(1+pop.info)]
  seq <- tmp[-(1:(1+pop.info))]
  seq <- seq[,sample(1:ncol(seq),n)]
  str <- cbind(names,seq)
  write.table(str,
              paste0(tools::file_path_sans_ext(infile),
                     "_subsample.str"),
              sep="\t",row.names=F,col.names=F,quote=F,append=F)
}

#mode function for checking k-means cluster assignments (via https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ratio of mean within-population distance to mean distance across all individuals
pc_st <- function(pc_matrix,populations){
  out <- c()
  for(k in unique(populations)){
    within_dist <- mean(dist(pc_matrix[populations==k,]))
    total_dist <- mean(dist(pc_matrix))
    ratio <- within_dist/total_dist
    out <- append(out,ratio)
  }
  return(1-mean(out))
}

# ratio if inter-individual to inter-population distance of population centers
# pc_st <- function(pc_matrix,populations){
#   pop_centers <- sapply(unique(populations),function(e){
#       sapply(1:3,function(f) mean(pc_matrix[populations==e,f]))
#   })
#   popdist <- mean(dist(pop_centers))
#   inddist <- mean(dist(pc_matrix))
#   return(1-inddist/popdist)
# }

#distance across populations with each axis scaled to 0-1
# pc_st <- function(pc_matrix,populations){
#   pop_centers <- sapply(unique(populations),function(e){
#       sapply(1:3,function(f) mean(pc_matrix[populations==e,f]))
#   })
#   popdist <- mean(dist(pop_centers))
#   inddist <- mean(dist(pc_matrix))
#   return(popdist/1.414214)
# }

#run k-means assignments and dapc x-validations on a designated infile (single-threaded)
cluster_multivar <- function(infile,pop.info=T,pop=NULL,nreps=1,satrapa=F){
  require(adegenet)
  tmp <- read.table(infile)
  nloci <- ncol(tmp)-(1+pop.info)
  nind <- nrow(tmp)/2
  if(pop.info==T){
    str <- read.structure(infile,n.ind=nind,n.loc=nloci,col.lab=1,col.pop=2,row.marknames=0,onerowperind=F,ask=F)
  } else {
    str <- read.structure(infile,n.ind=nind,n.loc=nloci,col.lab=1,col.pop=0,row.marknames=0,onerowperind=F,ask=F)
    str@pop <- factor(pop)
  }
  #run k-means and report the proportion of individuals that would be correctly grouped with their predefined populations
  kmeans_accuracy <- c();pcst <- c()
  for(i in 1:nreps){
    clust <- find.clusters(str,n.pca=length(str@pop),n.clust=3)$grp
    str_scaled <- scaleGen(str,NA.method="mean",center=T,scale=F)
    pc <- prcomp(str_scaled,center=F,scale=F)$x[,1:3] %>% data.frame()
    pc$clust <- unlist(clust)
    pcst[i] <- pc_st(pc[,1:3],str@pop)
    #ggplot(data=pc,aes(x=PC1,y=PC2,col=clust))+geom_point()+geom_point(data=pc)+stat_ellipse()
    
    popclust <- c(Mode(clust[str@pop==1]),Mode(clust[str@pop==2]),Mode(clust[str@pop==3]))
    if(length(unique(popclust)) < 3){ #if any pop lacks a unique modal cluster assignment, assign it one of the missing clusters
      missingclust <- c(1:3)[!c(1:3) %in% popclust]
      popclustct <- c(sum(as.integer(clust[str@pop==1]==popclust[1])),
                      sum(as.integer(clust[str@pop==2]==popclust[2])),
                      sum(as.integer(clust[str@pop==3]==popclust[3])))
      popclust[which.min(popclustct)] <- missingclust
    }
    kmeans_accuracy[i] <- sum(sum(as.integer(clust[str@pop==1]==popclust[1])),
                              sum(as.integer(clust[str@pop==2]==popclust[2])),
                              sum(as.integer(clust[str@pop==3]==popclust[3])))/length(str@pop)
  }

  #dapc using original (pre-defined) clusters, using n PC's w/highest success in first 10 reps. 
  str_noMD <- scaleGen(str,NA.method="mean",center=T,scale=F)
  xval <- xvalDapc(str_noMD,grp=str@pop,n.da=2,n.pca=3,center=F,scale=F,xval.plot=F,n.rep=nreps,training.set=0.5)
  xval <- xval$`Cross-Validation Results`[,2] %>% mean()

  if(!satrapa){
    if(grepl("subsample",infile)){
      if(!grepl("mac",infile)){
        mac <- 1
        sim <- infile %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
        out <- data.frame(kmeans=kmeans_accuracy,xval=xval,pcst=pcst,sim=sim,mac=mac,run=1:nreps)
        out
      } else{
        mac <- infile %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit("\\_") %>% unlist() %>% .[1] %>% as.numeric()
        sim <- infile %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
        out <- data.frame(kmeans=kmeans_accuracy,xval=xval,pcst=pcst,sim=sim,mac=mac,run=1:nreps)
        out
      }
    } else if(!grepl("subsample",infile)){
      if(!grepl("mac",infile)){
        mac <- 1
        sim <- infile %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
        out <- data.frame(kmeans=kmeans_accuracy,xval=xval,pcst=pcst,sim=sim,mac=mac,run=1:nreps)
        out
      } else{
        mac <- infile %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit("\\.") %>% unlist() %>% .[1] %>% as.numeric()
        sim <- infile %>% strsplit("sim") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
        out <- data.frame(kmeans=kmeans_accuracy,xval=xval,pcst=pcst,sim=sim,mac=mac,run=1:nreps)
        out
      }
    }
   } else if(satrapa){
    mac <- infile %>% strsplit("mac") %>% unlist() %>% .[2] %>% strsplit("\\.") %>% unlist() %>% .[1] %>% as.numeric()
    out <- data.frame(kmeans=kmeans_accuracy,xval=xval,pcst=pcst,mac=mac)
    out
  }
}

#get pc's for a set of structure infiles
get_pc_from_structure <- function(infiles,pop.info=T,pop=NULL,satrapa=F){
  require(adegenet)
  out <- list()
  for(i in 1:length(infiles)){
    tmp <- read.table(infiles[i])
    nloci <- ncol(tmp)-(1+pop.info)
    nind <- nrow(tmp)/2
    if(pop.info==T){
      str <- read.structure(infiles[i],n.ind=nind,n.loc=nloci,col.lab=1,col.pop=2,row.marknames=0,onerowperind=F,ask=F)
    } else {
      str <- read.structure(infiles[i],n.ind=nind,n.loc=nloci,col.lab=1,col.pop=0,row.marknames=0,onerowperind=F,ask=F)
      str@pop <- factor(pop)
    }
    str_scaled <- scaleGen(str,NA.method="mean",center=T,scale=F)
    pc <- prcomp(str_scaled,center=F,scale=F)$x[,1:3] %>% data.frame()
    pc$sim <- strsplit(infiles[i],"sim") %>% unlist() %>% strsplit("_|\\.") %>% unlist() %>% .[3] %>% as.numeric()
    if(satrapa){
      pc$mac <- strsplit(infiles[i],"mac") %>% unlist() %>% strsplit("_|\\.") %>% unlist() %>% .[4] %>% as.numeric()
    } else {
      pc$mac <- strsplit(infiles[i],"mac") %>% unlist() %>% strsplit("_|\\.") %>% unlist() %>% .[3] %>% as.numeric()
    }
    pc$clust <- find.clusters(str,n.pca=length(str@pop),n.clust=3)$grp
    out[[i]] <- pc
  }
  return(out)
}







