nreps <- 100

truemodelbetter <- logical(nreps)

for (i in 1:nreps) {
  # initialize things for calculations
  p <- vector("list",2)
  individ <- vector("list",2)
  individ_nominors <- vector("list",2)
  ztrue <- vector("list",2)
  qtrue <- vector("list",2)
  pzgivenq <- vector("list",2)
  pxgivenall <- vector("list",2)
  pxmaf <- vector("list",2)
  
  pfalse <- vector("list",2)
  zfalse <- vector("list",2)
  qfalse <- vector("list",2)
  
  ## Let's just say there are 100 loci in total, and 2 populations
  ## Each population has a single site frequency spectrum that repeats at every locus, for convenience
  ## There are 3 alleles in each pop at freq >= 0.1, two "minor" alleles at < 0.1
  
  # Allele freqs in pop 1 at all loci
  # p[[1]] <- c(0.65,0.15,0.1,0.09,0.01)
  # Allele freqs in pop 1 at all loci
  # p[[2]] <- c(0.6,0.2,0.1,0.01,0.09)
  
  p1[[1]] <- c(0.6,0.4) #diallelic common allele freqs for p1
  p1[[2]] <- c(0.02,0.98) #diallelic rare allele freqs for p1
  
  p2[[1]] <- c(0.70,0.30)
  p2[[2]] <- c(0.08,0.92)
  
  # Simulate individuals with minor alleles
  # individ[[1]] <- sapply(1:100,function(i){sample(x=1:5,size=1,replace=TRUE,prob=p[[1]])})
  # individ[[2]] <- sapply(1:100,function(i){sample(x=1:5,size=1,replace=TRUE,prob=p[[2]])})
  
  individ_common <- sapply(1:70,function(i){sample(x=0:1,size=1,replace=TRUE,prob=p1[[1]])})
  individ_rare <- sapply(1:30,function(i){sample(x=0:1,size=1,replace=TRUE,prob=p1[[2]])})
  individ[[1]] <- c(individ_common, individ_rare)
  
  individ_common <- sapply(1:70,function(i){sample(x=0:1,size=1,replace=TRUE,prob=p2[[1]])})
  individ_rare <- sapply(1:30,function(i){sample(x=0:1,size=1,replace=TRUE,prob=p2[[2]])})
  individ[[2]] <- c(individ_common, individ_rare)
  
  ztrue[[1]] <- rep(1,100)
  ztrue[[2]] <- rep(2,100)
  
  qtrue[[1]] <- c(1,0)
  qtrue[[2]] <- c(0,1)
  
  # Remove the last two loci (minor alleles)
  # p[[1]] <- p[[1]][1:3]
  # p[[2]] <- p[[2]][1:3]
  
  #p <- lapply(p,function(vec){vec/sum(vec)})
  
  #nominors <- intersect(which(individ[[1]] < 4),which(individ[[2]] < 4))
  
  individ_nominors[[1]] <- individ[[1]][1:70] # simulate applying a
  individ_nominors[[2]] <- individ[[2]][1:70]
  
  # Components of likelihood under the true frequencies
  pzgivenq[[1]] <- prod(qtrue[[1]][ztrue[[1]]])
  pzgivenq[[2]] <- prod(qtrue[[2]][ztrue[[2]]])
  
  #pxgivenall[[1]] <- prod(p[[1]][individ[[1]]])
  #pxgivenall[[2]] <- prod(p[[2]][individ[[2]]])
  
  # with all alleles
  pxgivenall[[1]] <- prod(p1[[1]][individ[[1]][1:70]])*prod(p1[[2]][individ[[1]][71:100]])
  pxgivenall[[2]] <- prod(p2[[1]][individ[[2]][1:70]])*prod(p2[[2]][individ[[2]][71:100]])
  
  # filtered
  pxmaf[[1]] <- prod(p1[[1]][individ_nominors[[1]]])
  pxmaf[[2]] <- prod(p2[[1]][individ_nominors[[2]]])
  
  # the log-likelihood of the simulating model
  lnLtruenominors <- sum(log(unlist(pxgivenall))) + sum(log(unlist(pzgivenq)))
  
  # for the filtered model 
  lnLtruenominors2 <- sum(log(unlist(pxmaf))) + sum(log(unlist(pzgivenq)))
  
  # build the incorrect, smear model
  # zfalse[[1]] <- ifelse(individ[[1]] <= 2,1,2)
  # zfalse[[2]] <- ifelse(individ[[2]] <= 2,1,2)
  
  zfalse[[1]] <- ifelse(individ[[1]] <= 0,1,0)
  zfalse[[2]] <- ifelse(individ[[2]] <= 0,1,0)
  
  # qfalse[[1]] <- c(length(which(zfalse[[1]] == 1))/length(nominors),1-length(which(zfalse[[1]] == 1))/length(nominors))
  # qfalse[[2]] <- c(length(which(zfalse[[2]] == 1))/length(nominors),1-length(which(zfalse[[2]] == 1))/length(nominors))
  
  qfalse[[1]] <- c(length(which(zfalse[[1]] == 1))/70,1-length(which(zfalse[[1]] == 1))/70)
  qfalse[[2]] <- c(length(which(zfalse[[2]] == 1))/70,1-length(which(zfalse[[2]] == 1))/70)
  
  
  lnLtruenominors <- sum(log(unlist(pxgivenall))) + sum(log(unlist(pzgivenq)))
  
  #avgoccurences <- c((p[[1]][1]+p[[2]][1])/2,(p[[1]][2]+p[[2]][2])/2)
  #avgoccurences <- avgoccurences/sum(avgoccurences)
  
  avgoccurences <- c((p1[[1]][1]+p2[[1]][1])/2,(p1[[2]][1]+p2[[2]][1])/2)
  avgoccurences <- avgoccurences/sum(avgoccurences)
  
  pfalse[[1]] <- sapply(individ[[1]],function(x){
    if ( x == 1 ) {
      return(avgoccurences[1])
    } else if (x == 3) {
      return(avgoccurences[2])
    } else {
      return(1)
    }
  })
  
  pfalse[[2]] <- sapply(individ[[2]],function(x){
    if ( x == 1 ) {
      return(avgoccurences[1])
    } else if (x == 3) {
      return(avgoccurences[2])
    } else {
      return(1)
    }
  })
  
  # likelihood components under the smear model
  pzgivenq[[1]] <- prod(qfalse[[1]][zfalse[[1]]])
  pzgivenq[[2]] <- prod(qfalse[[2]][zfalse[[2]]])
  
  pxgivenall[[1]] <- prod(pfalse[[1]][zfalse[[1]]][individ[[1]]])
  pxgivenall[[2]] <- prod(pfalse[[2]][zfalse[[2]]][individ[[2]]])
  
  # log-likelihood of smear model
  lnLfalsenominors <- sum(log(unlist(pxgivenall))) + sum(log(unlist(pzgivenq)))
  
  # compare log likelihoods
  truemodelbetter[i] <- lnLtruenominors > lnLfalsenominors
}

length(which(truemodelbetter))
