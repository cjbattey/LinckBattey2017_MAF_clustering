#multithreaded structure runs
setwd("/Users/cj/Dropbox/structure_simulations/")
library(foreach);library(doMC);library(data.table)
registerDoMC(cores=8)

#set options and file paths here (use full file paths. no ~ or .)
reps <- 1 #number of independent analyses per input file
structure_path <- "/Applications/structure/structure"
mainparams_path <- "/Users/cj/Dropbox/structure_simulations/str_params/k2/mainparams.txt"
extraparams_path <- "/Users/cj/Dropbox/structure_simulations/str_params/k2/extraparams.txt"
params_dir <- "/Users/cj/Dropbox/structure_simulations/str_params/"
str_in <- "/Users/cj/Dropbox/structure_simulations/str_in/k2/"
str_out <- "/Users/cj/Dropbox/structure_simulations/str_out/k2/"

#scan files in str_in for n loci
files <- list.files(str_in)
files <- grep(".str",files,value = T)
files <- grep("maf",files,value = T)
n_loci <- c()
for(i in files){
  tmp <- fread(paste0(str_in,i)) #can use fread() for significant speed increase if str files don't have extra whitespace columns
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


