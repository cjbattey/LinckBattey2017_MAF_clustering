#faststructure parallelization script
#setwd("/Users/cj/Dropbox/structure_simulations/satrapa/")
setwd("/media/burke/bigMac/Dropbox/structure_simulations/satrapa")
library(foreach);library(doMC);library(data.table);library(tidyr)
registerDoMC(cores=8)

#strip taxa names and population column, add six empty columns to structure input to match faststructure input reqs (srsly...)
files <- list.files("str_in",full.names = T)
str2faststr <- function(file){
  str <- read.table(file)
  #str <- str[,-c(1:2)] #use this row if there's a population column
  str <- str[,-1] #use this if no population column
  blank <- data.frame(matrix(nrow=nrow(str),ncol=6,data="faststructuremademeputthishere"))
  str <- cbind(blank,str)
  outname <- basename(file) %>% tools::file_path_sans_ext()
  write.table(str,paste0("./fstr_in/",outname,".str"),row.names = F,col.names = F)
}

foreach(i=files) %dopar% str2faststr(i)

#build list of commands to run faststructure in parallel
files <- list.files("fstr_in",full.names = T)
commands <- c()
nreps <- 10
for(i in files){
  for(j in 1:nreps){
  outname <- basename(i) %>% tools::file_path_sans_ext()
  outname <- paste0(outname,"_",j)
  command <- paste0("python ~/fastStructure/structure.py -K 3 --format=str --input=/media/burke/bigMac/Dropbox/structure_simulations/satrapa/",
                    tools::file_path_sans_ext(i),
                    " --output=/media/burke/bigMac/Dropbox/structure_simulations/satrapa/fstr_out/",
                    outname,
                    " --seed=",sample(1:1e6,1))
  commands <- append(commands,command)
  }
}

#run in parallel
foreach(i=commands) %dopar% system(i)

