#structure output combo & analysis
setwd("~/Dropbox/structure_simulations/")
source("./R_scripts/ggthemes.R");source("./R_scripts/structurePlot.R")
library(ggplot2);library(plyr);library(magrittr);library(data.table);library(ggjoy);library(gridExtra)

#index of files by maf
files <- list.files("satrapa/str_out/k3")
maf_index <- list()
#maf_levels <- c(0.0,0.0125,0.025,0.0375,0.05,0.1,0.25)
maf_levels <- c(1,2,3,4,5,8,20)
b <- 1
for(i in maf_levels){
  maf_index[[b]] <- grep(paste0("maf",i,"."),files,fixed=T,value=T)
  b <- b+1
}

#get summary statistics from structure files (simulated data)
df <- data.frame(matrix(nrow=length(maf_index)*length(maf_levels),ncol=13,data=1))
q <- list()
row <- 1
for(i in 1:length(maf_index)){
  maf_level <- maf_index[[i]]
  for(j in 1:length(maf_level)){
    tmp <- maf_level[j]
    str <- readLines(paste0("./str_out/k3/",tmp))
    alpha <- str[grep("alpha",str)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
    lnl <- str[grep("Mean value of ln likelihood",str)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
    v_lnl <- str[grep("Variance of ln likelihood",str)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
    popclustdist <- str[(grep("Inferred Clusters",str)+3):(grep("Inferred Clusters",str)+5)] %>% 
                      strsplit(.," * ") %>% lapply(function(e) e[3:5]) %>% unlist() %>% as.numeric()
    q_matrix <- str[(grep("Label",str)+1):(grep("Label",str)+30)] %>% 
                  strsplit(" * ") %>% 
                    lapply(function(e) c(id=e[3],pop=e[5],q1=e[7],q2=e[8],q3=e[9])) %>% 
                      do.call(rbind.data.frame,.)
    colnames(q_matrix) <- c("id","pop","q1","q2","q3")
    q_matrix$q1 <- as.numeric(as.character(q_matrix$q1))
    q_matrix$q2 <- as.numeric(as.character(q_matrix$q2))
    q_matrix$q3 <- as.numeric(as.character(q_matrix$q3))
    q[[row]] <- q_matrix
    df[row,] <- c(maf_levels[i],alpha,lnl,v_lnl,popclustdist)
    row <- row+1
  }
}
names(df) <- c("maf","alpha","lnl","v_lnl","pop1clust1","pop1clust2","pop1clust3","pop2clust1","pop2clust2","pop2clust3","pop3clust1","pop3clust2","pop3clust3")
names(q) <- unlist(maf_index)
df <- read.csv("~/Dropbox/structure_simulations/sim_str_out.csv")

#get summary statistics from structure files (empirical data)
df <- data.frame(matrix(nrow=length(maf_index)*length(maf_levels),ncol=4,data=1))
q <- list()
row <- 1
for(i in 1:length(maf_index)){
  maf_level <- maf_index[[i]]
  for(j in 1:length(maf_level)){
    tmp <- maf_level[j]
    str <- readLines(paste0("satrapa/str_out/k3/",tmp))
    alpha <- str[grep("alpha",str)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
    lnl <- str[grep("Mean value of ln likelihood",str)] %>% strsplit(" * ") %>% unlist() %>% .[7] %>% as.numeric()
    v_lnl <- str[grep("Variance of ln likelihood",str)] %>% strsplit(" * ") %>% unlist() %>% .[6] %>% as.numeric()
    q_matrix <- str[(grep("Label",str)+1):(grep("Label",str)+33)] %>% 
      strsplit(" * ") %>% 
        lapply(function(e) c(id=e[3],q1=e[6],q2=e[7],q3=e[8])) %>% 
          do.call(rbind.data.frame,.)
    colnames(q_matrix) <- c("id","q1","q2","q3")
    q_matrix$q1 <- as.numeric(as.character(q_matrix$q1))
    q_matrix$q2 <- as.numeric(as.character(q_matrix$q2))
    q_matrix$q3 <- as.numeric(as.character(q_matrix$q3))
    q_matrix$pop <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,2,1,1,1)
    q[[row]] <- q_matrix
    df[row,] <- c(maf_levels[i],alpha,lnl,v_lnl)
    row <- row+1
  }
}
names(df) <- c("maf","alpha","lnl","v_lnl")
names(q) <- unlist(maf_index)
df

#find distance between simulated populations in q-matrix space
#need a list of df's each with 1 row per population
mean_q_list <- lapply(q,function(e) {
  ddply(e,"pop",summarize,q1=mean(q1),q2=mean(q2),q3=mean(q3))[-1]
})
q_dist_list <- sapply(mean_q_list,function(e){
  mean(dist(e))
})
df$pop_clust_dist <- q_dist_list
#data.frame(q1=c(1,0,0),q2=c(0,1,0),q3=c(0,0,1)) %>% dist() %>% mean() #maximum possible distance = 1.414214
df$pop_clust_dist_scaled <- df$pop_clust_dist/1.414214
  
#plot population discrimination
df$log_alpha <- log(df$alpha)
meltdf <- melt(df,id.vars="maf") %>% subset(variable %in% c("pop_clust_dist","log_alpha"))
meltdf$maf <- factor(meltdf$maf)

# facet solution
barplot <- ggplot(meltdf,aes(x=maf,y=value))+theme_minimal()+facet_grid(~variable,scales="free")+
  xlab("Minimum Minor Allele Count")+ylab("Population Discrimination")+
  coord_flip()+geom_point(position = position_jitter(width=0.25),col="grey")+
  geom_boxplot(fill=NA,outlier.color = NA)

joyplot <- ggplot(data=meltdf,aes(x=value,y=maf))+theme_minimal()+facet_grid(~variable,scales="free")+
  #facet_grid(~variable,scales="free")+
  ylab("Minimum Minor Allele Count")+
  xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy()

# read in and plot empirical data
df <- read.csv("sim_str_out.csv")
df$maf <- as.factor(df$maf)
df$log_alpha <- log(df$alpha)

jp_alpha <- ggplot(df, aes(x = log_alpha, y = maf))+
  theme_minimal()+
  ylab("Minimum Minor Allele Count")+
  xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy()

jp_popdis <- ggplot(df, aes(x = pop_clust_dist, y = maf))+
  theme_minimal()+
  ylab("Minimum Minor Allele Count")+
  xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy()

# read in and plot empirical data
df2 <- read.csv("satrapa_str_out.csv")
df2$maf <- as.factor(df2$maf)
levels(df2$maf) <- c(0,1,2,3,4,8,20)
df2$log_alpha <- log(df2$alpha)

jp_alpha_2 <- ggplot(df2, aes(x = log_alpha, y = maf))+
  theme_minimal()+
  ylab("Minimum Minor Allele Count")+
  xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy()

jp_popdis_2 <- ggplot(df2, aes(x = pop_clust_dist, y = maf))+
  theme_minimal()+
  ylab("Minimum Minor Allele Count")+
  xlab("")+
  scale_y_discrete(expand = c(0.01, 0))+
  geom_joy()

grid.arrange(jp_alpha, jp_popdis, jp_alpha_2, jp_popdis_2, ncol=2)

#structure plots (simulated)
par(mfrow=c(7,1),mar=c(0,0,0.5,0))
structurePlot("./str_out/k3/sim03.str_maf20.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
structurePlot("./str_out/k3/sim03.str_maf8.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
structurePlot("./str_out/k3/sim03.str_maf5.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
structurePlot("./str_out/k3/sim03.str_maf4.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
structurePlot("./str_out/k3/sim03.str_maf3.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
structurePlot("./str_out/k3/sim03.str_maf2.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
structurePlot("./str_out/k3/sim03.str_maf1.str_01_f",colors=grey.colors(3),color.matching=T,outline="black")
par(mfrow=c(1,1))

#structure plots (empirical satrapa) -- file paths got all fucked up by textedit
par(mfrow=c(7,1),mar=c(0,0,0.5,0))
structurePlot("./satrapa/str_out/plotting/maf0.25_f.str_01_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
structurePlot("./satrapa/str_out/plotting/maf0.1_f.str_01_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
structurePlot("./satrapa/str_out/plotting/maf0.05_f.5_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
structurePlot("./satrapa/str_out/plotting/maf0.0375_f.str_01_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
structurePlot("./satrapa/str_out/plotting/maf0.025_f.str_01_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
structurePlot("./satrapa/str_out/plotting/maf0.0125_f.str_01_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
structurePlot("./satrapa/str_out/plotting/maf0_f.str_01_f",colors=grey.colors(3),color.matching = T,pop.info=T,outline="black")
par(mfrow=c(1,1))

#ggplot2 / grid structure plots
#plot highest-likelihood structure output for each maf (no color matching)
max_lnl <- ddply(df,"maf",summarize,max_lnl=max(lnl))
best_runs <- df[df$lnl %in% max_lnl$max_lnl,] %>% rownames() %>% as.numeric()
best_q <- q[best_runs]
best_q <- do.call(rbind,best_q)
best_q$maf <- lapply(c(1,2,3,4,5,8,20),function(e) rep(e,30)) %>% unlist()

meltq <- melt(best_q,id.vars=c("maf","id","pop"))
meltq$maf <- factor(meltq$maf,levels=c(20,8,5,4,3,2,1))

structure_plots <- ggplot(meltq,aes(x=id,y=value,fill=variable))+facet_grid(maf~.)+
  theme_minimal()+ylab("")+xlab("")+
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_manual(values=grey.colors(3))+
  geom_bar(stat="identity")

#side-by-side plots
layout <- rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),
                c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),c(1,2),
                c(1,NA))
grid.arrange(joyplot,structure_plots,layout_matrix=layout)


