## perform anova / tukey test on kmeans / dapc results

setwd("~/Dropbox/structure_simulations/")
library(ggplot2);library(gtable);library(gridExtra);library(grid);library(dplyr);library(reshape)

# simulated data

# read in structure output 
sim_str <- as.data.frame(read.csv("sim_str_out.csv"))
sim_str$maf <- as.factor(sim_str$maf)
aov_str_sim <- aov(pop_clust_dist_scaled ~ maf, sim_str)
tukeysimstr <- TukeyHSD(aov_str_sim)
a <- as.data.frame(tukeysimstr$maf)
a <- abs(a$diff)
b <- rep("structure", length(a))
df_str <- cbind.data.frame(a,b)

# read in kmeans output
sim_km  <- read.csv("clust_multivar_out.csv")

# create a new column to identify results by MAF level
files <- as.character(sim_km$file)
maf <- lapply(files, function(i){
  gsub(pattern = "(.*maf)(.*)(.str)", replacement = "\\2", x = i)
})
sim_km$maf <- as.factor(unlist(maf))

aov_km_sim <- aov(kmeans ~ maf, sim_km)
tukeysimkm <- TukeyHSD(aov_km_sim)
a <- as.data.frame(tukeysimkm$maf)
a <- abs(a$diff)
b <- as.factor(rep("kmeans", length(c)))
x <- cbind.data.frame(a,b)

aov_xval_sim <- aov(xval ~ maf, sim_km)
tukeysimxv <- TukeyHSD(aov_xval_sim)
a <- as.data.frame(tukeysimxv$maf)
a  <- abs(a$diff)
b <- as.factor(rep("DAPC", length(c)))
y <- cbind.data.frame(a,b)
df_pc <- rbind(x,y)
df_sim <- rbind(df_str,df_pc)
colnames(df_sim) <- c("difference","method")

# empirical data
reg_str <- as.data.frame(read.csv("satrapa_str_out.csv"))
reg_str$maf <- as.factor(reg_str$maf)
aov_str_reg <- aov(pop_clust_dist_scaled ~ maf, reg_str)
tukeyregstr <- TukeyHSD(aov_str_reg)
c <- as.data.frame(tukeyregstr$maf)
c <- abs(c$diff)
d <- rep("structure", length(c))
df2_str <- cbind.data.frame(c,d)

reg_km  <- read.csv("./satrapa/satrapa_cluster_multivar_out.csv")
reg_km$maf <- as.factor(reg_km$maf)
aov_km_reg <- aov(kmeans ~ maf, reg_km)
tukeyregkm <- TukeyHSD(aov_km_reg)
c <- as.data.frame(tukeyregkm$maf)
c <- abs(c$diff)
d <- as.factor(rep("kmeans", length(c)))
df2_km <- cbind.data.frame(c,d)

aov_xval_reg <- aov(xval ~ maf, reg_km)
tukeyregxv <- TukeyHSD(aov_xval_reg)
c <- as.data.frame(tukeyregxv$maf)
c <- abs(c$diff)
d <- as.factor(rep("DAPC", length(c)))
df2_xv <- cbind.data.frame(c,d)
df_reg <- rbind(df2_str,df2_km,df2_xv)
colnames(df_reg) <- c("difference","method")

# tukey test on sim difference of means across MAF filters between methods
sstr <- as.data.frame(tukeysimstr$maf)
sstr <- abs(sstr$diff)
skm <- as.data.frame(tukeysimkm$maf)
skm <- abs(skm$diff)
sxv <- as.data.frame(tukeysimxv$maf)
sxv <- abs(sxv$diff)
r1 <- as.data.frame(cbind(sstr, skm, sxv))
r1 <- melt(r1, value.name="diff",varnames=c('sstr','skm','sxv'))
r1aov <- aov(value~variable, r1)
TukeyHSD(r1aov)

# tukey test on empirical difference of means across MAF filters between methods
rstr <- as.data.frame(tukeyregstr$maf)
rstr <- abs(rstr$diff)
rkm <- as.data.frame(tukeyregkm$maf)
rkm <- abs(rkm$diff)
rxv <- as.data.frame(tukeyregxv$maf)
rxv <- abs(rxv$diff)
r2 <- as.data.frame(cbind(rstr, rkm, rxv))
r2 <- melt(r2, value.name="diff",varnames=c('rstr','rkm','rxv'))
r2aov <- aov(value~variable, r2)
TukeyHSD(r2aov)

# plot!
p1 <- ggplot(data = df_sim, aes(x = method, y = difference)) +
  geom_boxplot(aes(fill = method)) +
  theme_bw()+
  guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=1)+
  theme(panel.grid=element_blank(),axis.ticks = element_blank())+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank()) +
  ylab("Mean Tukey Difference") +
  annotate(geom = "text", y = max(df_sim$difference), x = 3.5, label="Simulated data",hjust=1)

p2 <- ggplot(data = df_reg, aes(x = method, y = difference)) +
  geom_boxplot(aes(fill = method)) +
  theme_bw()+
  guides(fill=FALSE)+
  scale_fill_grey(start=0.5, end=1)+
  theme(panel.grid=element_blank(),axis.ticks = element_blank())+
  xlab("Clustering method")+
  ylab("Mean Tukey Difference")+
  annotate(geom = "text", y = max(df_reg$difference), x = 3.5, label="Empirical data",hjust=1)

grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))



