#simple population growth viz to figure out Nm levels for diff growth rates
library(tidyr);library(tools);library(reshape);library(plyr);library(ggplot2);library(magrittr)

growth_rate <- c(1e-3,1e-4,1e-5,0)
migration_rate <- 5e-5
popsize <- c();nmigrants <- c();eqFst <- c()
df <- data.frame(generation=numeric(),popsize=numeric(),nmigrants=numeric(),eqFst=numeric(),growth=numeric())
for(j in growth_rates){
  for(i in 1:40000){
    if(i==1){
      n_prev <- 500000
    } else if(i>1){
      n_prev <- n_next
    }
    n_next <- (n_prev*j)+n_prev
    popsize[i] <- n_next
    nmigrants[i] <- n_next*migration_rate
    eqFst[i] <- 1/(4*n_next*migration_rate+1)
  }
  tmp <- data.frame(generation=1:40000,popsize,nmigrants,eqFst,growth=j)
  df <- rbind(df,tmp)
}
meltdf <- melt(df,id.vars=c("generation","growth"))
ggplot(subset(meltdf,variable=="eqFst"),aes(x=generation,y=value))+ #slow : (
  theme_minimal()+
  facet_grid(growth~variable,scales="free_y")+
  geom_line()

