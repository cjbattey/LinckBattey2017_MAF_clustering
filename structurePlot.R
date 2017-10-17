########## structurePlot: an R function for plotting Structure output
##Required input: path to structure output file (ending in "_f")
##                pop.info = TRUE/FALSE for presence of column of putative pop info in structure input. Defaults to TRUE.
##
##optional input: colors = vector of at least k colors. Required if k > 8.
##                pop.order = TRUE plots in order of populations as in structure input file optional column "putative pop info" 
##                spacing = Amount of white space left between bars, as a fraction of bar width. Defaults to 0. 
##                outline = color used for outline of bars (default is white).
##                cex=font scaling for labels. defaults to 0.75.
##                names=vector of sample ID's. Defaults to input. For unlabeled plots, set names=NULL. Do not use anything other than default or NULL. 
##
##experimental:   color.matching = TRUE will attach colors to putative populations (requires pop.info=TRUE) for consistent 
##plotting across multiple runs. The cluster attached to a given population is taken to be that with the highest 
##average ancestry across all individuals in that population. 
##If color.matching=TRUE, user-specified colors are assigned by population number - first color in the vector for pop 1, second for pop 2, etc.
##
##To make plots that look like the standard STRUCTURE plots from the GUI, set spacing=0 and outline='black'
#good colors: 
#c("gold","forestgreen","magenta3","orangered","cornflowerblue","orange","sienna","dodgerblue4")
#c(colors=brewer.pal(12,"Set3"))

structurePlot <- function(strOutput,pop.order=FALSE,q.order=FALSE,pop.info=TRUE,
                          colors=c("gold","forestgreen","magenta3","orangered","cornflowerblue","orange","sienna","dodgerblue4"),
                          color.matching=FALSE,spacing=0,outline='white',cex=0.75,names=df$sample) {
  require(plyr);require(reshape);require(RColorBrewer)
  infile <- readLines(strOutput,warn=FALSE)
  a <- grep("Inferred ancestry of individuals:",infile)+2
  b <- grep("Estimated Ln Prob of Data",infile)
  c <- grep("Estimated Allele Frequencies in each cluster",infile)-3
  k <- as.numeric(substr(infile[grep("populations assumed",infile)],4,5))
  lnml <- infile[b]
  df <- read.table(strOutput,skip=a-1,nrow=(c-a+1))
  df <- df[2:ncol(df)]
  if (pop.info=="TRUE"){
    colnames(df) <- c("sample","percent.missing","X",1:k)
    df <- df[!is.na(colnames(df))] #if 90% intervals are included
    n.pops <- max(as.numeric(as.character(df$pop)))
    if (color.matching=="TRUE"){
      palette <- c(rep("NA",k))
      for(i in 1:n.pops){
        d <- subset(df,pop==i)
        e <- colMeans(d[5:ncol(d)])
        f <- as.numeric(names(e[which(e==max(e))]))
        if (palette[f] == "NA"){
          palette[f] <- colors[i]
        }
      }
      palette[which(palette == "NA")] <- colors[which(colors%in%palette == FALSE)]
      colors <- palette
    }
    if (pop.order=="TRUE"){
      df <- arrange(df,pop)
      if(q.order==T){
        df <- arrange(df,pop,`1`)
      }
      a <- df[5:ncol(df)]
      barplot(t(a),axes=FALSE,col=colors,border=outline,names=names,cex.names=cex,las=2,cex.main=0.75,font.main=1,space=spacing,xpd=FALSE)
    }
    if(q.order==T){
      df <- arrange(df,pop,`1`)
    }
    a <- df[5:ncol(df)]
    barplot(t(a),axes=FALSE,col=colors,border=outline,names=names,cex.names=cex,las=2,cex.main=0.75,font.main=1,space=spacing,xpd=FALSE)
  } #end pop.info=T
  
  else {
    colnames(df) <- c("sample","percent.missing","X",1:k)
    if(q.order==T){
      df <- arrange(df,`1`)
    }
    a <- df[4:ncol(df)]
    n.pops <- k
    barplot(t(a),axes=FALSE,col=colors,border=outline,names=names,cex.names=cex,las=2,cex.main=0.75,font.main=1,space=spacing,xpd=FALSE)
  }
} 


