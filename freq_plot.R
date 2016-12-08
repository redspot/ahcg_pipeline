#!/usr/bin/env Rscript

#Rscript -e 'install.packages("ggplot2", dep = TRUE, repos="http://cran.rstudio.com/")'
library(ggplot2)

args<-commandArgs(TRUE)

fn=args[1]

cov <- read.delim(fn, header=FALSE)
names(cov) <- c("Chr", "Position", "Coverage")

options(bitmapType='cairo')
png(filename=args[2])
qplot(cov$Coverage,
      geom="histogram",
      binwidth=10,
      main="Distribution of Coverage",
      xlab="Coverage",
      fill=I("blue"),
      col=I("black"),
      alpha=I(.2),
      breaks=seq(0,350,5))

junk <- dev.off()
