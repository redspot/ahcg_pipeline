


library(ggplot2)
args<-commandArgs(TRUE)


cov <- read.delim(args[1], header=FALSE)
names(cov) <- c("Chr", "Position", "Coverage")

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


