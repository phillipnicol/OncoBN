

library(ggplot2)
library(grDevices)

## A. CBN to DP comparison


##B.

dpmcc <- read.csv("/Users/phillipnicol/Desktop/dpmcc_ga.csv")
dpmcc <- dpmcc[,-1]
caprimcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/caprimcc.csv")
caprimcc <- caprimcc[,-1]
df <- as.data.frame(cbind(eta, apply(dpmcc,2,median), apply(dpmcc, 2, function(x) quantile(x,0.25)),
                          apply(dpmcc, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("eta", "Recall", "firstQ", "thirdQ")

df2 <- as.data.frame(cbind(eta, apply(caprimcc,2,median), apply(caprimcc, 2, function(x) quantile(x,0.25)),
                           apply(caprimcc, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("eta", "Recall", "firstQ", "thirdQ")

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,20)
supp[1:10] <- "GA"
supp[11:20] <- "CAPRI"
df$Algorithm <- supp

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=2.5) + geom_line(size = 1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("maroon","skyblue"))
p <- p + xlab("Noise level"~(eta))
p <- p + ylab("MCC")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)

