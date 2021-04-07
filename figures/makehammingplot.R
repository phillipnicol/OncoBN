

library(ggplot2)
library(grDevices)

##A. 

dpmcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/dp_cbnHAMMING.csv")
cbnmcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/cbnHAMMING.csv")
dpmcc <- dpmcc[,-1]
dpmcc <- 1-dpmcc/400
cbnmcc <- cbnmcc[,-1]
cbnmcc[is.na(cbnmcc)] <- 0
cbnmcc <- 1-cbnmcc/400
eta <- c(0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2)
df <- as.data.frame(cbind(eta, apply(dpmcc,2,median), apply(dpmcc, 2, function(x) quantile(x,0.25)),
                          apply(dpmcc, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("eta", "Recall", "firstQ", "thirdQ")

df2 <- as.data.frame(cbind(eta, apply(cbnmcc,2,median), apply(cbnmcc, 2, function(x) quantile(x,0.25)),
                           apply(cbnmcc, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("eta", "Recall", "firstQ", "thirdQ")

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,20)
supp[1:10] <- "DP"
supp[11:20] <- "MC-CBN"
df$Algorithm <- supp

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=2.5) + geom_line(size = 1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("sienna2","darkblue"))
p <- p + xlab("Noise level"~(eta))
p <- p + ylab("Accuracy")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)
ggsave(filename="/Users/phillipnicol/Desktop/hd_1.png",width=7,height=3)


##B.

dpmcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/dpmccHAMMING.csv")
dpmcc <- dpmcc[,-1]
dpmcc <- 1-dpmcc/400
caprimcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/caprimccHAMMING.csv")
caprimcc <- caprimcc[,-1]
caprimcc <- 1-caprimcc/400
df <- as.data.frame(cbind(eta, apply(dpmcc,2,median), apply(dpmcc, 2, function(x) quantile(x,0.25)),
                          apply(dpmcc, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("eta", "Recall", "firstQ", "thirdQ")

df2 <- as.data.frame(cbind(eta, apply(caprimcc,2,median), apply(caprimcc, 2, function(x) quantile(x,0.25)),
                           apply(caprimcc, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("eta", "Recall", "firstQ", "thirdQ")

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,20)
supp[1:10] <- "DP"
supp[11:20] <- "CAPRI"
df$Algorithm <- supp

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=2.5) + geom_line(size = 1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("maroon","sienna2"))
p <- p + xlab("Noise level"~(eta))
p <- p + ylab("Accuracy")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
ggsave(filename="/Users/phillipnicol/Desktop/hd_2.png",width=7,height=3)
print(p)
