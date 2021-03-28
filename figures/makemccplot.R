

library(ggplot2)
library(grDevices)

## A. CBN to DP comparison


library(ggplot2)
dpmcc <- readRDS("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/dpmcc.RDS")
cbnmcc <- readRDS("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/cbnmcc.RDS")
cbnmcc[is.na(cbnmcc)] <- 0
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

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=2.5) + ylim(-0.5,1) + geom_line(size = 1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("sienna2","darkblue"))
p <- p + xlab("Noise level"~(eta))
p <- p + ylab("MCC")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)

##B.

dpmcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/dpmcc.csv")
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
supp[1:10] <- "DP"
supp[11:20] <- "CAPRI"
df$Algorithm <- supp

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=2.5) + geom_line(size = 1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("maroon","sienna2"))
p <- p + xlab("Noise level"~(eta))
p <- p + ylab("MCC")
p <- p + ylim(0,1)
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)


dpmcc_NO_FN <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/dpmcc2_NO_FN.csv")
dpmcc_NO_FN <- dpmcc_NO_FN[,-1]
caprimcc_NO_FN <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/caprimcc_NO_FN.csv")
caprimcc_NO_FN <- caprimcc_NO_FN[,-1]
df <- as.data.frame(cbind(eta, apply(dpmcc_NO_FN,2,median), apply(dpmcc_NO_FN, 2, function(x) quantile(x,0.25)),
                          apply(dpmcc_NO_FN, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("eta", "Recall", "firstQ", "thirdQ")

df2 <- as.data.frame(cbind(eta, apply(caprimcc_NO_FN,2,median), apply(caprimcc_NO_FN, 2, function(x) quantile(x,0.25)),
                           apply(caprimcc_NO_FN, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("eta", "Recall", "firstQ", "thirdQ")

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,20)
supp[1:10] <- "DP"
supp[11:20] <- "CAPRI"
df$Algorithm <- supp

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=2.5) + geom_line(size = 1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("maroon", "sienna2"))
p <- p + xlab("False Positive Rate")
p <- p + ylab("MCC")
p <- p + ylim(0,1)
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)


## C. Increasing net. complexity
dpmcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/dpmccIC.csv")
dpmcc <- dpmcc[,-1]
caprimcc <- read.csv("/Users/phillipnicol/Desktop/local_files/data/CPN_DP_DATA/troncoIC.csv")
caprimcc <- caprimcc[,-1]
df <- as.data.frame(cbind(c(1:8), apply(dpmcc,2,median), apply(dpmcc, 2, function(x) quantile(x,0.25)),
                          apply(dpmcc, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("avg_deg", "Recall", "firstQ", "thirdQ")
df$algorithm <- "DP"

df2 <- as.data.frame(cbind(c(1:8), apply(caprimcc,2,median), apply(caprimcc, 2, function(x) quantile(x,0.25)),
                           apply(caprimcc, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("avg_deg", "Recall", "firstQ", "thirdQ")
df2$algorithm <- "CAPRI"

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,16)
supp[1:8] <- "DP"
supp[9:16] <- "CAPRI"
df$Algorithm <- supp

p <- ggplot(df, aes(x = avg_deg, y = Recall, color = Algorithm)) + geom_point(size=2.5) + geom_line(size=1.25)
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_classic()
p <- p + scale_color_manual(values=c("maroon", "sienna2"))
p <- p + xlab("Average degree")
p <- p + ylab("MCC")
p <- p + ylim(0,1)
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)



