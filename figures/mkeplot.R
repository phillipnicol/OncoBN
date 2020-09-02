

library(ggplot2)
library(grDevices)

df <- as.data.frame(cbind(eta, apply(dpprec,2,median), apply(dpprec, 2, function(x) quantile(x,0.25)), 
                          apply(dpprec, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("eta", "Precision", "firstQ", "thirdQ")

df2 <- as.data.frame(cbind(eta, apply(cbnprec,2,median), apply(cbnprec, 2, function(x) quantile(x,0.25)), 
                           apply(cbnprec, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("eta", "Precision", "firstQ", "thirdQ")

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,20)
supp[1:10] <- "SM"
supp[11:20] <- "MC-CBN"
df$Algorithm <- supp


p <- ggplot(df, aes(x =eta, y = Precision, color=Algorithm)) + geom_point(size=1.75) + geom_line() 
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_linedraw()
p <- p + scale_color_manual(values=c("sienna2","maroon"))
p <- p + ylab("Precision")
p <- p + xlab("")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)




df <- as.data.frame(cbind(eta, apply(dprecall,2,median), apply(dprecall, 2, function(x) quantile(x,0.25)), 
                          apply(dprecall, 2, function(x) quantile(x,0.75))))
colnames(df) <- c("eta", "Recall", "firstQ", "thirdQ")

df2 <- as.data.frame(cbind(eta, apply(cbnrecall,2,median), apply(cbnrecall, 2, function(x) quantile(x,0.25)), 
                           apply(cbnrecall, 2, function(x) quantile(x,0.75))))
colnames(df2) <- c("eta", "Recall", "firstQ", "thirdQ")

df <- as.data.frame(rbind(df,df2))
supp <- rep(0,20)
supp[1:10] <- "SM"
supp[11:20] <- "MC-CBN"
df$Algorithm <- supp

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=1.75) + geom_line() 
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_linedraw()
p <- p + scale_color_manual(values=c("sienna2","maroon"))
p <- p + xlab("Noise level"~(eta))
p <- p + ylab("Recall")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)










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

p <- ggplot(df, aes(x =eta, y = Recall, color=Algorithm)) + geom_point(size=1.75) + geom_line() 
p <- p + geom_errorbar(aes(ymin=firstQ, ymax=thirdQ), width = 0.005)
p <- p + theme_linedraw()
p <- p + scale_color_manual(values=c("sienna2","maroon"))
p <- p + xlab("Noise level")
p <- p + ylab("Recall")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)




