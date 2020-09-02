
library(TRONCO) 


randomNoiseTest <- function() {
  edges <- c()
  tests <- c(10,15,20,25,30,40,50,60,70,80,90,100)
  for(test in tests) {
    D <- matrix(rbinom(400*test, 1, 0.2), nrow = 400, ncol = test)
    colnames(D) <- c(1:test)
    Dt <- import.genotypes(D)
    fit <- tronco.capri(Dt)
    edges <- c(edges, sum(fit$model$capri_bic$adj.matrix$adj.matrix.fit))
    
  }
  return(edges) 
}

M <- matrix(0, nrow = 10, ncol = 12)
for(i in 1:10) {
  M[i,] <- randomNoiseTest()
}


library(ggplot2) 

df <- as.data.frame(cbind(c(10,15,20,25,30,40,50,60,70,80,90,100), colMeans(M)))
colnames(df) <- c("n", "edges")


p <- ggplot(df,aes(x=n,y=edges)) + geom_point(size=2, color = "blue")
p <- p + ylab("(Mean) edges")
p <- p + theme_linedraw()
p <- p + xlab("n")
p <- p + theme(legend.title = element_text(size=20), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
               axis.text.x = element_text(size=15), axis.text.y = element_text(size=15),
               legend.text=element_text(size=15))
print(p)


