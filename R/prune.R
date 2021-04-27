

prune <- function(D, head, tail, pval) {
  c11 <- nrow(D[D[,head] == 1 & D[,tail] == 1,])
  c12 <- nrow(D[D[,head] == 1 & D[,tail] == 0,])
  c21 <- nrow(D[D[,head] == 0 & D[,tail] == 1,])
  c22 <- nrow(D[D[,head] == 0 & D[,tail] == 0,])
  
  C <- matrix(c(c11,c21,c12,c22), nrow = 2, ncol = 2, byrow = TRUE)
  return(fisher.test(C, alternative = "greater")$p.value < pval)
}