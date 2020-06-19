
bootstrapCPN <- function(df, model = "CBN", algorithm = "DP", k = 3, epsilon = 0.025, N = 100) {
  Edges <- matrix(0, nrow = 0, ncol = 2)
  freqs <- c()
  scores <- rep(0, N)
  for(i in 1:N) {
    bootstrap_df <- df[sample(1:nrow(df), size = nrow(df), replace = TRUE),]
    out <- fitCPN(bootstrap_df, model, algorithm, k, epsilon, verbose = FALSE)
    scores[i] <- out$score
    edge_pairs <- split(out$edgelist, as.integer((seq_along(out$edgelist) - 1) / 2))
    for(j in edge_pairs) {
      index <- which(apply(Edges, 1, function(x) all(x %in% j)))
      if(length(index) == 0) {Edges <- rbind(Edges, j); freqs <- c(freqs, 1)}
      else {freqs[index] <- freqs[index] + 1}
    }
  }
  out <- list()
  rownames(Edges) <- NULL
  Edges <- Edges[order(freqs, decreasing = TRUE),]
  freqs <- freqs[order(freqs, decreasing = TRUE)]
    
  out$Edges <- Edges
  out$counts <- freqs
  out$scores <- scores
  out$N <- N
  return(out)
}