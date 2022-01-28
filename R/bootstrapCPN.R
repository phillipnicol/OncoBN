#' @export
#' 
#' @title Bootstrap cancer progression networks 
#' 
#' @description This function performs a nonparameteric
#' bootstrap to quanity the uncertainty in the 
#' maximum likelihood network. 
#' 
#' @param df A dataframe or matrix of binary mutation data. The columns should be 
#' mutations and the rows should be patients. 
#' @param model Whether to fit a network with "conjuctive" interactions (CBN), or
#' "disjunctive" interactions (DBN). 
#' @param algorithm The algorithm used to infer the network. Default is
#' an optimal dynamic programming ("DP") algorithm, which is only feasible when the number of mutations is
#' smaller than 30. For larger datasets, use the genetic algorithm ("GA"), which 
#' gives an approximate solution. 
#' @param k The in-degree bound on the estimated network. 
#' @param epsilon Penalty term (for mutations not conforming to the estimated network).
#' @param N The number of bootstrap iterations.
#' 
#' @return A list with components 
#' \itemize{
#' \item \code{edgelist} - A vector of length 2*k listing
#' the vertices in the k inferred edges.
#' \item \code{score} - The maximum likelihood value.
#' } 
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
#' @examples 
#' ## Load the example data
#' data("example")
#' 
#' ## Fit a conjunctive Bayesian network (CBN)
#' out <- fitCPN(example, epsilon=0.01)
#' 
#' ## Plot the graph
#' bs <- boots
#' 
#' @references 
#' Nicol, PB, Coombes, KR, Deaver, C, Chkrebtii, O, Paul, S, Toland, AE, Asiaee, A. 
#' Oncogenetic network estimation with disjunctive Bayesian networks. Comp. Sys. Onco. 2021; 1:e1027.
#'  https://doi.org/10.1002/cso.21027 
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