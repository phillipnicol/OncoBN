#' @title Fit a Bayesian network model to cancer mutation data 
#' 
#' @description Fit a Bayesian network model to cancer evolution data. 
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
#' Default 
#' @param ngen For genetic algorithm only: number of generations. 
#' @param popsize For genetic algorithm only: initial population size. 
#' @param verbose Whether or not to print to the console. 
#' 
#' @return A list with components 
#' \itemize{
#' \item \code{edgelist} - A vector of length 2*k listing
#' the vertices in the k inferred edges.
#' \item \code{score} - The maximum likelihood value.
#' \item \code{graph} - The maximum likelihood graph saved as 
#' an igraph object.
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
#' plotCPN(out)
#' 
#' @references 
#' Nicol, PB, Coombes, KR, Deaver, C, Chkrebtii, O, Paul, S, Toland, AE, Asiaee, A. 
#' Oncogenetic network estimation with disjunctive Bayesian networks. Comp. Sys. Onco. 2021; 1:e1027.
#'  https://doi.org/10.1002/cso.21027 
fitCPN <- function(df, model = "CBN", algorithm = "DP", k = 3, 
                   epsilon = colMeans(df)/2, ngen=100, popsize=100, verbose = TRUE) {
  input <- list()
  input$df <- t(data.matrix(df))
  input$dims <- dim(df)
  input$options <- c(k, epsilon)
  input$model <- model
  input$verbose <- verbose
  
  if(!any(model == "CBN", model == "DBN")) {
    stop("Invalid model selected.")
  }
  
  #call Rcpp function
  if(algorithm == "DP") {
    results <- dp(input)
    #add wildtype to dataframe
    df <- cbind(df, 1)
    colnames(df)[ncol(df)] <- "WT" 
  }
  else if(algorithm == "GA") {
    results <- GA(df, leaky = epsilon, MAX_GEN=ngen, N=popsize,
                  model = model, suppress = !verbose) 
    print(results[[1]])
    #add wildtype to dataframe
    df <- cbind(1, df)
    colnames(df)[1] <- "WT" 
  }
  else {
    stop("Invalid algorithm selected.")
  }
  
  #format return list
  out <- list()
  out$edgelist <- sapply(results[[1]], function(x) colnames(df)[x])
  out$score <- results[[2]]
  out$graph <- igraph::make_directed_graph(out$edgelist)
  return(out)
}

inferTheta <- function(df, fit) {
  df <- cbind(df, 1)
  colnames(df)[ncol(df)] <- "WT" 
  
  g <- make_directed_graph(fit$edgelist) 
  parent_list <- list()
  theta_list <- list()
  for(i in 1:(ncol(df) - 1)) {
    parents <- incident(g, colnames(df)[i], mode = "in")
    parent_list[[i]] <- ends(g, parents)[,1]

    p_on <- as.data.frame(df[, parent_list[[i]] ])
    theta_list[[i]] <- length(which(df[rowSums(p_on) > 0, i] == 1))/nrow(as.data.frame(p_on[rowSums(p_on) > 0, ]))
    if(theta_list[[i]] > 1) {
      stop(theta_list[[i]])
    }
  }
  out <- list()
  out$parents <- parent_list
  out$thetas <- theta_list
  return(out)
}

Likelihood <- function(df, graph.fit, epsilon) {
  df <- cbind(df, 1)
  colnames(df)[ncol(df)] <- "WT" 
  likelihood <- 0
  
  for(i in 1:(ncol(df) - 1)) {
    p_on <- as.data.frame(df[,graph.fit$parents[[i]] ])
    
    c_and_p <- length(which(df[rowSums(p_on) > 0, i] == 1))
    no_c_and_p <- nrow(df[rowSums(p_on) > 0, ]) - c_and_p
    
    c_and_no_p <- nrow(df)*colMeans(df)[i] - c_and_p
    
    theta_i <- graph.fit$thetas[[i]]
    if(theta_i == 1) {
      likelihood = likelihood + c_and_no_p*log(epsilon)      
    }
    else{
      likelihood = likelihood + c_and_p*log(theta_i) + no_c_and_p*log(1-theta_i)
      likelihood = likelihood + c_and_no_p*log(epsilon)      
    }
  }
  return(likelihood)
}

