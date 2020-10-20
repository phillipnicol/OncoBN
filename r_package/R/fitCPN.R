
fitCPN <- function(df, model = "CBN", algorithm = "DP", k = 3, epsilon = 0.025, verbose = TRUE) {
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
    results <- GA(df, leaky = epsilon, model = model, suppress = !verbose) 
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

