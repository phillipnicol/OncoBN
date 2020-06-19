
validateCPN <- function(df, model = "CBN", algorithm = "DP", epsilon = 0.025, summaryPlot = TRUE) {
  validation <- c(1:ceiling(ncol(df)/2))
  for(i in validation) {
    out <- fitCPN(df, model = model, algorithm = algorithm, k = i, epsilon = epsilon, verbose = FALSE)
    validation[i] <- out$score
  }
  if(summaryPlot) {
    plot(validation, main = "Scree plot", xlab = "in-degree bound", ylab = "Network score", pch = 16)
  }
  out <- list()
  out$validation <- validation
  return(out)
}

cvalidateCPN <- function(df, nfolds = 5, model = "CBN", algorithm = "DP") {
  df < - df[sample(nrow(df)),]
  folds <- cut(seq(1,nrow(ColorectalCancer)),breaks=5,labels=FALSE)
  
  eps_seq <- seq(exp(-100), min(colMeans(df))/2, length.out = 12)

  CV <- matrix(0, nrow = length(eps_seq), ncol = ceiling(ncol(df)/2))
  print(dim(CV))
  
  for(i in 1:nfolds) {
    test_indices <- which(folds == i)
    train_df <- df[-test_indices,]
    test_df <- df[test_indices,]
    
    for(j in 1:length(eps_seq)) {
      for(k in 1:4) {
        out.fit <- fitCPN(train_df, model = model, algorithm = algorithm, epsilon = eps_seq[j], k = k, verbose = FALSE)
        out.graph <- inferTheta(train_df, out.fit)
        CV[j,k] = CV[j,k] + prediction(test_df, out.graph$parents, out.graph$thetas, eps_seq[j])
      }
    }
  }
  
  out <- list()
  rownames(CV) <- eps_seq
  out$cv.matrix <- CV/nfolds
  return(out)
}

prediction <- function(df, parents, thetas, epsilon) {
  df <- cbind(df, 1)
  colnames(df)[ncol(df)] <- "WT" 
  
  correct <- 0
  
  for(i in 1:(ncol(df) - 1)) {
    for(j in 1:nrow(df)) {
      if(any(df[j,parents[[i]] ] == 1)) {
        prediction <- rbinom(1,1,thetas[[i]])
      }
      else {
        prediction <- rbinom(1,1,epsilon)
      }
      
      if(prediction == df[j,i]) {
        correct <- correct + 1
      }
    }
  }
  return(correct/(nrow(df)*(ncol(df)-1)))
}