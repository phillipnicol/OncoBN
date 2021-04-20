
GA <- function(Data, N=100, MAX_GEN=100, p_mat_c=0.5, p_cx=0.5, p_ba=0.01, p_ea=0.05, p_perm_mut=0.01, leaky, model, suppress = FALSE) {
  penalty <- "hard"
  Data <- cbind(1,Data)
  colnames(Data)[1] <- "WT"
  
  n <- ncol(Data)
  num_samples <- nrow(Data)
  
  if(!suppress) {
    cat("Preparing Data... ...\n")
  }
  
  Data <- as.matrix(Data)
  Unique_Data <- matrix(0, nrow = 1, ncol = n)
  Unique_Data[1,] = Data[1,]
  counts <- 1
  ###Given Data, Find Unique Data Values. Extract counts
  for(i in 2:nrow(Data))
  {
    for(j in 1:nrow(Unique_Data))
    {
      if(identical(Data[i,], Unique_Data[j,]) == TRUE)
      {
        counts[j] = counts[j] + 1
        break
      }
      if(j == nrow(Unique_Data))
      {
        Unique_Data = rbind(Unique_Data, Data[i,])
        counts[j+1] = 1
      }
    }
  }
  
  mode(Unique_Data) <- "integer"
  
  Node_counts <- list()
  for(i in 2:n)
  {
    Node_counts[[i]] = which(Data[,i] == 1)
  }
  
  ##### INITIALIZE POPULATION ######
  
  if(!suppress)
  {
    cat("Initializing Population ... ... \n")   
  }
  
  ## Run CAPRI
  gen <- TRONCO::import.genotypes(Data)
  capri <- TRONCO::tronco.capri(gen)
  
  init_perm <- order(colSums(capri$adj.matrix.prima.facie))
  init_mat <- capri$adj.matrix.prima.facie[init_perm,init_perm]
  
  ### Initialize Binary matrix chromosomes
  Mats <- list()
  for(i in 1:N)
  {
    Mats[[i]] <- init_mat 
  }
  
  # Sort in decreasing order
  freqs <- colMeans(Data)
  initial_perm <- order(freqs, decreasing=TRUE)
  
  ### The other chromosome will be random permutations
  Perms <- list()
  for(i in 1:N)
  {
    Perms[[i]] = initial_perm
    Perms[[i]] = Perm_Repair(Mats[[i]], Perms[[i]])
  }
  
  ### For each Graph, infer theta
  ### Calculate fitness for each member of the population
  fitness <- c(1:N)
  for(i in 1:N)
  {
    p = invPerm(Perms[[i]])
    Temp = Mats[[i]][p,]
    theta = Infer_Theta(Temp[,p], Node_counts, counts)
    fitness[i] = Fitness(Mats[[i]], Perms[[i]], theta, leaky, Unique_Data, counts, model)
  }
  
  ###Save the maximum value of the graph
  key = which(fitness == max(fitness))
  Max_fitness = max(fitness)
  Max_perm = Perms[[key[1]]]
  Max_edge = Mats[[key[1]]]
  max_gen = 1
  
  ###Initialize Data Frame to store information
  Results <- matrix(0, nrow = MAX_GEN, ncol = 2)
  Results <- as.data.frame(Results)
  
  if(!suppress)
  {
    cat("Done! Starting GA ... ... \n")    
  }
  
  ##### Genetic Algorithm #####
  generation = 1
  while(generation <= MAX_GEN)
  {
    if(generation %% 10 == 0)
    {
      if(!suppress)
      {
        cat("Generation ", generation, " ... ... \n")
      }
    }
    
    ###Create copies of the lists we currently have
    Mats_c = Mats
    Perms_c = Perms
    fitness_c = fitness
    
    dart_board = Create_Dart_Board(fitness, N)
    
    ### Selection
    for(i in 1:(N/2))
    {
      selected = Selection(dart_board)
      
      Children = DAGrR(Mats[[selected[1]]], Mats[[selected[2]]], Perms[[selected[1]]],
                       Perms[[selected[2]]], fitness[selected[1]], fitness[selected[2]],
                       p_mat_c, p_cx, p_ba, p_ea, p_perm_mut, penalty)
      
      Mats_c[[2*i - 1]] = Children[[1]]
      Mats_c[[2*i]] = Children[[3]]
      Perms_c[[2*i - 1]] = Children[[2]]
      Perms_c[[2*i]] = Children[[4]]
      
      if(identical(Children[[1]], Mats[[selected[1]]]) & 
         identical(Children[[2]], Perms[[selected[1]]]))
      {
        fitness_c[2*i - 1] = fitness[selected[1]]
      }
      else
      {
        p = invPerm(Children[[2]])
        Temp = Children[[1]][p,]
        theta = Infer_Theta(Temp[,p], Node_counts, counts)
        fitness_c[2*i - 1] = Fitness(Children[[1]], Children[[2]], theta, leaky, Unique_Data,
                                       counts, model)
      }
      
      if(identical(Children[[3]], Mats[[selected[2]]]) & 
         identical(Children[[4]], Perms[[selected[2]]]))
      {
        fitness_c[2*i] = fitness[selected[2]]
      }
      else
      {
        p = invPerm(Children[[4]])
        Temp = Children[[3]][p,]
        theta = Infer_Theta(Temp[,p], Node_counts, counts)
        fitness_c[2*i] = Fitness(Children[[3]], Children[[4]], theta, leaky, Unique_Data, counts, model)
      }
    }
    
    Mats = Mats_c
    Perms = Perms_c
    fitness = fitness_c
    
    Results[generation, 1] = max(fitness)
    Results[generation, 2] = mean(fitness)
    
    if(max(fitness) > Max_fitness)
    {
      key = which(fitness == max(fitness))
      Max_fitness = max(fitness)
      Max_perm = Perms[[key[1]]]
      Max_edge = Mats[[key[1]]]
      max_gen = generation
    }
    generation = generation + 1
  }
  
  p = invPerm(Max_perm)
  Temp = Max_edge[p,]
  theta <- Infer_Theta(Temp[,p], Node_counts, counts)
  
  G <- igraph::graph_from_adjacency_matrix(t(Temp[,p]), mode = "directed")
  
  Return_List <- list()
  print(G)
  Return_List[[1]] <- as.vector(t(igraph::get.edgelist(G)))
  Return_List[[2]] <- Max_fitness
  return(Return_List)
}


##### Genetic Operators

### Function Matrix Crossover
### Takes two topologically sorted DAGs and produces 2 baby DAGs

Matrix_Crossover <- function(M1, M2)
{
  n <- nrow(M1)
  ###Define children
  C1 = M1
  C2 = M2
  
  for(i in 3:n)
  {
    if(i %% 2 == 0)
    {
      C1[i,] = 0
      x <- which(M1[i,] != 0)
      
      C1[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C1) == FALSE)
        {
          C1[i, j] = 1
        }
      }
      
      C2[i,] = 0
      x <- which(M2[i,] != 0)
      
      C2[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C2) == FALSE)
        {
          C2[i, j] = 1
        }
      }
    }
    else
    {
      C1[i,] = 0
      x <- which(M2[i,] != 0)
      
      C1[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C1) == FALSE)
        {
          C1[i, j] = 1
        }
      }
      
      C2[i,] = 0
      x <- which(M1[i,] != 0)
      
      C2[i,x[1]] = 1
      
      x <- x[-1]
      for(j in x)
      {
        if(Redundant(j, i, C2) == FALSE)
        {
          C2[i, j] = 1
        }
      }
    }
  }
  
  Children = list()
  Children[[1]] = C1
  Children[[2]] = C2
  return(Children)
}

Edge_Mutation <- function(M, hard_penalty)
{
  Candidate_Changes <- list()
  count <- 1
  n <- nrow(M)
  
  ### Consider all possible edge additions
  for(i in 3:n)
  {
    if(sum(M[i,]) < 3 & hard_penalty)
    {
      for(j in 2:(i-1))
      {
        if(M[i, j] == 0 & Redundant(j, i, M) == FALSE & Overpower(j, i, M) == FALSE)
        {
          M_new = M
          M_new[i, j] = 1
          Candidate_Changes[[count]] = M_new
          count = count + 1
        }
      }
    }
    else
    {
      for(j in 2:(i-1))
      {
        if(M[i, j] == 0 & Redundant(j, i, M) == FALSE & Overpower(j, i, M) == FALSE)
        {
          M_new = M
          M_new[i, j] = 1
          Candidate_Changes[[count]] = M_new
          count = count + 1
        }
      }     
    }
  }
  
  if(rbinom(1, 1, 0.5) == 1 & length(Candidate_Changes) != 0)
  {
    key = sample(1:length(Candidate_Changes), 1, replace = FALSE)
    return(Candidate_Changes[[key]])
  }
  
  
  ### Consider all possible edge removals
  for(i in 2:n)
  {
    parents = which(M[i,] != 0)
    
    if(length(parents) > 1)
    {
      for(j in parents)
      {
        M_new = M
        M_new[i, j] = 0
        Candidate_Changes[[count]] <- M_new
        count = count + 1
      }
    }
  }
  
  if(length(Candidate_Changes) == 0)
  {
    return(M)
  }
  
  key = sample(1:length(Candidate_Changes), 1, replace = FALSE)
  return(Candidate_Changes[[key]])
}

Branch_Adjust <- function(M, node)
{
  ### Remove ties with other branches
  children <- c(node, which(M[,node] != 0))
  
  repeat
  {
    children_orig = children
    
    for(i in children)
    {
      children <- c(children, which(M[,i] != 0))
    }
    
    children = unique(children)
    
    if(length(children_orig) == length(children))
    {
      break
    }
  }
  
  set <- c(1:(node - 1))
  parents = which(M[node,] != 0)
  set = set[-parents]
  if(length(set) == 0)
  {
    return(M)
  }
  
  for(i in children)
  {
    for(j in 1:nrow(M))
    {
      if(M[i,j] == 1 & (j %in% children == FALSE))
      {
        M[i,j] = 0
      }
    }
  }
  
  ### Move the branch somewhere
  new_node = sample(set, 1, replace = FALSE)
  M[node, new_node] = 1
  
  return(M)
}


Redundant <- function(parent_node, child_node, M)
{
  G = graph_from_adjacency_matrix(t(M), mode = "directed")
  Paths = all_simple_paths(G, 1, parent_node)
  
  if(length(Paths) == 0)
  {
    return(TRUE)
  }
  for(i in 1:length(Paths))
  {
    Paths[[i]] = as_ids(Paths[[i]])
  }
  
  parents <- which(M[child_node,] != 0)
  
  count <- 0
  for(i in 1:length(Paths))
  {
    if(sum(parents %in% Paths[[i]]) > 0)
    {
      count = count + 1
    }
  }
  
  if(count == length(Paths))
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

Overpower <- function(parent_node, child_node, M)
{
  vec = child_node
  
  repeat
  {
    vec_copy = vec
    for(i in vec)
    {
      vec = c(vec, which(M[i,] != 0))
    }
    vec = vec[-c(1:length(vec_copy))]
    vec = unique(vec)
    
    if(length(vec) == 1 & vec[1] == parent_node)
    {
      return(TRUE)
    }
    else if(length(vec) == 1 & vec[1] == 1)
    {
      break
    }
  }
  
  return(FALSE)
}



CX <- function(p1, p2, n)
{
  Children = matrix(1, nrow = 2, ncol = n)
  p1 = p1[-1]
  p2 = p2[-1]
  p1 = p1 - 1
  p2 = p2 - 1
  
  ###Identify cycles in first parent
  cycle = p1[1]
  cycle_over = FALSE
  
  while(cycle_over == FALSE)
  {
    key = which(p1 == cycle[length(cycle)])
    if(p2[key] != p1[1])
    {
      cycle = c(cycle, p2[key])
    }
    else
    {
      cycle_over = TRUE
    }
  }
  for(i in 1:(n-1))
  {
    Children[1, i+1] = ifelse(p1[i] %in% cycle, p1[i], p2[i])
  }
  
  Children[1, 2:n] = Children[1,2:n] + 1
  
  ###Rinse and repeat
  
  cycle = p2[1]
  cycle_over = FALSE
  
  while(cycle_over == FALSE)
  {
    key = which(p2 == cycle[length(cycle)])
    if(p1[key] != p2[1])
    {
      cycle = c(cycle, p1[key])
    }
    else
    {
      cycle_over = TRUE
    }
  }
  
  for(i in 1:(n-1))
  {
    Children[2, i+1] = ifelse(p2[i] %in% cycle, p2[i], p1[i])
  }
  
  Children[2, 2:n] = Children[2, 2:n] + 1
  
  return(Children)
}



Permutation_Mutation <- function(perm_chromosome, n)
{
  keys = sample(2:n, 2, replace = FALSE)
  
  perm_chromosome <- replace(perm_chromosome, keys, perm_chromosome[c(keys[2], keys[1])])
  
  return(perm_chromosome)
}


Perm_Repair <- function(G, p)
{
  n <- nrow(G)
  set <- c(2:n)
  
  while(length(set) != 0)
  {
    row = set[1]
    eq_rows = which(apply(G, 1, function(x) return(all(x == G[row,]))))
    
    if(length(eq_rows) != 1)
    {
      eq_cols = which(apply(G, 2, function(x) return(all(x == G[,row]))))
      eq_cols = intersect(eq_cols, eq_rows)
      eq_rows = eq_cols
      
      if(length(eq_rows) > 0)
      {
        ###Sort
        extract = p[eq_rows]
        extract = extract[order(extract, decreasing = FALSE)]
        p[eq_rows] = extract
      }
    }
    
    set = setdiff(set, eq_rows)
  }
  
  return(p)
}



DAGrR <- function(m1, m2, p1, p2, f1, f2, p_mc, p_cx, p_ba, p_ea, p_perm_mut, penalty)
{
  hard <- FALSE
  if(penalty == "hard")
  {
    hard = TRUE
  }
  
  ### Define return
  Children <- list()
  n <- nrow(m1)
  
  recomb = TRUE
  
  ### Check for equivalence
  ### THIS IS NOT CURRENTLY ACCURATE
  ### TO DO: IMPROVE THIS CHECK
  if(f1 == f2)
  {
    recomb = FALSE
  }
  
  child_m1 = m1
  child_m2 = m2
  child_p1 = p1
  child_p2 = p2
  
  if(rbinom(1, 1, p_mc) == 1 & recomb == TRUE)
  {
    Children = Matrix_Crossover(child_m1, child_m2)
    child_m1 = Children[[1]]
    child_m2 = Children[[2]]
  }
  
  if(rbinom(1, 1, p_cx) == 1 & recomb == TRUE)
  {
    child_p = CX(p1, p2, nrow(m1))
    child_p1 = child_p[1,]
    child_p2 = child_p[2,]
  }
  
  if(rbinom(1, 1, p_ba) == 1)
  {
    child_m1 = Branch_Adjust(child_m1, sample(3:n, 1, replace = FALSE))
  }
  
  if(rbinom(1, 1, p_ba) == 1)
  {
    child_m2 = Branch_Adjust(child_m1, sample(3:n, 1, replace = FALSE))
  }
  
  if(rbinom(1, 1, p_ea) == 1)
  {
    child_m1 = Edge_Mutation(child_m1, hard)
  }
  
  if(rbinom(1, 1, p_ea) == 1)
  {
    child_m2 = Edge_Mutation(child_m2, hard)
  }
  
  if(rbinom(1, 1, p_perm_mut) == 1)
  {
    child_p1 = Permutation_Mutation(child_p1, nrow(m1))
  }
  
  if(rbinom(1, 1, p_perm_mut) == 1)
  {
    child_p2 = Permutation_Mutation(child_p2, nrow(m1))
  }
  
  child_p1 = Perm_Repair(child_m1, child_p1)
  child_p2 = Perm_Repair(child_m2, child_p2)
  
  Children[[1]] = child_m1
  Children[[2]] = child_p1
  Children[[3]] = child_m2
  Children[[4]] = child_p2
  return(Children)
}




Fitness <- function(m, p, theta, leaky, Unique_Data, counts, model)
{
  ###Permute the matrix
  p = invPerm(p)
  Temp = m[p,]
  m = Temp[,p]
  
  val = 0
  for(i in 1:nrow(Unique_Data))
  {
    val = val + counts[i]*log(GA_Likelihood(Unique_Data[i,], m, theta, leaky, model))
  }
  
  return(val)
}

GA_Likelihood <- function(x, G, theta, leaky, model)
{
  N <- nrow(G)
  val <- 1
  for(i in 2:N)
  {
    G_t = G[i,]
    if(((min(1, sum(G_t * x)) == 1) && model == "DBN") || ((min(1, sum(G_t * x)) == sum(G_t)) && model == "CBN"))
    {
      if(x[i] == 1)
      {
        val = val*theta[i]
      }
      else
      {
        val = val*(1 - theta[i])
      }
    }
    else
    {
      if(x[i] == 1)
      {
        val = val*leaky
      }
      else
      {
        val = val*(1 - leaky)
      }
    }
  }
  return(val)
}

Create_Dart_Board <- function(fitness, N)
{
  fitness = fitness - min(fitness)
  dart_board = c(1:N)
  sum_vec = sum(fitness)
  
  if(sum_vec == 0)
  {
    dart_board = seq(0, 1, 1/N)
    dart_board = dart_board[-1]
    return(dart_board)
  }
  
  for(i in 1:length(fitness))
  {
    dart_board[i] = (sum(fitness[1:i])) / sum_vec
  }
  
  return(dart_board)
}

Selection <- function(dart_board)
{
  indices <- c(1,2)
  
  dart_1 <- runif(1, 0, 1)
  
  for(i in 2:length(dart_board))
  {
    if(dart_1 < dart_board[i])
    {
      indices[1] = i
      break
    }
  }
  
  dart_2 <- runif(1, 0, 1)
  
  for(i in 2:length(dart_board))
  {
    if(dart_2 < dart_board[i])
    {
      indices[2] = i
      break
    }
  }
  
  return(indices)
}


Generate_Tree_DAG <- function(n)
{
  M <- matrix(0, nrow = n, ncol = n)
  M[2, 1] = 1
  
  for(i in 3:n)
  {
    key = sample(1:(i - 1), 1, replace = FALSE)
    M[i, key] = 1
  }
  
  return(M)
}


Infer_Theta <- function(G, Node_counts, counts)
{
  n <- nrow(G)
  theta <- c(1:n)
  
  for(i in 2:n)
  {
    parents = Parents(G, i)
    if(1 %in% parents)
    {
      denom = sum(counts)
      num = length(Node_counts[[i]])
    }
    else
    {
      set = Node_counts[[parents[1]]]
      for(j in parents)
      {
        set = union(set, Node_counts[[j]])
      }
      denom = length(set)
      set = intersect(set, Node_counts[[i]])
      num = length(set)
    }
    
    theta[i] = num/denom
    
    if(num == 0)
    {
      theta[i] = 1/sum(counts)
    }
    if(denom == 0)
    {
      theta[i] = 1/sum(counts)
    }
    if(theta[i] == 1)
    {
      theta[i] = (sum(counts) - 1)/sum(counts)
    }
  }
  return(theta)
}



Infer_Spontaneous <- function(G, Node_counts, counts)
{
  n <- nrow(G)
  spont <- c(1:n)
  spont[1] = 0
  
  for(i in 2:n)
  {
    parents = Parents(G, i)
    if(1 %in% parents)
    {
      spont[i] = 0
    }
    else
    {
      set = setdiff(c(1:sum(counts)), Node_counts[[parents[1]]])
      for(j in parents)
      {
        set = intersect(set, setdiff(c(1:sum(counts)), Node_counts[[j]]))
      }
      set = intersect(set, Node_counts[[i]])
      num = length(set)
      denom = length(Node_counts[[i]])
      
      if(denom == 0 | num == 0)
      {
        spont[i] = 0.01
      }
      else
      {
        spont[i] = num/denom
      }
      
      if(spont[i] == 1)
      {
        spont[i] = 0.99
      }
    }
  }
  return(spont)
}

Parents <- function(graph, node)
{
  return(which(graph[node,] != 0))
}


#### CODE SPECIFICALLY FOR ME MODEL
EM <- function(G, theta_init, CPT, Unique_Data, counts, CONVERGENCE)
{
  
  n <- length(theta_init)
  difference <- replicate(n, 1)
  counter <- 1
  vec = convert_mat(as.integer(n), as.integer(G))
  theta <- theta_init
  
  while(max(difference) > CONVERGENCE)
  {
    likelihood_vec = likelihood(theta, vec)
    theta_prior = theta
    P_y = as.vector(t(CPT %*% likelihood_vec))
    
    for(j in 2:n)
    {
      parents = Parents(G, j)
      if(1 %in% parents)
      {
        denom = sum(counts)
        set_num = indices_denom(as.integer(n), as.integer(j))
        a <- as.vector(t(CPT[,set_num] %*% likelihood_vec[set_num]))
        num = sum(counts*(P_y)^{-1}*a)
      }
      else
      {
        set_denom = indices_denom(as.integer(n), as.integer(parents))
        set_num = indices_num(as.integer(n), as.integer(parents), as.integer(j))
        set_denom = unique(set_denom)
        set_num = unique(set_num)
        a <- as.vector(t(CPT[,set_num] %*% likelihood_vec[set_num]))
        b <- as.vector(t(CPT[,set_denom] %*% likelihood_vec[set_denom]))
        num = sum(counts*(P_y)^{-1}*a)
        denom = sum(counts*(P_y)^{-1}*b)
      }
      
      theta[j] = num/denom
    }
    
    difference = abs(theta - theta_prior)
    counter = counter+1
  }
  
  likelihood_vec = likelihood(theta, vec)
  P_y = as.vector(t(CPT %*% likelihood_vec))
  Return_list <- list()
  Return_list[[1]] = sum(counts*log(P_y))
  Return_list[[2]] = theta
  return(Return_list)
}


### Code copied from
###https://stackoverflow.com/questions/12088080/how-to-convert-integer-number-into-binary-vector

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}


Dist <- function(A, B)
{
  if(length(A) == 0)
  {
    return(0)
  }
  
  return(length(which(A != B)))
}

Y_Given_X <- function(hidden_data, true_data, e, d)
{
  hidden_data = hidden_data[-1]
  true_data = true_data[-1]
  
  hidden_data_1 <- hidden_data[which(hidden_data == 1)]
  hidden_data_0 <- hidden_data[which(hidden_data == 0)]
  
  distance_1 <- Dist(hidden_data_1, true_data[which(hidden_data == 1)])
  distance_0 <- Dist(hidden_data_0, true_data[which(hidden_data == 0)])
  
  prob = e^{distance_0}*(1-d)^{length(hidden_data_0) - distance_0}
  prob = prob*d^{distance_1}*(1-e)^{length(hidden_data_1) - distance_1}
  
  return(prob)        
}


