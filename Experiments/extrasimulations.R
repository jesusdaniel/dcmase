# Extra simulations
#########################################################################
### same B same theta, changing B
# K = parameters[3]
# p = parameters[1]
# q = parameters[2]
# B = (p-q)*diag*K + q 
extrasimulation1 <- function(parameters, seed = 1989) {
  set.seed(seed)
  
  K = parameters[3]
  p = parameters[1]
  q = parameters[2]
  m <- parameters[4]
  n <- parameters[5]
  
  B = (p-q)*diag(K) + q
  
  
  ave_deg <- 2*log(n)
  # Generate network parameters ------------------------------------------------
  #n <- 450
  Z <- kronecker(diag(K), rep(1, n/K))
  
  theta <- rep(1, n)#
  theta <- runif(n, min = 0.1, max = 1)
  theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  
  
  degree_corrections <- lapply(1:m, function(i) theta)
  B_matrices <- lapply(1:m, function(i) B)
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}

# 
# m = 10
# plot_adjmatrix(Adj_list[[1]])
#colSums(Adj_list[[1]])




extrasimulation2 <- function(parameters, seed = 1989) {
  set.seed(seed)
  
  K = parameters[3]
  p = parameters[1]
  q = parameters[2]
  m <- parameters[4]
  n <- parameters[5]
  
  B = (p-q)*diag(K) + q
  
  
  # Generate network parameters ------------------------------------------------
  #n <- 450
  Z <- kronecker(diag(K), rep(1, n/K))
  
  
  P = Z %*% (B/p) %*% t(Z)
  
  theta = (rexp(n)) + 0.4
  Theta = diag(theta/max(theta))
  
  P = Theta %*% P %*% Theta

  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}

