library(Matrix)
library(igraph)
library(mclust)
source("R/Codes_Spectral_Matrix_Paul_Chen_AOS_2020.r")
source("R/comdet-dcmase.R")
source("R/comdetmethods.R")
source("R/dcmase.R")
source("R/SpectralMethods.R")
source("R/run_graph_tool.R")
source("R/make_ggplot.R")


run_simulations <- function(sim_setting, parameters, repetitions = 20) {
  library(parallel)
  cl = makeCluster(10)
  clusterEvalQ(cl = cl, source("Experiments/run_all_methods.R"))
  clusterEvalQ(cl = cl, source("Experiments/extrasimulations.R"))
  clusterExport(cl = cl, varlist = c("sim_setting", "parameters"),envir = environment()) 
  results <- parLapply(cl, 1:repetitions, function(seed) {
    generate_data <- sim_setting(parameters, seed)
    run_all_methods(generate_data$Adj_list, generate_data$truecom)
  })
  
  df_res <- data.frame(Reduce(rbind, results))
  rownames(df_res) <- 1:repetitions
  stopCluster(cl)
  return(df_res)
}

iterate_parameters <- function(sim_setting, parameters_list, param_iter, 
                               repetitions = 20) {
  df_res <- lapply(1:length(parameters_list),  function(i) {
    cat("Running parameter ", param_iter[[i]], "...\n", sep = "")
    sim_res <- run_simulations(sim_setting, parameters = parameters_list[[i]], repetitions)
    sim_res$parameter <- param_iter[[i]]
    return(sim_res)
  })
  return(Reduce(rbind, df_res))
}



run_all_methods <- function(Adj_list, truecoms) {
  K <- length(unique(truecoms))
  
  # Note: to run graph-tool, uncomment the following lines and comment the next
  # methods_to_run <- c("dcmase", "ave_spherical", "sq-bias-adjusted",
  #                    "mase-spherical","lmfo", "graph-tool")
  methods_to_run <- c("dcmase", "ave_spherical", "sq-bias-adjusted",
                      "mase-spherical","lmfo")
  results <- sapply(methods_to_run, function(method) {
    print(method)
  
    classError(comdetmethods(Adj_list, K, method = method), truecoms)$errorRate
    }
    )
  # Note: to run graph-tool, uncomment the following lines and comment the next
  #names(results) <- c("DC-MASE", "Sum A", "S-A^2-Bias-adj",
  #                    "MASE", "OLMF", "graph-tool")
  names(results) <-  c("DC-MASE", "Sum A", "S-A^2-Bias-adj",
                          "MASE", "OLMF")
  return(results)
}



#########################################################################
### Simulation 1: same B same theta
simulation1 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  #theta <- runif(n, min = 0.05, max = 1)
  theta <- rexp(n) + 0.2
  theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  
  B <- 0.06*diag(K) + 0.04
  
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

### Simulation 2: same theta, different B
simulation2 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  
  #theta <- runif(n, min = 0.05, max = 1)
  theta <- rexp(n) + 0.2
  theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  degree_corrections <- lapply(1:m, function(i) theta)
  
  B_matrices <- lapply(1:m, function(i) {
    p <- runif(1)
    q <- runif(1)
    (p-q) *diag(K) + q
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}





### Simulation 3: different B, different theta
simulation3 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  #degree_corrections <- lapply(1:m, function(i) runif(n, min = sqrt(0.05), max = 1)^2)
  degree_corrections <- lapply(1:m, function(i) {
    theta <- runif(n, min = 0.05, max = 1)
    theta <- rexp(n) + 0.2
    as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  })
  B_matrices <- lapply(1:m, function(i) {
    p <- runif(1)
    q <- runif(1)
    (p-q) *diag(K) + q
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}


### Simulation 4: same B, different theta
simulation4 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  ave_deg <- 10
  # Generate network parameters ------------------------------------------------
  n <- 150
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  
  degree_corrections <- lapply(1:m, function(i) {
    theta <- runif(n, min = 0.05, max = 1)
    theta <- rexp(n) + 0.2
    as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  })
  B_matrices <- lapply(1:m, function(i) {
    B <- 0.06*diag(K) + 0.04
    B
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  return(list(Adj_list = Adj_list, truecom = as.vector(Z %*% 1:K)))
}






# Simulation 6
### Simulation 6: changing high degree
simulation6 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  # Set network parameters ------------------------------------------------
  n <- parameters[2] #150
  K <- parameters[3] #3
  
  ave_deg <- 10
  Z <- kronecker(diag(K), rep(1, n/K))
  
  theta1 <- rep(c(rep(0.15, 0.5*n/K), rep(0.8, 0.5*n/K)), K)
  theta2 <- 0.95 - theta1
  theta1 <- as.vector(theta1 / (Z%*%crossprod(Z,theta1)/(n/K)))
  theta2 <- as.vector(theta2 / (Z%*%crossprod(Z,theta2)/(n/K)))               
  degree_corrections <- lapply(1:m, function(i) if(i%%2==1){
    theta1
  }else{
    theta2
  })
  
  B_matrices <- lapply(1:m, function(i) {
    B <- 0.06*diag(K) + 0.04
    B
  })
  
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  #plot_adjmatrix(Adj_list[[1]])
  #colSums(Adj_list[[1]])
  truecom <- as.vector(Z %*% 1:K)
  
  results <- list(Adj_list = Adj_list, truecom = truecom)
  return(results)
}





# Simulation 7
### Simulation 7: changing high degree and random B
simulation7 <- function(parameters, seed = 1989) {
  set.seed(seed)
  m <- parameters[1]
  # Set network parameters ------------------------------------------------
  n <- parameters[2] #150
  K <- parameters[3] #3
  
  ave_deg <- 10
  Z <- kronecker(diag(K), rep(1, n/K))
  
  theta1 <- rep(c(rep(0.15, 0.5*n/K), rep(0.8, 0.5*n/K)), K)
  theta2 <- 0.95 - theta1
  theta1 <- as.vector(theta1 / (Z%*%crossprod(Z,theta1)/(n/K)))
  theta2 <- as.vector(theta2 / (Z%*%crossprod(Z,theta2)/(n/K)))               
  degree_corrections <- lapply(1:m, function(i) if(i%%2==1){
    theta1
  }else{
    theta2
  })
  
  B_matrices <- lapply(1:m, function(i) {
    p <- runif(1)
    q <- runif(1)
    ((p-q) *diag(K) + q)
  })
  
  # Sample graphs --------------------------------------------------------------
  Adj_list <- lapply(1:m, function(i) {
    P <- tcrossprod((degree_corrections[[i]] * Z) %*% B_matrices[[i]], degree_corrections[[i]] * Z)
    P <- (ave_deg*n/sum(P)) * P
    sample_from_P(P)
  })
  #plot_adjmatrix(Adj_list[[1]])
  #colSums(Adj_list[[1]])
  truecom <- as.vector(Z %*% 1:K)
  
  results <- list(Adj_list = Adj_list, truecom = truecom)
  return(results)
}


