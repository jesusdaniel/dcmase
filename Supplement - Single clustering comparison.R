#######################################
# Supplementary simulations for
# "Joint Spectral Clustering in
# Multilayer Degree Corrected Blockmodles"
#######################################
# This code compares different embedding
# methods in a single network

#######################################
# Simulate 3 communities
library(mclust)
library(ScorePlus)

iterate_parameters <- function(sim_setting, parameters_list, param_iter, 
                               repetitions = 20) {
  df_res <- lapply(1:length(parameters_list),  function(i) {
    cat("Running parameter ", param_iter[[i]], "...\n", sep = "")
    sim_res <- run_simulations_onegraph(sim_setting, parameters = parameters_list[[i]], repetitions)
    sim_res$parameter <- param_iter[[i]]
    return(sim_res)
  })
  return(Reduce(rbind, df_res))
}


run_simulations_onegraph <- function(sim_setting, parameters, repetitions = 20) {
  library(parallel)
  cl = makeCluster(10)
  clusterEvalQ(cl = cl, library(mclust))
  clusterEvalQ(cl = cl, library(SpectralGraphInference))
  clusterEvalQ(cl = cl, library(ScorePlus))
  clusterExport(cl = cl, varlist = c("sim_setting", "parameters", "run_all_methods",
                                     "spherical_sc", "spherical_sc_unweighted",
                                     "SCORE"),envir = environment()) 
  #results <- parLapply(cl, 1:repetitions, function(seed) {
  results <- lapply(1:repetitions, function(seed) {
    generate_data <- sim_setting(parameters, seed)
    run_all_methods(generate_data$A, generate_data$truecom)
  })
  
  df_res <- data.frame(Reduce(rbind, results))
  rownames(df_res) <- 1:repetitions
  stopCluster(cl)
  return(df_res)
}

spherical_sc <- function(A, K) {
  V <- g.ase(A, K)$X
  normsV <- apply(V, 1, function(x) sqrt(sum(x^2)))
  Y <- V/(normsV + 1*(normsV==0))
  kmeans(Y, K, nstart = 100)$cluster
}

spherical_sc_unweighted <- function(A, K) {
  V <- eig_embedding(A, K)
  normsV <- apply(V, 1, function(x) sqrt(sum(x^2)))
  Y <- V/(normsV + 1*(normsV==0))
  kmeans(Y, K, nstart = 100)$cluster
}

run_all_methods <- function(A, truecoms) {
  K <- length(unique(truecoms))
  results <- c(adjustedRandIndex(spherical_sc(A, K), truecoms),
               adjustedRandIndex(spherical_sc_unweighted(A, K), truecoms),
               adjustedRandIndex(SCORE(A, K)$labels, truecoms),
               adjustedRandIndex(SCOREplus(A, K)$labels, truecoms),
               adjustedRandIndex(SCORE(A, K,threshold = Inf)$labels, truecoms))
  
  names(results) <- c("Spherical", "Spherical unscaled",
                      "SCORE", "SCORE+", "SCORE (no thresh.)")
  return(results)
}


## eigenvector size, change B11
simulation1 <- function(params, seed) {
  set.seed(seed)
  ave_deg <- params[1]
  ratioB11 <- params[2] 
  offdiag <- params[3]
  degree_distr <- params[4]
  n <- 300
  K <- 3
  Z <- kronecker(diag(K), rep(1, n/K))
  B <- matrix(offdiag, K, K) + (1-offdiag) * diag(K)
  B[1,1] = ratioB11
  
  if(degree_distr == 1) theta <- rexp(n) + 0.05
  if(degree_distr == 2) theta <- runif(n, min = 0.1, 1)
  if(degree_distr == 3) theta <- 0.1 + 0.9*(runif(n) > 0.95)
  
  theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  
  P <- tcrossprod((theta * Z) %*% B, theta * Z)
  P <- (ave_deg*n/sum(P)) * P
  A <- sample_from_P(P)
  truecoms <- Z %*% 1:K
  return(list(A = A, truecom = truecoms))
}

# DIF COMMUNITY SIZES
simulation2 <- function(params, seed) {
  set.seed(seed)
  ave_deg <- params[1]
  offdiag <- params[2]
  ratioB11 <- 1
  n <- 300
  K <- 4
  
  commprobs <- params[3:(K+2)]
  truecoms <- sample(1:K, n,replace = T, prob = commprobs)
  Z <- matrix(0, n, K)
  Z[cbind(1:n, truecoms)] = 1
  
  B <- matrix(offdiag, K, K) + (1-offdiag) * diag(K)
  #B[2,1] = B[1,2] = 0.4
  
  theta <- runif(n, min = 0.1, max = 1)
  #theta <- as.vector(theta / (Z%*%crossprod(Z,theta)/(n/K)))
  
  P <- tcrossprod((theta * Z) %*% B, theta * Z)
  P <- (ave_deg*n/sum(P)) * P
  A <- sample_from_P(P)
  truecoms <- Z %*% 1:K
  return(list(A = A, truecom = truecoms))
}

# diff degree distribution
simulation3 <- function(params, seed) {
  set.seed(seed)
  ave_deg <- params[1]
  offdiag <- params[2]
  degpower <- params[3]
  ratioB11 <- 1
  n <- 300
  K <- 3
  
  Z <- kronecker(diag(K), rep(1, n/K))
  
  B <- matrix(offdiag, K, K) + (1-offdiag) * diag(K)
  B[1,1] = ratioB11
  
  theta <- runif(n, min = 0.1, max = 1)^degpower
  
  P <- tcrossprod((theta * Z) %*% B, theta * Z)
  P <- (ave_deg*n/sum(P)) * P
  A <- sample_from_P(P)
  truecoms <- Z %*% 1:K
  return(list(A = A, truecom = truecoms))
}

make_ggplot_res <- function(results, parameter_name, xbreaks = c(1, 5, 10, 15, 20, 25)) {
  require(ggplot2)
  require(reshape2)
  require(scales)
  results_melted <- melt(results, id.vars = "parameter")
  names(results_melted) <- c("parameter", "Method", "ARI")
  
  resdf = data_summary(results_melted, "ARI", c("parameter", "Method"))
  
  p<- ggplot(resdf, aes(x=parameter, y = ARI)) + 
    geom_line(aes(color=Method, linetype=Method)) +
    geom_point(aes(color=Method, shape=Method))+
    #geom_errorbar(aes(ymin=ARI-2*se, ymax=ARI+2*se, color = Method), width=0.2) +
    ylim(c(0,1)) +
    scale_x_continuous(breaks= xbreaks) +
    xlab(parameter_name) +theme_bw()
  p
}

make_ggplot_multiple <- function(different_scenarios, parameter_name, xbreaks = c(1, 5, 10, 15, 20, 25)) {
  require(ggplot2)
  require(reshape2)
  require(scales)
  results_melted <- melt(different_scenarios , id.vars = c("parameter", "type"))
  names(results_melted) <- c("parameter", "Type", "Method", "ARI")
  
  resdf = data_summary(results_melted, "ARI", c("parameter", "Method", "Type"))
  
  p<- ggplot(resdf, aes(x=parameter, y = ARI)) + 
    geom_line(aes(color=Method, linetype=Method)) +
    geom_point(aes(color=Method, shape=Method))+
    #geom_errorbar(aes(ymin=ARI-2*se, ymax=ARI+2*se, color = Method), width=0.2) +
    ylim(c(0,1)) +
    #scale_x_continuous(breaks= xbreaks) +
    xlab(parameter_name) +theme_bw() +
    facet_wrap("Type", scales = "free", strip.position = "bottom")
  p
}

#######################################
# Simulations
#######################################

num_replications = 500

#######################################
# Change B_11
param_iter = seq(1,7, 0.5)
parameters_list <- lapply(param_iter, function(x) c(15, x, 0.1, 2))
parameters_list <- parameters_list
results_simulation1 <- iterate_parameters(sim_setting = simulation1, parameters_list, param_iter, num_replications)

make_ggplot_res(results_simulation1, "B_11", xbreaks = param_iter)



#######################################
# Changecommunity sizes
param_iter = seq(1/4, 0.49, 0.03)
parameters_list <- lapply(param_iter, function(x) c(15, 0.1, c(x, x,1/2 - x, 1/2 - x)))

results_simulation2 <- iterate_parameters(sim_setting = simulation2, parameters_list, param_iter, num_replications)
make_ggplot_res(results_simulation2,  "Largest community size (K=4)", xbreaks = param_iter)


#######################################
# Change community degree distribution
param_iter = seq(0, 8, 0.5)
parameters_list <- lapply(param_iter, function(x) c(15, 0.1, x))
results_simulation3 <- iterate_parameters(sim_setting = simulation3, parameters_list, param_iter, num_replications)

make_ggplot_res(results_simulation3,  "Degree power", xbreaks = param_iter)



#######################################
# Change off-diag
param_iter = seq(0.05, 0.7, 0.05)
parameters_list <- lapply(param_iter, function(x) c(15, 1, x, 2))
parameters_list <- parameters_list
results_simulation4 <- iterate_parameters(sim_setting = simulation1, parameters_list, param_iter, num_replications)
make_ggplot_res(results_simulation4,  "Off-diagonal elements of B", xbreaks = param_iter)

#######################################
# Combine
#######################################
results_simulation1$type = "Community magnitudes"
results_simulation2$type = "Largest community size (K=4)"
results_simulation3$type = "Degree distribution power"
results_simulation4$type = "Between-community connectivity"

results_all <- rbind(results_simulation1, results_simulation2, results_simulation3, results_simulation4)
results_all$type <- factor(results_all$type, c("Degree distribution power",
                                               "Largest community size (K=4)",
                                               "Community magnitudes", 
                                               "Between-community connectivity"))

#png("Simulation-singlenetwork-comparisons.png", width = 1800, height = 1200, res = 200)
make_ggplot_multiple(results_all, "")
#dev.off()
