#'
#' Function to perform community detection on multilayer network data
#' @references
#'
#' @author Jes\'us Arroyo <jarroyo@tamu.edu>
#' @export
comdetmethods <- function(Adj_list, K, method = "dcmase") {
  n <- ncol(Adj_list[[1]])
  
  #################################################################################
  if(method == "dcmase") {
    community_memberships <- comdet_dcmase(Adj_list, K, "kmeans")$community_memberships
  }
  #################################################################################
  if(method == "dcmase-gmm") {
    community_memberships <- comdet_dcmase(Adj_list, K, "gmm")$community_memberships
  }
  #################################################################################
  if(method == "dcmase-unscaled") {
    community_memberships <- comdet_dcmase(Adj_list, K,scaled = FALSE, "kmeans")$community_memberships
  }
  #################################################################################
  if(method == "mase") {
    require(SpectralGraphInference)
    mase.res <- mase(Adj_list, d = K, scaled.ASE = FALSE, diag.augment = FALSE)
    community_memberships <- kmeans(mase.res$V, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "mase-gmm") {
    require(SpectralGraphInference)
    require(mclust)
    mase.res <- mase(Adj_list, d = K, scaled.ASE = FALSE, diag.augment = FALSE)
    community_memberships <- Mclust(mase.res$V, K, verbose = FALSE)$classification
  }
  #################################################################################
  if(method == "mase-spherical") {
    require(SpectralGraphInference)
    require(mclust)
    mase.res <- mase(Adj_list, d = K, scaled.ASE = FALSE, diag.augment = FALSE)
    V <- mase.res$V
    rownorms <- apply(V, 1,function(x) sqrt(sum(x^2)))
    rownorms <- rownorms + 1*(abs(rownorms) < 1e-15)
    community_memberships <- kmeans(V/rownorms, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "mase-score") {
    community_memberships <- comdet_dcmase(Adj_list, K,row.normalization = "score", "kmeans")$community_memberships
  }
  #################################################################################
  if(method == "ave") {
    Abar <- Reduce("+", Adj_list)
    V <- eig_embedding(Abar, K)
    community_memberships <- kmeans(V, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "ave_spherical") {
    Abar <- Reduce("+", Adj_list)
    V <- eig_embedding(Abar, K)
    rownorms <- apply(V, 1,function(x) sqrt(sum(x^2)))
    rownorms <- rownorms + 1*(abs(rownorms) < 1e-15)
    
    community_memberships <- kmeans(V / rownorms, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "score") {
    Abar <- Reduce("+", Adj_list)
    V <- eig_embedding(Abar, K)
    V <- V[, 2:K, drop = F] / (V[,1] + 1*(abs(V[,1]) < 1e-15))
    community_memberships <- kmeans(V, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "sq-bias-adjusted") {
    # https://arxiv.org/pdf/2003.08222.pdf
    Asq_bias_removed <- lapply(Adj_list, function(A) {
      d = colSums(A) # degrees
      A2 <- crossprod(A)
      diag(A2) <- diag(A2) - d
      A2
    })
    Asq_bias_sum <- Reduce("+", Asq_bias_removed)
    V <- eig_embedding(Asq_bias_sum, K)
    rownorms <- apply(V, 1,function(x) sqrt(sum(x^2)))
    rownorms <- rownorms + 1*(abs(rownorms) < 1e-15)
    community_memberships <- kmeans(V/rownorms, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "omni") {
    if(length(Adj_list)==1) {
      omnibar = g.ase(Adj_list[[1]], K, diag.augment = F)$X
    } else{
      omni <- OMNI_matrix(Adj_list, K)
      omnibar <- Reduce("+", lapply(omni, function(x) x$X)) 
    }
    community_memberships <- kmeans(omnibar, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "lmfo") {
    # https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-1/Spectral-and-matrix-factorization-methods-for-consistent-community-detection-in/10.1214/18-AOS1800.short
    community_memberships <- lmfo(lapply(Adj_list, as.matrix), n, K)
  }
  #################################################################################
  if(method == "speck") {
    # # https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-1/Spectral-and-matrix-factorization-methods-for-consistent-community-detection-in/10.1214/18-AOS1800.short
    community_memberships <- speck(Adj_list, n, K)
  }
  #################################################################################
  if(method == "graph-tool") {
    # Note: this method requires graph-tool library to be installed https://graph-tool.skewed.de/
    community_memberships <- run_graph_tool(Adj_list, K)
  }
  #################################################################################
  if(method == "score-onenetwork") {
    m <-length(Adj_list)
    V <- eig_embedding(Adj_list[[sample(m,1)]], K)
    V <- V[, 2:K, drop = F] / (V[,1] + 1*(abs(V[,1]) < 1e-15))
    community_memberships <- kmeans(V, K, nstart = 100)$cluster
  }
  #################################################################################
  if(method == "spherical-onenetwork") {
    m <-length(Adj_list)
    V <- eig_embedding(Adj_list[[sample(m,1)]], K)
    rownorms <- apply(V, 1,function(x) sqrt(sum(x^2)))
    rownorms <- rownorms + 1*(abs(rownorms) < 1e-15)
    community_memberships <- kmeans(V/rownorms, K, nstart = 100)$cluster
  }
  
  return(community_memberships)
}



