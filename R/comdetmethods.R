#'
#' Function to perform community detection on multilayer network data
#' @references
#'
#' @author Jes\'us Arroyo <jarroyo@tamu.edu>
#' @export
comdetmethods <- function(Adj_list, K, method = "mase", parameters) {
  n <- ncol(Adj_list[[1]])
  if(method == "ave_spherical") {
    Abar <- Reduce("+", Adj_list)
    V <- eig_embedding(Abar, K)
    rownorms <- apply(V, 1,function(x) sqrt(sum(x^2)))
    community_memberships <- kmeans(V/rownorms, K, nstart = 100)$cluster
    #require(Gmedian)
    #adjustedRandIndex(community_memberships, truecom)
    #community_memberships <- kGmedian(V/rownorms, K, nstart = 100)$cluster
    #community_memberships <- Mclust(V/rownorms, K)$classification
  }
  if(method == "ave") {
    Abar <- Reduce("+", Adj_list)
    V <- eig_embedding(Abar, K)
    community_memberships <- kmeans(V, K, nstart = 100)$cluster
  }
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
    community_memberships <- kmeans(V, K, nstart = 100)$cluster
    rownorms <- apply(V, 1,function(x) sqrt(sum(x^2)))
    community_memberships <- kmeans(V/rownorms, K, nstart = 100)$cluster
  }
  
  if(method == "score") {
    Abar <- Reduce("+", Adj_list)
    V <- eig_embedding(Abar, K)
    V <- V[, 2:K, drop = F] / V[,1]
    community_memberships <- kmeans(V, K, nstart = 100)$cluster
  }
  if(method == "mase") {
    require(SpectralGraphInference)
    mase.res <- mase(Adj_list, d = K, scaled.ASE = FALSE, diag.augment = FALSE)
    community_memberships <- kmeans(mase.res$V, K, nstart = 100)$cluster
  }
  if(method == "mase-gmm") {
    require(SpectralGraphInference)
    require(mclust)
    mase.res <- mase(Adj_list, d = K, scaled.ASE = FALSE, diag.augment = FALSE)
    community_memberships <- Mclust(mase.res$V, K, verbose = FALSE)$classification
  }
  if(method == "dcmase") {
    community_memberships <- comdet_dcmase(Adj_list, K, "kmeans")$community_memberships
  }
  if(method == "dcmase-gmm") {
    community_memberships <- comdet_dcmase(Adj_list, K, "gmm")$community_memberships
  }
  if(method == "laplacian-mase") {
    community_memberships <- speck(Adj_list, n, K)
  }
  if(method == "lmfo") {
    community_memberships <- lmfo(lapply(Adj_list, as.matrix), n, K)
  }
  return(community_memberships)
}



