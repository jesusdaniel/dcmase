#######################################
# Supplementary material
# Comparison of scaled vs unscaled
# eigenvector normalization
#######################################

library(SpectralGraphInference)

sample_graph <- function(n, K, offdiag) {
  truecoms <- sample(1:K, n,replace = T, prob = c(5/11, 5/11, 1/22, 1/22))
  Z <- matrix(0, n, K)
  Z[cbind(1:n, truecoms)] = 1
  B <- matrix(offdiag, K, K) + (1-offdiag) * diag(K)
  theta <- runif(n, min = 0.8, max = 1)
  P <- tcrossprod((theta * Z) %*% B, theta * Z)
  A <- sample_from_P(P)
  list(A = A, Z = Z)
}

spherical_spectral_clustering <- function(A, K, truecoms) {
  # Embeddings------------------------------------------------------------------
  X <- g.ase(A, K)$X
  U = eig_embedding(A, K)
  
  # Row-normalization ----------------------------------------------------------
  normsU <- apply(U, 1, function(x) sqrt(sum(x^2)))
  Unorm <- U/(normsU + 1*(normsU==0))
  normsX <- apply(X, 1, function(x) sqrt(sum(x^2)))
  Xnorm <- X/(normsX + 1*(normsX==0))
  
  # Clustering -----------------------------------------------------------------
  Unorm.estcoms = kmeans(Unorm, K, nstart = 100)$cluster
  Xnorm.estcoms = kmeans(Xnorm, K, nstart = 100)$cluster
  error_unscaled = mclust::classError(truecoms, Unorm.estcoms)$errorRate
  error_scaled = mclust::classError(truecoms, Xnorm.estcoms)$errorRate  
  return(list(Unorm = Unorm, Xnorm = Xnorm,
              Unorm.estcoms = Unorm.estcoms, Xnorm.estcoms = Xnorm.estcoms,
              error_unscaled = error_unscaled, error_scaled = error_scaled))
}

set.seed(123)
offdiag <- 0.3
n <- 800
K <- 4
data = sample_graph(n, K, offdiag)
A = data$A
Z = data$Z
truecoms <- Z %*% 1:K
clusres = spherical_spectral_clustering(A, K, truecoms)

#png("Experiments SCORE vs spherical/unscaled-sim.png")
pairs(clusres$Unorm, col = truecoms, pch = truecoms)
#dev.off()
#png("Experiments SCORE vs spherical/scaled-sim.png")
pairs(clusres$Xnorm, col = truecoms, pch = truecoms)
#dev.off()


off_diag_vals = seq(0.1, 0.8, 0.1)
results = matrix(0, 2, length(off_diag_vals))
u = 1
for(offdiag in off_diag_vals) {
  print(u)
  results[, u] = rowMeans(sapply(1:10, function(i) {
    data = sample_graph(n, K, offdiag)
    A = data$A
    Z = data$Z
    truecoms <- Z %*% 1:K
    clusres = spherical_spectral_clustering(A, K, truecoms)
    return(c(clusres$error_unscaled, clusres$error_scaled))
  }))
  u = u+1
}
library(xtable)
xtable(rbind(off_diag_vals, results))

