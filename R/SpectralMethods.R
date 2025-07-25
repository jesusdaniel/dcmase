#'
#' Function to perform multiple adjacency spectral embedding
#'
#' @param Adj_list a list of adjacency matrices with the same size n x n
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically
#' @param d_vec vector of individual embedding dimensions for each graph. If NA, it is the same as d.
#' @param scaled.ASE logical. When true, the scaled adjacency spectral embedding is used.
#' @param diag.augment logical. When true, the diagonal augmentation method for ASE is performed.
#' @param elbow_graph number of elbow selected in Zhu & Ghodsi method for the scree plot of each individual graph eigenvalues.
#' @param elbow_mase number of elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of MASE.
#' @param show.scree.results when TRUE, the histogram of the estimated d for each graph, and the scree plot of the singular values of  the graph is shown if d is not specified.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#'
#' @return A list containing a matrix V of size n x d, with the
#' estimated invariant subspace, and a list R with the individual score parameters for each graph (size d x d each).
#'
#' @references
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' @export
mase <- function(Adj_list, d = NA, d_vec = NA,
                 scaled.ASE = TRUE, diag.augment = TRUE,
                 elbow_graph = 1, elbow_mase = 2,
                 show.scree.results = FALSE,
                 par = FALSE, numpar = 12) {
  if(is.na(d_vec)) {
    d_vec = rep(d, length(Adj_list))
  }
  # running in parallel
  if(par) {
    cl <- parallel::makeCluster(numpar)
    parallel::clusterEvalQ(cl, source("R/loadAll.R"))
    parallel::clusterExport(cl = cl, varlist = list("ase", "eig_embedding", "getElbows", "Adj_list",
                                                    "elbow_graph", "d_vec", "diag.augment"), envir = environment())
    if(scaled.ASE) {
      latpos.list <- parallel::parLapply(cl = cl, 1:length(Adj_list), function(i)
        ase(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }else{
      latpos.list <- parallel::parLapply(cl = cl, 1:length(Adj_list), function(i)
        eig_embedding(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }
    parallel::stopCluster(cl)
  } else {
    if(scaled.ASE) {
      latpos.list <- lapply(1:length(Adj_list), function(i)
        ase(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }else{
      latpos.list <- lapply(1:length(Adj_list), function(i)
        eig_embedding(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph))
    }
  }
  V_all  <- Reduce(cbind, latpos.list)
  jointsvd <- svd(V_all)
  if(is.na(d)) {
    if(show.scree.results) {
      hist(sapply(latpos.list, ncol), main = "Estimated d for each graph")
    }
    d = getElbows(jointsvd$d, plot = show.scree.results)[elbow_mase]
  }
  V = jointsvd$u[, 1:d, drop = FALSE]
  R <- project_networks(Adj_list, V)
  return(list(V = V, sigma = jointsvd$d, R = R))
}



#'
#' Function to perform graph adjacency spectral embedding (ASE)
#'
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to try when d is not provided. Default is sqrt(ncol(A)).
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#' @return A matrix with n rows and d columns containing the estimated latent positions
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#'
#' @export
ase <- function(A, d = NA, d.max = sqrt(ncol(A)), diag.augment = TRUE, elbow = 1) {
  # Diagonal augmentation
  if(diag.augment & sum(abs(Matrix::diag(A))) == 0) {
    deg = Matrix::colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.na(d)) {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    V <- eig$vectors[,selected.eigs, drop = F]
    D <- diag(sqrt(abs(eig$values[selected.eigs])), nrow = d)
    X <- V %*% D
    return(X)
  } else {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), k = d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    return(X)
  }
}

#' Function to perform graph generalised adjacency spectral embedding
#'
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to choose (if d is given, this number is ignored)
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#'
#' @return A list containing a matrix with n rows and d columns representing the estimated latent positions, and the estimated
#' indefinite orthogonal projection matrix
#'
#' @references Rubin-Delanchy, Patrick, et al. "A statistical interpretation of spectral embedding:
#' the generalised random dot product graph." arXiv preprint arXiv:1709.05506 (2017).
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' @export
g.ase <- function(A, d = NA, d.max = ncol(A), diag.augment = T, elbow = 1) {
  if(is.na(d)) {
    if(diag.augment & sum(abs(Matrix::diag(A))) == 0) {
      deg = Matrix::colSums(A)
      n = ncol(A)
      diag(A) = deg / (n-1)
    }
    eigv <- rARPACK::eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eigv$values), decreasing = TRUE)
    #d = getElbow_GMM(vals)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eigv$values) >= vals[d])
    X <- eigv$vectors[,selected.eigs, drop = F] %*% diag(sqrt(abs(eigv$values[selected.eigs])), nrow = d)
    D <- sign(eigv$values[selected.eigs])
  } else{
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), k =  d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    D <- sign(eig$values)
  }
  return(list(X=X, D=D))
}

#'
#' Obtain expected adjacency matrix from generalised RDPG
#' @param g_ase output of g.ase
#'
#' @return A matrix with expected adjacency
#' @export
P_from_g.ase <- function(g_ase) {
  g_ase$X %*% diag(g_ase$D, ncol = ncol(g_ase$X)) %*% t(g_ase$X)
}

#'
#' Function to compute the graph unscaled adjacency spectral embedding (top eigenvectors)
#'
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NULL, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to choose (if d is given, this number is ignored)
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#'
#' @return A list containing a matrix with n rows and d columns representing the estimated latent positions, and the estimated
#' indefinite orthogonal projection matrix
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' @export
eig_embedding <- function(A, d = NA, d.max = ncol(A), diag.augment = FALSE, elbow = 1) {
  n <- ncol(A)
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = Matrix::colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.na(d)) {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eig$values), decreasing = TRUE)#[1:sqrt(n)]
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    eig <- eig$vectors[,selected.eigs, drop = FALSE]
  }else {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), k = d)$vectors
  }
  return(eig)
}

#'
#' Function to estimated the score matrices of a list of graphs given the common invariant subspace V
#'
#' @param Adj_list list of adjacency matrices, of size n x n
#' @param V common invariant subspace. A matrix of size n x d.
#' @return A list containing the score matrices
#'
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
#' @export
project_networks <- function(Adj_list, V) {
  lapply(Adj_list, function(A) Matrix::crossprod(Matrix::crossprod(A, V), V))
}




graph_distance <- function(A, B) {
  return(sum(abs(A-B)))
}

cluster_positions <- function(V, K) {
  clus <- Mclust(V, G = K)
  return(clus$classification)
}


sample_from_P<-function (P) 
{
  n = ncol(P)
  A = Matrix::Matrix(0, n, n)
  A[upper.tri(A)] <- 1 * (runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A = A + Matrix::t(A)
  return(A)
}