#'
#' Function to perform degree-corrected multiple adjacency spectral embedding
#'
#' @param Adj_list a list of adjacency matrices with the same size n x n
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically via scree plot.
#' @param d_vec vector of individual embedding dimensions for each graph. If NA, this is the same as d.
#' @param row.normalization parameter that indicates the option to normalize the rows of each embedding. Default is "2" (currently, this is the only one implemented).
#' @param diag.augment logical. When true, the diagonal augmentation method for ASE is performed.
#' @param elbow_graph number of elbow selected in Zhu & Ghodsi method for the scree plot of each individual graph eigenvalues.
#' @param elbow_dcmase number of elbow selected in Zhu & Ghodsi method for the scree plot of the singular values of the concatenated spectral embeddings of MASE.
#' @param show.scree.results when TRUE, the histogram of the estimated d for each graph, and the scree plot of the singular values of  the graph is shown if d is not specified.
#' @param par whether to run in parallel (TRUE/FALSE)
#' @param numpar number of clusters for parallel implmentation
#'
#' @return A matrix of size n x d, with the degree corrected joint embedding.
#'
#' @references
#'
#' @author Jes\'us Arroyo <jarroyo@tamu.edu>
#' @export
dcmase <- function(Adj_list, d = NA, d_vec = NA,
                    row.normalization = "2", diag.augment = FALSE,
                 elbow_graph = 1, elbow_dcmase = 2,
                 show.scree.results = FALSE,
                 par = FALSE, numpar = 12) {
  if(is.na(d_vec)) {
    d_vec = rep(d, length(Adj_list))
  }
  # running in parallel--------------------------------
  if(par) {
    stop("Not yet implemented")
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
    if(row.normalization=="2") {
      latpos.list <- lapply(1:length(Adj_list), function(i) {
        U <- eig_embedding(Adj_list[[i]], d = d_vec[i], diag.augment = diag.augment, elbow = elbow_graph)
        unorms <- apply(U, 1, function(x) sqrt(sum(x^2)))
        U <- U / ifelse(unorms==0, 1, unorms)
        U
      })
    }else{
      stop("Not yet implemented")
    }
  }
  V_all  <- Reduce(cbind, latpos.list)
  jointsvd <- svd(V_all)
  if(is.na(d)) {
    if(show.scree.results) {
      hist(sapply(latpos.list, ncol), main = "Estimated d for each graph")
    }
    d = getElbows(jointsvd$d, plot = show.scree.results)[elbow_dcmase]
  }
  V = jointsvd$u[, 1:d, drop = FALSE]
  #print(plot(jointsvd$d))
  return(V)
}
