# Spectral methods

ase <- function (A, d = NA, d.max = sqrt(ncol(A)), diag.augment = TRUE, 
                 elbow = 1) {
  if (diag.augment & sum(abs(Matrix::diag(A))) == 0) {
    deg = Matrix::colSums(A)
    n = ncol(A)
    diag(A) = deg/(n - 1)
  }
  if (is.na(d)) {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x = abs(eig$values), decreasing = TRUE)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    V <- eig$vectors[, selected.eigs, drop = F]
    D <- diag(sqrt(abs(eig$values[selected.eigs])), nrow = d)
    X <- V %*% D
    return(X)
  }
  else {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), k = d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    return(X)
  }
}

g.ase <- function (A, d = NA, d.max = ncol(A), diag.augment = T, elbow = 1) 
{
  if (is.na(d)) {
    if (diag.augment & sum(abs(Matrix::diag(A))) == 0) {
      deg = Matrix::colSums(A)
      n = ncol(A)
      diag(A) = deg/(n - 1)
    }
    eigv <- rARPACK::eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x = abs(eigv$values), decreasing = TRUE)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eigv$values) >= vals[d])
    X <- eigv$vectors[, selected.eigs, drop = F] %*% diag(sqrt(abs(eigv$values[selected.eigs])), 
                                                          nrow = d)
    D <- sign(eigv$values[selected.eigs])
  }
  else {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), k = d)
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    D <- sign(eig$values)
  }
  return(list(X = X, D = D))
}
eig_embedding <- function (A, d = NA, d.max = ncol(A), diag.augment = FALSE, elbow = 1) 
{
  n <- ncol(A)
  if (diag.augment & sum(abs(Matrix::diag(A))) == 0) {
    deg = Matrix::colSums(A)
    n = ncol(A)
    Matrix::diag(A) = deg/(n - 1)
  }
  if (is.na(d)) {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x = abs(eig$values), decreasing = TRUE)
    d = getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eig$values) >= vals[d])
    eig <- eig$vectors[, selected.eigs, drop = FALSE]
  }
  else {
    eig <- rARPACK::eigs(as(A, "dgeMatrix"), k = d)$vectors
  }
  return(eig)
}