#' @title Multinomial classifier
#'
#' @description
#' \code{MultinomialClassifier}
#'
#' @details
#'
#' @param X
#' @param catalogue.profile
#' @param catalogue.celltype
#' @param catalogue.prior
#' @param genes.use
#' @param cells.use
#' @param BPPARAM
#'
#' @return
#'
#' @export
#'
#' @import BiocParallel
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
MultinomialClassifier <- function(X,
                                  catalogue.profile,
                                  catalogue.celltype,
                                  catalogue.prior = NULL,
                                  X.cluster = NULL,
                                  genes.use = NULL,
                                  cells.use = NULL,
                                  BPPARAM = SerialParam()){

  if (is.vector(X)){
    X <- as.matrix(X, ncol=1)
  }

  if (is.null(genes.use)) {
    genes.use <- rownames(X)
  }else {
    genes.use <- genes.use
  }
  genes.use <- intersect(genes.use,
                         intersect(rownames(X),
                                   rownames(catalogue.profile)))

  if (is.null(cells.use)) {
    cells.use <- colnames(X)
  }else {
    cells.use <- cells.use
  }

  P <- catalogue.profile[genes.use, ]
  P <- apply(P, 2, function(z) z/sum(z))
  X <- X[genes.use, cells.use]
  if (!is.null(X.cluster)) {
    X.cluster <- X.cluster[cells.use]
    X <- tapply(1:ncol(X), X.cluster, function(i) rowMeans(X[,i]))
    X <- do.call("cbind", X)
  }

  # likelihood
  L <- bplapply(1:ncol(X), function(i, X, P, MultiLikelihood) MultiLikelihood(X[,i], P, use.log = TRUE, BPPARAM = BiocParallel::SerialParam()), X, P, MultiLikelihood, BPPARAM = BPPARAM)
  L <- do.call("rbind", L)
  L <- t(L)
  L <- bplapply(1:ncol(L), function(i, L) L[,i] - max(L[,i]), L, BPPARAM = BPPARAM)
  L <- do.call("rbind", L)
  L <- t(L)
  colnames(L) <- cells.use

  # posterior for typical cells
  if (is.null(catalogue.prior)){
    P <- bplapply(1:ncol(L), function(i, L, post) post(L[,i], NULL, TRUE), L, Posterior, BPPARAM = BPPARAM)
  }else{
    P <- bplapply(1:ncol(L), function(i, L, post, prior) post(L[,i], prior, TRUE), L, Posterior, catalogue.prior, BPPARAM = BPPARAM)
  }
  P <- do.call("rbind", P)
  P <- data.frame(t(P))
  rownames(P) <- colnames(catalogue.profile)
  colnames(P) <- cells.use
  P[is.na(P)] <- 0

  # posterior of cell types
  names(catalogue.celltype) <- colnames(catalogue.profile)
  P$type <- bplapply(rownames(P), function(x, celltype) celltype[x], catalogue.celltype, BPPARAM = BPPARAM)
  if (is.list(P$type)) {
    P$type <- unlist(P$type)
  }
  Z <- aggregate(. ~ type, P, sum)
  n <- Z[, which(names(Z) %in% "type")]
  Z <- Z[, -which(names(Z) %in% "type")]
  rownames(Z) <- n
  colnames(Z) <- cells.use

  # assign a cell to the type, of which the posterior is max
  C <- bplapply(1:ncol(Z), function(i, Z) which.max(Z[,i]), Z, BPPARAM = BPPARAM)
  C <- unlist(C)
  names(C) <- cells.use
  S <- bplapply(1:ncol(Z), function(i, Z) Z[which.max(Z[,i]),i], Z, BPPARAM = BPPARAM)
  S <- unlist(S)
  names(S) <- cells.use

  # results
  if (is.null(X.cluster)) {
    results <- list(celltype = as.factor(n[C]),
                    score = S,
                    score.matrix = Z)
  }else{
    results <- list(celltype = as.factor(n[C[X.cluster]]),
                    score = S[X.cluster],
                    score.matrix = Z[, X.cluster])
  }

  return(results)
}
