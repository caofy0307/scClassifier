#' @title Dirichlet multinomial classifier
#'
#' @description
#' \code{DirichletMultinomialClassifier} computes the predictive probability for the given cells following
#' Multinomial distribution model with the Dirichlet conjugate prior.
#'
#' @details
#'
#' @param X
#' @param catalogue.profile
#' @param catalogue.celltype
#' @param catalogue.prior
#' @param X.cluster
#' @param genes.use
#' @param cells.use
#' @param do.training
#' @param use.log
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
DirichletMultinomialClassifier <- function(X,
                                           catalogue.profile,
                                           catalogue.celltype,
                                           catalogue.prior = NULL,
                                           X.cluster = NULL,
                                           genes.use = NULL,
                                           cells.use = NULL,
                                           do.training = T,
                                           BPPARAM = SerialParam()){
  # normalize genes
  if (is.null(genes.use)){
    genes.use <- rownames(X)
  }
  genes.use <- intersect(genes.use,
                         intersect(rownames(X), rownames(catalogue.profile)))

  # normalize cells
  if (is.null(cells.use)){
    cells.use = colnames(X)
  }

  if (is.vector(X)){
    X <- as.matrix(X, ncol = 1)
  }

  X <- X[genes.use, cells.use, drop = F]
  if (!is.null(X.cluster)) {
    X.cluster <- X.cluster[cells.use]
    X <- tapply(1:ncol(X), X.cluster, function(i) rowMeans(X[,i]))
    X <- do.call("cbind", X)
  }

  # normalize parameters
  catalogue.profile <- catalogue.profile[genes.use, ]
  if (!is.factor(catalogue.celltype)){
    catalogue.celltype <- as.factor(catalogue.celltype)
  }

  if (do.training){
    catalogue.profile <- dirichlet_fit(catalogue.profile, catalogue.celltype)
    catalogue.celltype <- factor(levels(catalogue.celltype))
  }
  catalogue.profile[catalogue.profile<=0] <- 1e-2 * min(catalogue.profile[catalogue.profile>0])
  names(catalogue.celltype) <- colnames(catalogue.profile)

  # likelihood
  L <- bplapply(1:ncol(X),
                function(i, X, p, f) {
                  apply(p, 2, function(q) f(X[,i], q, TRUE))
                },
                X,
                catalogue.profile,
                DirichletMultinomialClassifier1,
                BPPARAM = BPPARAM)
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
