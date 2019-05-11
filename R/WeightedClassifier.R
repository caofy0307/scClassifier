#' @title Weighted classifier
#'
#' @description
#' \code{WeightedClassifier}
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
#' @param debatch
#' @param do.training
#' @param latent.classifier
#' @param distance.method
#' @param neighbor.size
#' @param neighbor.scale
#' @param gamma
#' @param tau
#' @param BPPARM
#' @param verbose
#'
#' @return
#'
#' @import BiocParallel
#' @import scran
#'
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
WeightedClassifier <- function(X,
                               catalogue.profile,
                               catalogue.celltype,
                               catalogue.prior = NULL,
                               X.cluster = NULL,
                               genes.use = NULL,
                               cells.use = NULL,
                               debatch = FALSE,
                               do.training = FALSE,
                               latent.classifier = c(DirichletMultinomialClassifier, MultinomialClassifier),
                               distance.method = c("euclidean", "KL", "manhattan"),
                               neighbor.size = 15,
                               neighbor.scale = 0.3,
                               gamma = 0.01,
                               tau = 5,
                               do.pca = TRUE,
                               pca.params = 15,
                               BPPARAM = SnowParam(),
                               verbose = TRUE)
{
  distance.method <- match.arg(distance.method)
  latent.classifier <- match.arg(latent.classifier)
  neighbor.scale <- min(neighbor.scale, 1)
  neighbor.scale <- max(neighbor.scale, 0)

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

  X <- X[genes.use, cells.use]
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

  # invoking latent classifier
  Z <- do.call(latent.classifier,
               list(X = X,
                    catalogue.profile = catalogue.profile,
                    catalogue.celltype = catalogue.celltype,
                    catalogue.prior = catalogue.prior,
                    do.training = do.training,
                    BPPARAM = BPPARAM))

  # invoking diffusion kernel averaging
  L <- DiffusionKernelAveraging(X = X,
                                Z = t(Z$score.matrix),
                                neighbor.size = neighbor.size,
                                neighbor.scale = neighbor.scale,
                                tau = tau,
                                gamma = gamma,
                                do.pca = do.pca,
                                pca.params = pca.params)
  L <- t(L)

  # wrap results
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
