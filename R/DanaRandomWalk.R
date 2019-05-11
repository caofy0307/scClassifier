#' @title Dana Pe'er's random walking algorithm
#'
#' @description
#' \code{DanaRandomWalk}
#'
#' @details
#'
#' @param X
#' @param Z
#' @param labeled
#' @param unlabeled
#' @param genes.use
#' @param neighbor.size
#' @param do.pca
#' @param pca.param
#' @param thresh
#' @param BPPARAM
#'
#' @return
#'
#' @export
#'
#' @import BiocParallel
#' @import rsvd
#'
#' @author Xiaohan Chen <1356957916@qq.com>
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'

DanaRandomWalk <- function(X,
                           Z,
                           labeled = NULL,
                           unlabeled = NULL,
                           genes.use = NULL,
                           neighbor.size = 15,
                           do.pca = TRUE, pca.params = 15,
                           thresh = 0.9,
                           BPPARAM = SerialParam())
{

  if (!is.null(genes.use)) {
    X <- X[genes.use, ]
  }

  if (!is.null(thresh)) {
    Z.max <- apply(Z, 1, max)
    labeled <- rownames(Z)[Z.max >= thresh]
    unlabeled <- setdiff(rownames(Z), labeled)
  }else if (is.null(labeled) && is.null(unlabeled)) {
    labeled <- rownames(Z)
    unlabeled <- rownames(Z)
  }else if (!is.null(labeled)) {
    unlabeled <- setdiff(rownames(Z), labeled)
  }else if (!is.null(unlabeled)) {
    labeled <- setdiff(rownames(Z), unlabeled)
  }
  if ((length(labeled) == nrow(Z)) || (length(labeled) == 0)) {
    labeled <- rownames(Z)
    unlabeled <- rownames(Z)
  }

  if (do.pca) {
    Y <- rpca(t(X), k = pca.params)
    X <- t(Y$x)
    rownames(X) <- paste("PC", 1:pca.params, sep = "")
  }

  AdaptiveAnisotropicKernel <- function(d.xy, s.x, s.y) {
    s.xy <- s.x + s.y
    k <- exp(-0.5 * d.xy ** 2 / s.xy)
    k / sqrt(2 * pi * s.xy)
  }

  # kernel matrix
  L <- as.matrix(dist(t(X), method = "euclidean"))
  S <- bplapply(1:ncol(L), function(i, X, k, kf) {
    s.x <- sort(X[i,])[k+1]
    sd.x <- sd(X[i,])
    Y <- sapply(i:ncol(X), function(j){
      s.y <- sort(X[j,])[k+1]
      sd.y <- sd(X[j,])
      if (X[i,j] < s.x) {
        return(kf(X[i,j],sd.x,sd.y))
      }else{
        return(0)
      }
    })
  }, L, neighbor.size, AdaptiveAnisotropicKernel,
  BPPARAM = BPPARAM)
  for (i in 1:ncol(L)) {
    L[i, i:ncol(L)] <- S[[i]]
  }
  L <- L + t(L)
  diag(L) <- 0.5 * diag(L)
  rownames(L) <- colnames(X)
  colnames(L) <- colnames(X)

  # adjacency matrix
  A <- bplapply(1:nrow(L), function(i, L, size){
    x <- L[i, ]
    names(x) <- colnames(L)
    if (size < length(x)) {
      s <- names(x[-i])[order(x[-i], decreasing = T)]
      s <- s[(size+1):length(s)]
      x[s] <- 0
    }
    return(x)
  }, L, neighbor.size,
  BPPARAM = BPPARAM)
  A <- do.call("rbind", A)
  rownames(A) <- rownames(L)



  # markovian probability
  A <- A / rowSums(A)
  Q <- A[unlabeled, unlabeled]
  R <- A[unlabeled, labeled]
  # O <- solve(diag(rep(1, length(unlabeled))) - Q)
  O <- MASS::ginv(diag(rep(1, length(unlabeled))) - Q)
  BC <- O %*% R
  PP <- BC %*% Z[labeled, ]

  P <- Z
  P[unlabeled, ] <- PP

  results <- list(celltype = factor(apply(P, 1, function(x) names(x)[which.max(x)])),
                  score = apply(P, 1, max),
                  score.matrix = t(P))

  return(results)

}
