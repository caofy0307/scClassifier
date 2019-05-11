#' @title Dana Pe'er's random walking algorithm
#'
#' @description
#' \code{DanaRandomWalk}
#'
#' @details
#'
#' @param X
#' @param Z
#' @param genes.use
#' @param neighbor.size
#' @param do.pca
#' @param pca.param
#' @param alpha
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

AkalinRandomWalk <- function(X,
                             Z,
                             genes.use = NULL,
                             neighbor.size = 15,
                             do.pca = TRUE, pca.params = 15,
                             alpha = 0.1,
                             BPPARAM = SnowParam())
{

  if (!is.null(genes.use)) {
    X <- X[genes.use, ]
  }

  if (do.pca) {
    Y <- rpca(t(X), k = pca.params)
    X <- t(Y$x)
    rownames(X) <- paste("PC", 1:pca.params, sep = "")
  }

  # kernel matrix
  L <- as.matrix(dist(t(X), method = "euclidean"))

  # adjacency matrix
  A <- apply(L, 1, function(x){
    s <- sort(abs(x))[neighbor.size+1]
    y <- exp(-(x/sd(x))**2)
    i <- abs(x) > s
    y[i] <- 0
    return(y)
  })
  diag(A) <- 0

  A <- apply(A, 2, function(x) x/sum(x))
  P <- (1-alpha) * MASS::ginv(diag(rep(1, ncol(A))) - alpha * A) %*% Z

  results <- list(celltype = factor(apply(P, 1, function(x) names(x)[which.max(x)])),
                  score = apply(P, 1, max),
                  score.matrix = t(P))

  return(results)

}
