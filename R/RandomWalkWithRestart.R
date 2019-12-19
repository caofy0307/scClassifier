#' @title Random walking with restart algorithm
#'
#' @description
#' \code{RandomWalkWithRestart}
#'
#' @details
#'
#' @param X
#' @param Z
#' @param genes.use
#' @param neighbor.size
#' @param do.pca
#' @param pca.param
#' @param tau
#' @param alpha
#' @param BPPARAM
#'
#' @return
#'
#' @export
#'
#' @import BiocParallel
#' @import rsvd
#' @import expm
#'
#' @author Xiaohan Chen <1356957916@qq.com>
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'

RandomWalkWithRestart <- function(X,
                                  Z,
                                  genes.use = NULL,
                                  neighbor.size = 15,
                                  do.pca = TRUE, pca.params = 15,
                                  tau = Inf,
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
    y <- exp(-(x/s)**2)
    i <- abs(x) > s
    y[i] <- 0
    return(y)
  })
  diag(A) <- 0

  A <- apply(A, 2, function(x) x/sum(x))
  if (is.infinite(tau)) {
    D <- (1-alpha) * MASS::ginv(diag(rep(1, ncol(A))) - alpha * A)
  }else{
    D <- alpha * A + diag(rep(1-alpha, ncol(A)))
    D <- D %^% tau
  }
  P <- D %*% Z

  results <- list(celltype = factor(apply(P, 1, function(x) names(x)[which.max(x)])),
                  score = apply(P, 1, max),
                  score.matrix = t(P))

  return(results)

}
