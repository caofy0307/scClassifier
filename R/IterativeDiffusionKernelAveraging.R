#' @title Iterative diffusion kernel averaging
#'
#' @description
#' \code{IterativeDiffusionKernelAveraging}
#'
#' @details
#'
#' @param X
#' @param Z
#' @param genes.use
#' @param neighbor.size
#' @param neighbor.scale
#' @param gamma
#' @param tau
#' @param do.pca
#' @param pca.param
#' @param thresh
#' @param pct.good
#'
#' @return
#'
#'
#' @author Xiaohan Chen <1356957916@qq.com>
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
IterativeDiffusionKernelAveraging <- function(X,
                                              Z,
                                              genes.use = NULL,
                                              neighbor.size = 15, neighbor.scale = 0.3,
                                              gamma = 0.01, tau = 15,
                                              do.pca = TRUE, pca.params = 15,
                                              thresh = 0.8,
                                              pct.good = 0.99)
{

  num.cells <- ncol(X)
  Z <- apply(Z, 1, function(x) x/sum(x))
  Z <- t(Z)
  X.z <- apply(Z, 1, max)
  num.good <- sum(X.z >= thresh)
  do.iterative = FALSE
  if (num.good < pct.good * num.cells){
    do.iterative = TRUE
  }

  if (do.iterative){
    good <- X.z >= thresh
    Y <- Z
    Y[!good, ] <- 0
    Y <- DiffusionKernelAveraging(X = X,
                                  Z = Y,
                                  genes.use = genes.use,
                                  neighbor.size = neighbor.size,
                                  neighbor.scale = neighbor.scale,
                                  gamma = gamma,
                                  tau = tau,
                                  do.pca = do.pca,
                                  pca.params = pca.params)
    Y <- t(Y$score.matrix)
    Y <- apply(Y, 1, function(x) x/sum(x))
    Y <- t(Y)
    Z[!good, ] <- Y[!good, ]

    X.z <- apply(Z, 1, max)
    num.good <- sum(X.z >= thresh)
    if (num.good < pct.good * num.cells){
      do.iterative = TRUE
    }
  }

  # results
  results <- list(celltype = factor(apply(Z, 1, function(x) names(x)[which.max(x)])),
                  score = apply(Z, 1, max),
                  score.matrix = t(Z))
  return(results)

}
