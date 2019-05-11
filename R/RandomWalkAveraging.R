#' @title Random walk averaging
#'
#' @description
#' \code{RandomWalkAveraging}
#'
#' @details
#'
#' @param X
#' @param Z
#' @param genes.use
#' @param distance.method
#' @param neighbor.size
#' @param neighbor.scale
#' @param gamma
#' @param tau
#' @param alpha
#' @param do.pca
#' @param pca.param
#' @param to.diffuse
#' @param diffuse.thresh
#'
#' @return
#'
#' @export
#'
#' @import expm
#'
#' @author Xiaohan Chen <1356957916@qq.com>
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
RandomWalkAveraging <- function(X,
                                Z,
                                genes.use = NULL,
                                distance.method = c("euclidean", "manhattan"),
                                neighbor.size = 15, neighbor.scale = 0.3,
                                tau = 15,
                                alpha = 0.5,
                                do.pca = TRUE, pca.params = 15,
                                to.diffuse = NULL,
                                diffuse.thresh = NULL)
{
  distance.method = match.arg(distance.method)

  if (!is.null(genes.use)) {
    X <- X[genes.use, ]
  }

  if (do.pca) {
    X_t <- t(X)
    pca_result <- rsvd::rpca(A = X_t, k = pca.params)
    L <- as.matrix(dist(pca_result$x, method = distance.method))
  } else {
    L <- as.matrix(dist(t(X), method = distance.method))
  }

  SS <- apply(L, 1, function(x){
    s <- sort(abs(x))[neighbor.size+1]
    xs <- sd(x)
    y <- exp(-(x/xs)**2)
    i <- abs(x) > s
    y[i] <- exp(-(x[i]/(neighbor.scale*xs))**2)
    return(y)
  })
  diag(SS) <- 0
  SS <- t(SS) + SS
  SS <- SS / rowSums(SS)

  # diffuse tau times
  if (!is.infinite(tau)) {
    D <- SS %^% tau
  }else{
    D <- (1-alpha) * MASS::ginv(diag(rep(1, ncol(SS))) - alpha * t(SS))
    D <- t(D)
  }
  rownames(D) <- rownames(Z)
  colnames(D) <- rownames(Z)

  # diffuse identity
  Z.max <- NULL
  noto.diffuse <- NULL
  if (!is.null(diffuse.thresh)) {
    if (is.null(Z.max)) {
      Z.max <- apply(Z, 1, max)
    }
    to.diffuse <- Z.max >= diffuse.thresh
  }
  if (is.logical(to.diffuse)) {
    to.diffuse <- rownames(Z)[to.diffuse]
  }
  if (is.null(to.diffuse) ||
      (length(to.diffuse) == nrow(Z)) ||
      (length(to.diffuse) == 0)) {
    to.diffuse <- rownames(Z)
    noto.diffuse <- rownames(Z)
  }else{
    noto.diffuse <- setdiff(rownames(Z), to.diffuse)
  }

  R <- Z
  R[noto.diffuse, ] <- t(D[to.diffuse, noto.diffuse]) %*% Z[to.diffuse, ]
  R <- R/rowSums(R)

  # results
  results <- list(celltype = factor(apply(R, 1, function(x) names(x)[which.max(x)])),
                  score = apply(R, 1, max),
                  score.matrix = t(R))
  return(results)

}
