#' @title Metacell
#'
#' @description
#' \code{Metacell}
#'
#' @details
#'
#' @param X
#' @param y
#' @param metacell.size
#' @param catalogue.prior
#' @param metacells.per.celltype
#'
#' @return
#'
#' @export
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'

Metacell <- function(X, y, metacell.size = 100, metacells.per.celltype = 30) {
  if (!is.factor(y)) {
    y = factor(y, levels = unique(y))
  }

  n.celltypes = nlevels(y)
  Z = sapply(levels(y), function(i) {
    x = X[, y==i, drop = F]
    x = x[, sample(colnames(x), metacell.size * metacells.per.celltype, replace = T)]
    z = sapply(1:metacells.per.celltype, function(j) {
      rowSums(x[, ((j-1)*metacell.size+1):(j*metacell.size)])
    })
    colnames(z) = paste(i, 1:metacells.per.celltype, sep="_")
    rownames(z) = rownames(X)
    return(z)
  }, simplify = F)
  Z = do.call("cbind", Z)

  res = list(metacells = Z, celltypes = rep(levels(y), each = metacells.per.celltype))
  return(res)
}
