#' @title SingleRDEGenes
#'
#' @description
#' \code{SingleRDEGenes}
#'
#' @details
#'
#' @param X
#' @param y
#' @param size
#'
#' @return
#'
#' @export
#'
#' @import SingleR
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
SingleRDEGenes <- function(X, y, size = 500){

  X <- medianMatrix(X, y)

  n = round(500 * (2/3)^(log2(c(ncol(X)))))
  genes.use = unique(unlist(unlist(lapply(1:ncol(X),
                                          function(j) {
                                            lapply(1:ncol(X), function(i) {
                                              s = sort(X[, j] - X[, i], decreasing = T)
                                              s = s[s > 0]
                                              names(s)[1:min(n, length(s))]
                                            })
                                          }))))[-1]
  return(genes.use[1:min(length(genes.use),size)])
}
