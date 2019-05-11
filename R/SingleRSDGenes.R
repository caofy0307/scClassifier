#' @title SingleRSDGenes
#'
#' @description
#' \code{SingleRSDGenes}
#'
#' @details
#'
#' @param X
#' @param y
#' @param sd.thresh
#' @param size
#'
#' @return
#'
#' @export
#'
#' @import matrixStats
#' @import SingleR
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
SingleRSDGenes <- function(X, y, sd.thresh = 1, size = 500){

  X <- medianMatrix(X, y)
  sd = rowSds(as.matrix(X))
  genes.use = rownames(X)[sd > sd.thresh]
  sd = sd[sd > sd.thresh]
  genes.use = genes.use[order(sd, decreasing = T)]
  return(genes.use[1:min(size,length(genes.use))])
}
