#' @title SeuratHVGenes
#'
#' @description
#' \code{SeuratHVGenes}
#'
#' @details
#'
#' @param X
#' @param x.low.cutoff
#' @param x.high.cutoff
#' @param y.cutoff
#' @param size
#'
#' @return
#'
#' @export
#'
#' @import Seurat
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
SeuratHVGenes <- function(X, x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.15, size = 500){

  sc <- CreateSeuratObject(X)
  sc <- NormalizeData(sc)
  sc <- FindVariableGenes(sc,
                          x.low.cutoff = x.low.cutoff,
                          x.high.cutoff = x.high.cutoff,
                          y.cutoff = y.cutoff)
  genes.use <- sc@var.genes
  return(genes.use[1:min(size, length(genes.use))])
}
