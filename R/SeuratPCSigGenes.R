#' @title SeuratPCSigGenes
#'
#' @description
#' \code{SeuratPCSigGenes}
#'
#' @details
#'
#' @param X
#' @param x.low.cutoff
#' @param x.high.cutoff
#' @param y.cutoff
#' @param num.replicate
#' @param num.cores
#' @param do.par
#' @param pcs.use
#' @param pval.cut
#' @param max.per.pc
#' @param size
#'
#' @return
#'
#' @export
#'
#' @import Seurat
#' @import parallel
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
SeuratPCSigGenes <- function(X,
                             x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.15,
                             num.replicate = 100, num.cores = detectCores() - 1, do.par = TRUE,
                             pcs.use = 1:5, pval.cut = 1e-10, max.per.pc = 200,
                             size = 500){

  sc <- CreateSeuratObject(X)
  sc <- NormalizeData(sc)
  sc <- FindVariableGenes(sc,
                          x.low.cutoff = x.low.cutoff,
                          x.high.cutoff = x.high.cutoff,
                          y.cutoff = y.cutoff)
  sc <- ScaleData(sc)

  sc <- RunPCA(sc, do.print = F)
  # JackStraw
  sc <- JackStraw(sc, num.replicate = num.replicate, num.cores = num.cores, do.par = do.par)
  # # projection
  sc <- ProjectPCA(sc, do.print = F)
  # genes that are significantly correlated with PCs
  sig.genes <- PCASigGenes(sc, pcs.use = pcs.use, pval.cut = pval.cut, max.per.pc = max.per.pc)

  genes.use <- sig.genes
  return(genes.use[1:min(size, length(genes.use))])
}
