#' @title SeuratDEGenes
#'
#' @description
#' \code{SeuratDEGenes}
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
#' @param pc.pval
#' @param max.per.pc
#' @param dims.use
#' @param k.param
#' @param resolution
#' @param test.use
#' @param min.pct
#' @param min.diff.pct
#' @param logfc.threshold
#' @param max.per.cluster
#' @param de.pval
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
SeuratDEGenes <- function(X,
                          x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 0.15,
                          num.replicate = 100, num.cores = detectCores() - 1, do.par = TRUE,
                          pcs.use = 1:5, pc.pval = 1e-10, max.per.pc = 200,
                          dims.use = 1:10, k.param = 30, resolution = 0.8,
                          test.use = "wilcox", min.pct = 0.2, min.diff.pct = -Inf, logfc.threshold = 0.2, max.per.cluster = 30, de.pval = 0.05,
                          size = 500){

  sc <- CreateSeuratObject(X)
  sc <- NormalizeData(sc)

  message("Find highly variable genes...")
  sc <- FindVariableGenes(sc,
                          x.low.cutoff = x.low.cutoff,
                          x.high.cutoff = x.high.cutoff,
                          y.cutoff = y.cutoff,
                          do.plot = F)
  var.genes <- sc@var.genes
  sc <- ScaleData(sc, display.progress = F)

  message("Find PCA significant genes...")
  sc <- RunPCA(sc, do.print = F)
  # JackStraw
  sc <- JackStraw(sc, num.replicate = num.replicate, num.cores = num.cores, do.par = do.par, display.progress = F)
  # # projection
  sc <- ProjectPCA(sc, do.print = F)
  # genes that are significantly correlated with PCs
  sig.genes <- PCASigGenes(sc, pcs.use = pcs.use, pval.cut = pc.pval, max.per.pc = max.per.pc)
  # Run PCA again but use significant genes
  sc <- RunPCA(sc, pc.genes = sig.genes, do.print = F)

  message("Find differentially expressed genes...")
  sc <- FindClusters(object = sc,
                     reduction.type = "pca",
                     dims.use = dims.use,
                     k.param = k.param,
                     resolution = resolution,
                     print.output = 0,
                     save.SNN = FALSE,
                     force.recalc = TRUE)


  markers <- FindAllMarkers(object = sc,
                            genes.use = unique(c(var.genes, sig.genes)),
                            test.use = test.use,
                            only.pos = TRUE,
                            min.pct = min.pct,
                            min.diff.pct = min.diff.pct,
                            logfc.threshold = logfc.threshold)
  markers %>% group_by(cluster) %>% top_n(max.per.cluster, avg_logFC) -> topMarkers
  topMarkers <- topMarkers[topMarkers$p_val_adj<de.pval,] %>% arrange(desc(avg_logFC))
  topGenes <- unique(topMarkers$gene)
  de.genes <- topGenes

  genes.use <- de.genes
  return(genes.use[1:min(size, length(genes.use))])
}
