#' @title Implement Jaitin et al. classifier
#'
#' @description
#' \code{JaitinClassifier}
#'
#' @details
#'
#' @param X
#' @param catalogue.profile
#' @param catalogue.celltype
#' @param X.cluster
#' @param genes.use
#' @param BPPARAM
#'
#' @return
#'
#' @export
#'
#' @import BiocParallel
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'

JaitinClassifier <- function(X,
                             catalogue.profile,
                             catalogue.celltype,
                             X.cluster = NULL,
                             genes.use = NULL,
                             BPPARAM = SnowParam())
{
  if (is.vector(X)) {
    X <- as.matrix(X, ncol = 1)
  }

  if (is.null(genes.use)) {
    genes.use <- rownames(X)
  }
  genes.use <- intersect(genes.use, rownames(X))
  genes.use <- intersect(genes.use, rownames(catalogue.profile))

  if (!is.factor(catalogue.celltype)) {
    catalogue.celltype <- factor(catalogue.celltype)
  }
  names(catalogue.celltype) <- colnames(catalogue.profile)

  X <- X[genes.use, ]
  catalogue.profile <- catalogue.profile[genes.use, ]

  # this function is to perform classification for one cluster
  classifyone <- function(x){
    # cluster model
    if (is.matrix(x)) {
      m <- rowSums(x)/sum(x)
    }else {
      m <- x/sum(x)
    }

    # compute likelihood of typical cells against the cluster model
    loglik <- sapply(colnames(catalogue.profile), function(p) MultiLikelihood(catalogue.profile[, p], m, use.log = TRUE, BPPARAM = BPPARAM))

    # compute posterior probability
    pp <- Posterior(loglik-max(loglik), use.log = TRUE)

    # pick the typical cell achiving the maximal likelihood
    typical.cell <- colnames(catalogue.profile)[which.max(loglik)]
    cell.type <- catalogue.celltype[typical.cell]
    loglik2 <- sapply(levels(catalogue.celltype), function(l) max(loglik[catalogue.celltype==l]))
    pp2 <- sapply(levels(catalogue.celltype), function(l) max(pp[catalogue.celltype==l]))

    return(list(celltype = cell.type, score = loglik[which.max(loglik)], loglik = loglik2, post = pp2))
  }

  if (is.null(X.cluster)) {
    X.type <- classifyone(X)
  }else {
    if (!is.factor(X.cluster)) {
      X.cluster <- factor(X.cluster)
    }
    X.type <- lapply(levels(X.cluster),
                     function(l) classifyone(X[, X.cluster==l]))
  }

  return(X.type)
}
