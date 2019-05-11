#' @title Pseudo Single Cell Generator
#'
#' @description
#' \code{PseudoSingleCellGenerator}
#'
#' @details
#'
#' @param X
#' @param catalogue.profile
#' @param catalogue.celltype
#' @param genes.use
#' @param sizes A numerical value or a vector. A numerical value indicates the total size of simulating cells. A vector gives the simulation size for each cell type.
#' @param method
#' @param lambda
#' @param epsilon
#' @param maxit
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

PseudoSingleCellGenerator <- function(X,
                                      catalogue.profile,
                                      catalogue.celltype,
                                      genes.use = NULL,
                                      sizes = 1000,
                                      simulator = c("scSimulator", "Multinomial"),
                                      method = c("bi", "qt"),
                                      lambda = 0.01,
                                      epsilon = 0.01,
                                      maxit = 100,
                                      BPPARAM = SnowParam())
{
  simulator = match.arg(simulator)
  simulator = tolower(simulator)
  if (!(simulator %in% c("scsimulator", "multinomial"))) {
    stop("Try to use a invalid simulation method!")
  }
  if (simulator == "multinomial") {
    BPPARAM = SerialParam()
  }

  method = match.arg(method)
  # normalize data
  if (is.vector(X)) {
    X <- as.matrix(X, ncol = 1)
  }

  # normalize genes
  if (is.null(genes.use)) {
    genes.use <- rownames(X)
  }else {
    genes.use <- genes.use
  }
  genes.use <- intersect(genes.use,
                         intersect(rownames(X),
                                   rownames(catalogue.profile)))
  X <- X[genes.use, ]
  catalogue.profile <- catalogue.profile[genes.use, ]

  catalogue.celltype <- factor(catalogue.celltype, levels = unique(catalogue.celltype))
  names(catalogue.celltype) <- colnames(catalogue.profile)


  n = nlevels(catalogue.celltype)
  if(length(sizes) < n) {
    s = floor(sizes / n)
    sizes = c(rep(s, n-1), sizes - (n-1)*s)
  }

  X = X[, sample(colnames(X), sum(sizes), replace = T)]
  Z = sapply(1:n, function(i){
    x = X[, (sum(sizes[1:i])-sizes[i]+1) : sum(sizes[1:i])]
    ct = levels(catalogue.celltype)[i]
    y = catalogue.profile[, catalogue.celltype == ct, drop = F]
    z = sizes[i]
    if (z>0) {
      Y = bplapply(1:z, function(j, x, y, method, lambda, epsilon, maxit, simulator){
        w = magrittr::mod(j, ncol(x))
        if (w==0) {
          w = ncol(x)
        }
        k = magrittr::mod(j, ncol(y))
        if (k==0) {
          k = ncol(y)
        }
        if (simulator == "scsimulator") {
          scSimulator::scSimulator(x[, w], y[, k], method, lambda, epsilon, maxit)
        }else{
          rmultinom(1, sum(x[, w]), y[, k])[, 1]
        }
      }, x, y, method, lambda, epsilon, maxit, simulator,
      BPPARAM = BPPARAM)
      Y = do.call("cbind", Y)
      rownames(Y) = rownames(catalogue.profile)
      return(Y)
    }
  }, simplify = F)
  Z = do.call("cbind", Z)
  colnames(Z) = paste("Sim", 1:sum(sizes), sep = "")

  res = list(scdata = Z, celltype=rep(levels(catalogue.celltype), sizes))
  return(res)
}
