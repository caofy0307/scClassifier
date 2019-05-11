#' @title Multinomial model likelihood
#'
#' @description
#' \code{MultiLikelihood}
#'
#' @details
#'
#' @param X
#' @param p
#' @param use.log
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

MultiLikelihood <- function(x, p, use.log = FALSE, BPPARAM = SerialParam()){

  p[p<1e-6] <- 1e-6

  if (is.vector(p)){
    p <- p/sum(p)
    likeli <- dmultinom(x, prob = p, log = use.log)
  }else{
    likeli <- bplapply(1:ncol(p), function(i){
      q <- p[, i]
      q <- q/sum(q)
      dmultinom(x, prob = q, log = use.log)
    }, BPPARAM = BPPARAM) %>% do.call("cbind", .) %>% as.matrix
  }
  return(likeli)
}
