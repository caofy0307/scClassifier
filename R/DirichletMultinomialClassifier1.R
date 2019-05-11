#' @title Dirichlet multinomial classifier for one data
#'
#' @description
#' \code{DirichletMultinomialClassifier1} computes the predictive probability for the given cells following
#' Multinomial distribution model with the Dirichlet conjugate prior.
#'
#' @details
#'
#' @param x
#' @param param
#' @param use.log
#'
#' @return
#'
#' @export
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
DirichletMultinomialClassifier1 <- function(x,
                                            param,
                                            use.log = FALSE) {
  a <- lgamma(x+param)
  b <- lgamma(sum(x+param))
  c <- lgamma(sum(param))
  d <- lgamma(param)
  e <- sum(a)-b+c-sum(d)
  if (!use.log) {
    e <- exp(e)
  }
  return(e)
}
