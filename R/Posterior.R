#' @title Posterior probability
#'
#' @description
#' \code{Posterior}
#'
#' @details
#'
#' @param likeli
#' @param prior
#' @param use.log
#' @param BPPARAM
#'
#' @return
#'
#' @export
#'
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'

Posterior <- function(likeli, prior = NULL, use.log = FALSE){

  # joint probability
  if (is.null(prior)){
    if (use.log) joint <- exp(likeli)
    else joint <- likeli
  }else{
    priot <- log(prior)
    if (!use.log) likeli <- log(likeli)

    joint <- likeli + prior
    joint <- exp(joint)
  }

  # posterior
  posterior <- joint/sum(joint)

  return(posterior)
}
