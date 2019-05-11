#' @title Fit Dirichlet distribution to cell type profile
#' 
#' @description 
#' \code{dirichlet_fit} estimates the parameters of cell type-associated Dirichlet distribution
#' 
#' @details 
#' 
#' @param target.profiles cell type-associated profile matrix
#' @param target.identity a factor of cell type identity
#' 
#' @return 
#' 
#' @export
#' 
#' @import MCMCprecision
#' 
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#' 
#' @examples 
#' 
dirichlet_fit <- function(target.profile, target.identity){
  if (length(target.identity) == 1){
    x <- target.profile / colSums(target.profile)
    x <- t(x)
    # y <- dirichlet.mle(x)
    y <- fit_dirichlet(x)
    y <- y$alpha
    colnames(y) <- target.identity
  }else{
    if (!is.factor(target.identity)){
      target.identity <- as.factor(target.identity)
    }
    
    y <- tapply(colnames(target.profile), target.identity, function(i){
      if (length(i) == 1) {
        r <- list(alpha = target.profile[, i])
      }else{
        x <- target.profile[, i] / colSums(target.profile[, i, drop = FALSE])
        x <- t(x)
        r <- fit_dirichlet(x)
      }
      return(r$alpha)
    })
    
    y <- matrix(unlist(y), ncol = nlevels(target.identity), byrow = FALSE)
    colnames(y) <- levels(target.identity)
  }
  
  rownames(y) <- rownames(target.profile)
  
  return(y)
}