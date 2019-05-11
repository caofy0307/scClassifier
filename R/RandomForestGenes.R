#' @title RandomForestGenes
#'
#' @description
#' \code{RandomForestGenes}
#'
#' @details
#'
#' @param X
#' @param y
#' @param size
#' @param do.par
#'
#' @return
#'
#' @export
#'
#' @import caret
#' @import doMC
#' @import parallel
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
RandomForestGenes <- function(X, y, size = 500, do.par = TRUE){

  registerDoMC(cores = detectCores()-1)

  control <- rfeControl(functions=rfFuncs, method="cv", number=10, allowParallel = do.par)
  results <- rfe(t(X), y, sizes = size, rfeControl = control)
  genes.use <- predictors(results)
  return(genes.use[1:min(size, length(genes.use))])
}
