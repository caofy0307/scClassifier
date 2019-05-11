#' @title XGBoost classifier
#'
#' @description
#' \code{XGBoostClassifier}
#'
#' @details
#'
#' @param X
#' @param catalogue.profile
#' @param catalogue.celltype
#' @param catalogue.prior
#' @param genes.use
#' @param cells.use
#' @param objective
#' @param probability
#' @param booster
#' @param eta
#' @param nthread
#' @param max.depth
#' @param nrounds
#' @param do.scale
#' @param importance.thresh
#' @param importance.quantile
#' @param BPPARAM
#'
#' @return
#'
#' @export
#'
#' @import parallel
#' @import xgboost
#' @import ranger
#' @import scran
#' @import Seurat
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
XGBoostClassifier <- function(X,
                              catalogue.profile,
                              catalogue.celltype,
                              catalogue.prior = NULL,
                              X.cluster = NULL,
                              genes.use = NULL,
                              cells.use = NULL,
                              objective = c("multi:softmax", "multi:softprob", "binary:logistic"),
                              probability = FALSE,
                              booster = c("gbtree", "gblinear"),
                              eta = 1,
                              nthread = detectCores(),
                              max.depth = 2,
                              nrounds = 1000,
                              do.scale = F,
                              importance.thresh = NULL,
                              importance.quantile = 0.8,
                              BPPARAM = SerialParam())
{

  objective = match.arg(objective)
  booster = match.arg(booster)
  # normalize genes
  if (is.null(genes.use)) {
    genes.use <- rownames(X)
  }else {
    genes.use <- genes.use
  }
  genes.use <- intersect(genes.use,
                         intersect(rownames(X),
                                   rownames(catalogue.profile)))

  # normalize cells
  if (is.null(cells.use)) {
    cells.use <- colnames(X)
  }else {
    cells.use <- cells.use
  }

  # normalize data
  if (is.vector(X)){
    X <- as.matrix(X, ncol=1)
  }
  X <- X[genes.use, cells.use]

  if (!is.null(X.cluster)) {
    X.cluster <- X.cluster[cells.use]
    X <- tapply(1:ncol(X), X.cluster, function(i) rowMeans(X[,i]))
    X <- do.call("cbind", X)
  }

  # normalize catalogue
  catalogue.profile <- catalogue.profile[genes.use, ]

  # batch correction
  sc <- CreateSeuratObject(cbind(X, catalogue.profile))
  sc@meta.data$type <- c(colnames(X), as.character(catalogue.celltype))
  sc@meta.data$batch <- c(rep(1, ncol(X)), rep(2, ncol(catalogue.profile)))
  sc <- NormalizeData(sc)

  debatch <- fastMNN(as.matrix(sc@data[,sc@meta.data$batch == 1]),
                     as.matrix(sc@data[,sc@meta.data$batch == 2]),
                     BPPARAM = BPPARAM)
  colnames(debatch$corrected) <- paste("Var", 1:ncol(debatch$corrected), sep = "")
  debatch.X <- debatch$corrected[sc@meta.data$batch == 1,]
  rownames(debatch.X) <- sc@cell.names[sc@meta.data$batch == 1]
  debatch.catalogue <- debatch$corrected[sc@meta.data$batch == 2,]
  rownames(debatch.catalogue) <- sc@cell.names[sc@meta.data$batch == 2]

  if (do.scale){
    debatch.X <- scale(debatch.X)
    debatch.catalogue <- scale(debatch.catalogue)
  }

  # use random forest to select informative features
  vars.use <- colnames(debatch.catalogue)
  if (is.null(importance.thresh)) {
    vars.use <- RandomForestGenes(t(debatch.catalogue), factor(sc@meta.data[sc@meta.data$batch==2, "type"]))
  }else{
    rf.train <- as.data.frame(debatch.catalogue)
    rf.train$y <- sc@meta.data[sc@meta.data$batch==2, "type"]
    rf.model <- ranger(y~., data = rf.train, importance = "impurity")
    importance.thresh <- min(importance.thresh, quantile(rf.model$variable.importance, importance.quantile))
    vars.use <- vars.use[rf.model$variable.importance >= importance.thresh]
  }
  debatch.catalogue <- debatch.catalogue[, vars.use]
  debatch.X <- debatch.X[, vars.use]

  # training xgboost
  if (probability) {
    objective = "multi:softprob"
  }

  catalogue.label <- as.numeric(factor(catalogue.celltype)) - 1
  catalogue.num_label <- length(unique(catalogue.label))
  classifier <- xgb.train(data = xgb.DMatrix(data = debatch.catalogue, label = catalogue.label),
                          num_class = catalogue.num_label,
                          booster = booster, eta = eta, max.depth = max.depth,
                          nthread = nthread, nrounds = nrounds, objective = objective)

  # predict
  y <- predict(classifier, debatch.X)
  y <- plyr::mapvalues(y,
                       from  = levels(factor(catalogue.label)),
                       to = levels(factor(catalogue.celltype)))

  # results
  if (is.null(X.cluster)) {
    if (probability) {
      results <- list(celltype = factor(apply(y, 1, function(x) names(x)[which.max(x)])),
                      score = apply(y, 1, max),
                      score.matrix = y)
    }else{
      results <- list(celltype = factor(y),
                      score = NULL,
                      score.matrix = NULL)
    }
  }else{
    if (probability) {
      C <- apply(y, 1, function(x) names(x)[which.max(x)])
      S <- apply(y, 1, max)
      results <- list(celltype = factor(catalogue.celltype[C[X.cluster]]),
                      score = S[X.cluster],
                      score.matrix = y[, X.cluster])
    }else{
      results <- list(celltype = factor(y[X.cluster]),
                      score = NULL,
                      score.matrix = NULL)
    }
  }

  results$model <- classifier
  results$train <- list(x = debatch.catalogue, y = catalogue.celltype)
  results$test <- list(x = debatch.X)

  return(results)
}
