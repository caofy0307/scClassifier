
VotingClassifier <- function(X, model, celltypes = NULL, method = c("multinomial", "naivebayes", "dmc", "dmc2"),
                             genes.use = NULL, cells.use = NULL,
                             cells.cluster = NULL, ...){

  method <- match.arg(method)

  voting <- function(x){
    # cell-wise classification
    if (method == "multinomial"){
      Z <- multinomial_classifier(x, prob = model, celltypes = celltypes, genes.use = genes.use, cells.use = cells.use, ...)
    }else if (method == "dmc"){
      z <- dmc(x, model, celltypes, ...)
      k <- apply(z, 1, function(x) names(x)[which.max(x)])
      Z <- list(celltypes = levels(celltypes),
                population = k)
    }else if (method == "dmc2"){
      z <- dmc2(x, model, celltypes, ...)
      k <- apply(z, 1, function(x) names(x)[which.max(x)]) %>% as.factor
      Z <- list(celltypes = levels(celltypes),
                population = k)
    }

    # voting
    V <- rep(0, length(Z$celltypes))
    names(V) <- Z$celltypes

    U <- table(as.factor(Z$population))
    V[names(U)] <- U

    V <- sort(V, decreasing = TRUE)

    # results
    structure(list(population = names(V[1]), voting = V), class = "voter")
  }

  if (is.null(cells.cluster)) {
    results <- voting(X)
  } else {
    results <- tapply(cells.cluster, cells.cluster,
                      function(z) voting(X[, names(z)])
                      )
  }

  return(results)
}
