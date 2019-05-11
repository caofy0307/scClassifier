#' @title Receiver operating characteristics of a single cell classifier
#'
#' @description
#' \code{roc} calculates the receiver operating characteristics (ROC) for a single cell classifier.
#'
#' @details
#'
#' @param x a matrix of data.frame of cell data
#' @param y a vector or data.frame of cell identity
#' @param method single cell classifier
#' @param genes.use a vector of genes to be used
#' @param pathways.use a list of pathways to be used
#' @param nfolds the number of cross-validation folds
#' @param foldid a vector of indicators assigning cells to folds
#' @param plot plot ROC curve
#' @param errbar.color error bar color
#' @param errbar.width error bar width
#' @param point.size point size
#' @param point.shape point shape
#' @param point.fill point fill color
#' @param line.color line color
#' @param line.width line width
#' @param ... further arguments
#'
#' @return a list of ROC results and a ggplot2 object
#'
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#'
#' @author Feng Zeng <zengfeng@xmu.edu.cn>
#'
#' @examples
#'
#'
roc <- function(x, y, method = c("multinomial", "naivebayes"),
                genes.use = NULL, pathways.use = NULL,
                nfolds = 5, foldid = NULL,
                plot.roc = FALSE, errbar.color = "black", errbar.width = 0.015,
                point.size = 2, point.shape = 21, point.fill = "pink",
                line.color = "grey85", line.width = 1, ...){
  method <- match.arg(method)
  cells.n <- ncol(x)

  if (!is.null(genes.use)) x <- x[genes.use, ]
  x <- as.matrix(x)

  if (is.null(foldid)) {
    fold.size <- ceiling(cells.n / nfolds)
    foldid <- rep(1:nfolds, fold.size)[1:cells.n]
    foldid <- sample(foldid)
  }

  foldid <- as.factor(foldid)
  nfolds <- nlevels(foldid)

  fpr.seq <- seq(0,1,0.01)
  # label should be 0 (negative) or 1 (positive)
  # score should be vector
  tprfpr <- function(label, score){
    score <- round(score, 4)
    thresh <- sort(score, decreasing = TRUE) %>% unique
    thresh <- c(1, thresh, 0)
    #if (thresh[length(thresh)] >0 ) thresh <- c(thresh, 0)

    count <- function(i) {
      cutoff <- thresh[i]
      pred.label <- ifelse(score >= cutoff, 1, 0)
      freq <- table(pred.label, label)
      if ("1" %in% rownames(freq)) {
        tpr <- freq["1","1"]/sum(freq[,"1"])
        fpr <- freq["1","0"]/sum(freq[,"0"])
      }else {
        tpr <- 0.0
        fpr <- 0.0
      }

      fpr <- fpr + i*1e-6
      return(c(fpr=fpr, tpr=tpr, cutoff=cutoff))
    }

    ft <- sapply(1:length(thresh), count)
    return(approx(ft["fpr",], ft["tpr",], fpr.seq, rule = 2)$y)
  }

  # macro-average method to generate tpr-fpr for multi-class
  macroave <- function(label, score, level){
    foo <- function(l){
      real.y <- ifelse(label==l, 1, 0)
      pred.y <- unlist(score[l,])
      tprfpr(real.y, pred.y)
    }

    ave <- sapply(level, foo)
    ave <- apply(ave,1,mean)
    return(c(0,ave))
  }

  foldval <- function(id){
    ind <- which(foldid==id)

    # split data into training and test subsets
    train.x <- x[,-ind]
    if (any(is.data.frame(y), is.array(y))) {
      train.y <- y[-ind, ]
    }else {
      train.y <- y[-ind]
    }

    test.x <- x[,ind]
    if (any(is.data.frame(y), is.array(y))) {
      test.y <- y[ind, ]
    }else {
      test.y <- y[ind]
    }

    if (method == "multinomial"){
      # train model
      prob <- apply(train.x, 2, function(z) z/sum(z))
      celltypes <- train.y
      pred.y <- multinomial_classifier(round(test.x), prob = prob, celltypes = celltypes)

      # label
      label <- test.y$population
      score <- pred.y$post
      level <- levels(as.factor(test.y$population))

      return(macroave(label, score, level))
    }
    # other classification methods
  }

  results <- sapply(1:nfolds, foldval)
  results <- as.data.frame(results)
  colnames(results) <- paste0("tpr", 1:nfolds, seq="")
  results$fpr <- c(0,fpr.seq)

  g <- NULL
  if (plot.roc){
    z <- results
    z$ind <- rownames(z) %>% as.numeric
    z <- melt(z, id=c("ind","fpr"))
    colnames(z) <- c("ind", "fpr", "fold", "tpr")
    z <- summarySE(z, measurevar = "tpr", groupvars = c("ind", "fpr"))

    g <- ggplot(z, aes(x = fpr,y = tpr)) +
      geom_errorbar(aes(ymin = tpr-se,ymax = tpr+se), width = errbar.width, color = errbar.color) +
      geom_line(size = line.width, color = line.color) +
      geom_point(size = point.size, shape = point.shape, fill = point.fill) +
      xlab("1-specificity") +
      ylab("sensitivity")
    g
  }

  structure(list(roc = results,
                 gplot = g),
            class = "roc")
}

#' @title Summarizes data
#'
#' @description
#' \code{summarySE} gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
#'
#' @param data a data frame.
#' @param measurevar the name of a column that contains the variable to be summariezed
#' @param groupvars a vector containing names of columns that contain grouping variables
#' @param na.rm a boolean that indicates whether to ignore NA's
#' @param conf.interval the percent range of the confidence interval (default is 95%)
#'
#' @return
#'
#' @export
#'
#' @import plyr
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
