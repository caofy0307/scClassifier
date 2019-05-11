#' @title compute roc curve 
#' 
#' @description 
#' \code{roc2} computes roc curve coordinates for the given prediction result
#' 
#' @details 
#' 
#' @param labels known labels
#' @param predicts predicted lables
#' @param scores prediction score
#' @param step step size
#' 
#' @return 
#' 
#' @export
#' 
#' @import pROC
#' 
#' @author Feng Zeng <zengfeng@mail.xmu.edu.cn>
#' 
#' @examples 
#' 
roc2 <- function(labels, predicts, scores, step = 0.1, smooth = TRUE){
  
  if (!is.factor(labels)){
    labels <- as.factor(labels)  
  }
  
  if (!is.factor(predicts)){
    predicts <- as.factor(predicts)
  }
  
  sen <- NULL
  spe <- NULL
  auc <- NULL
  cil <- NULL
  cih <- NULL
  s <- seq(0,1,step)
  for (y in levels(predicts)) {
    x <- factor(labels %in% y, levels = c(TRUE, FALSE))
    if (all(x==FALSE)){
      aa <- rep(0, length(step))
      bb <- s
      cc <- 0
      dd <- 0
      ee <- 0
    } else if (all(x==TRUE)){
        aa <- rep(1, length(step))
        bb <- s
        cc <- 1
        dd <- 1
        ee <- 1
    }else {
      r <- pROC::roc(x, scores[, y], smooth = smooth)
      R <- pROC::coords(r, s, "specificity")
      
      aa <- R["sensitivity", ]
      bb <- R["specificity", ]
      ci <- as.numeric(pROC::ci(r))
      cc <- ci[1]
      dd <- ci[2]
      ee <- ci[3]
    }
    
    if (is.null(sen)){
      sen <- aa
      spe <- bb
      cil <- cc
      auc <- dd
      cih <- ee
    }else{
      sen <- cbind(sen, aa)
      spe <- cbind(spe, bb)
      cil <- c(cil, cc)
      auc <- c(auc, dd)
      cih <- c(cih, ee)
    }
  }
  
  sen <- as.data.frame(sen)
  colnames(sen) <- 1:ncol(sen)
  
  spe <- as.data.frame(spe)
  colnames(spe) <- 1:ncol(spe)
  
  auc <- as.data.frame(auc)
  colnames(auc) <- 1:ncol(auc)
  
  cil <- as.data.frame(cil)
  colnames(cil) <- 1:ncol(cil)
  
  cih <- as.data.frame(cih)
  colnames(cih) <- 1:ncol(cih)
  
  return(list(sensitivity = sen, specificity = spe, auc = auc, cil = cil, cih = cih))
  
}