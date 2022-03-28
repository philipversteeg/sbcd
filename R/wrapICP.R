# Copyright (c) 2021-2022, Philip Versteeg (p.j.j.p.versteeg@gmail.com). All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
suppressMessages(library(InvariantCausalPrediction))

simpleicp <- function(
  data, 
  contextVars, 
  alpha, 
  verbose=FALSE){
  
  # parse data
  p <- ncol(data) - length(contextVars)
  pAll <- ncol(data) 
  X <- as.matrix(data[,1:p]) # only system data here
  result <- matrix(0., nrow=ncol(data), ncol=ncol(data))

  # compile experiment indices from context, merge if multiple
  if(length(contextVars) == 1) {
    ExpInd <- data[,contextVars]
  } else {
    ExpInd <- rowSums(as.matrix(data[,contextVars]) %*% c(0:(metadata$pContext - 1))) 
  }

  # compute ICP for each target variable
  for (target in 1:p){
    Ytarget <- X[, target]
    Xtarget <- X[,-target]
    icp <- ICP(Xtarget, Ytarget, ExpInd=ExpInd, alpha=alpha, showCompletion=(verbose == 1), showAcceptedSets=FALSE, stopIfEmpty=TRUE)
    if(! icp$modelReject){
      if (any(icp$pvalues <= alpha)) {
        result[c(1:p)[-target][which(icp$pvalues<=alpha)], target] <- -log(icp$pvalues[which(icp$pvalues<=alpha)])
      }
    } 
  }
  return(result)
}