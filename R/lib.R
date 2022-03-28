# Copyright (c) 2021-2022, Philip Versteeg (p.j.j.p.versteeg@gmail.com). All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
suppressMessages(library(RBGL))
suppressMessages(library(expm))

checkGraphConnected <- function(adjacency) {
  return(all((rowSums(adjacency != 0) > 0) | (colSums(adjacency != 0) > 0))) # this might be faster?
}

graphColliders <- function(adjacency) {
  tmp <- colSums(adjacency != 0) >= 2
  if (any(tmp)) {
    return(which(tmp))
  } else {
    return(c())
  }
}

randomGraph <- function(p, eps, pSel=1, acyclic=TRUE, requireNumColliders=0, requireSimplyConnected=FALSE, seed=1) {
  max.iters <- 1E5
  pSysObs <- p
  pSysConf <- 0
  pContext <- 1

  set.seed(seed)

  system.time(
    repeat {
      ATmp <- matrix(0,pSysObs + pSysConf,pSysObs + pSysConf)
      V <- matrix(as.numeric(runif(pSysObs*pSysObs)<eps),pSysObs,pSysObs)
      ATmp[1:pSysObs,1:pSysObs] <- V 

      # remove (self)loops or keep
      if (acyclic) {
        ATmp <- ATmp - lower.tri(ATmp,diag=TRUE) * ATmp
      } else {
        doneCylicity <- FALSE  
        while( !doneCylicity ) {
          ATmp <- ATmp - diag(diag(ATmp))
          doneCylicity <- sum((diag(expm(ATmp[1:pSysObs,1:pSysObs]))-rep(1,pSysObs))^2) < 1e-10
        }
      }

      # check for potential context variables
      potConVars <- intersect(which(colSums(ATmp) == 0), which(rowSums(ATmp) != 0))
      if (length(potConVars) != 0) {
        # check for certain conditions, then terminiate
        if (requireNumColliders > 0 & requireSimplyConnected) {
          if (length(graphColliders(ATmp)) > requireNumColliders & checkGraphConnected(ATmp)) {
            break
          } 
        } 
        if (requireNumColliders > 0 & !requireSimplyConnected) {
          if (length(graphColliders(ATmp)) > requireNumColliders) {
            break
          }
        }
        if (!requireNumColliders > 0 & requireSimplyConnected) {
          if (checkGraphConnected(ATmp)) {
            break
          }
        }
        if (!requireNumColliders > 0 & !requireSimplyConnected){
          break
        }
      }

      if (max.iters == 0) {
        cat('Max number of samples reached!')
        stop('Max number of samples reached!')
      }
      max.iters <- max.iters - 1 
    }
  ) -> time.loop
  cat('Sampled graph in', 1E5 - max.iters, 'iterations.\n')
  
  conVars <- sample(c(potConVars), 1)
  # and rotate ATmp to make it the last variable
  ATmp[,c(conVars,p)] <- ATmp[,c(p,conVars)]
  ATmp[c(conVars,p),] <- ATmp[c(p,conVars),]

  # sample selection bias variables 
  coliderVars <- graphColliders(ATmp) 
  anc<-sign(expm(ATmp))
  coliderVarsDesc <- which(colSums(anc[coliderVars,,drop=FALSE]) > 0)
  sysSel <- sample(c(coliderVarsDesc), 1) 
  if (pSel > 1) {
    sysSel <- c(sysSel, sample(c(1:pSysObs)[-which(c(1:pSysObs) == sysSel)],pSel - 1, replace=FALSE))
  }
  return(list(A=ATmp, sysSel=sysSel))
}

extractPosNeg <- function(pred, true, conVars) {
  stopifnot(all(dim(pred) == dim(true))) 
  tr <- true
  diag(tr) <- NA
  tr[,conVars] <- NA
  tr[conVars,] <- NA
  ind.pos <- which(tr == 1)
  ind.neg <- which(tr == 0)
  ret <- list(
    pos = pred[ind.pos],
    neg = pred[ind.neg]
  )
  return(ret)
}

plotPRCurveMarker <- function(pr.result, col, value, lwd, pch) {
  index.score <- which(pr.result$curve[,3] >= value)[1]
  if (!(index.score %in% c(dim(pr.result$curve)[1], (dim(pr.result$curve)[1]-1), (dim(pr.result$curve)[1]-2)))) { # check that the marker is not at the end
    marker.xy <- pr.result$curve[index.score,c(1,2)]
    points(x=marker.xy[1], y=marker.xy[2], pch=pch, col=col, lwd=lwd) #, cex=2)
  }
}

createSuffStatOracle <- function(A, sysSel, bias) {
  B <- matrix(0,nrow=nrow(A)+1, ncol=nrow(A)+1)
  B[1:p, 1:p] <- A
  B[sysSel,p+1] <- which(sysSel != 0)
  B <- as(B,'graphNEL')
  S <- NULL
  if (bias) S <- p+1
  return(list(g=B, jp=johnson.all.pairs.sp(B), S=S))
}