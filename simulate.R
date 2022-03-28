# Copyright (c) 2021-2022, Philip Versteeg (p.j.j.p.versteeg@gmail.com). All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
suppressMessages(library(PRROC))

# helper functions 
# methods
source('R/lib.R')
source('R/wrapYStructures.R') # std ystruc, over all vars
source('R/wrapLCD.R')
source('R/wrapICP.R')

compareTypes <- c(
  # 'anc', 
  'pattern',  
  NULL
)

# unique name for fig
figBase <- '' 
removePlotTail <- TRUE

# defaults
N <- 10000
nseed <- 5
verbose.method <- 0
experiments <- c('fixed', 'small', 'large', 'small_alt', 'large_alt')
exp <- 3

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  if (suppressWarnings(is.na(as.integer(args[1])))) {
    exp <- which(experiments == args[1])
  } else {
    exp <- as.integer(args[1])
  }
} 
if (length(args) > 1) {
  nseed <- as.integer(args[2])
} 
if (length(args) > 2) {
  N <- as.integer(args[3])
} 
if (length(args) > 3) { 
  if (args[4] %in% c('p','pat', 'pattern', 'patterns', 'oracle')) {
    compareTypes <- c('pattern', NULL)
  } else if (args[4] %in% c('a', 'anc', 'ancestors', 'an')) {
    compareTypes <- c('anc', NULL)
  } else if (args[4] %in% c('s', 'stat', 'stats')) {
    compareTypes <- c(NULL)
  } else {
    stop('No valid compareType!')
  }
} 
expName <- experiments[exp]

cat('Call script with args: Rscript selectionTest.R [#exp] [#seed] [#samples] [ancestors/patterns]\n')
cat('*** EXPERIMENT:', expName, '***\n')
cat('*** NSEED:', nseed, '***\n')
cat('*** N:', N, '***\n')
cat('*** PREDICTION TARGET:', compareTypes[1], '***\n')


methods <- c(
  'yst.full.gsn',
  'ystext.full.gsn',
  'lcd.gsn',
  'icp', 
  NULL
)
methodsPrettyName <- list(
  'yst.full.gsn'='YSt',
  'ystext.full.gsn'='YSt-Ext',
  'icp'='ICP',
  'lcd.gsn'='LCD',
  NULL
)

if ('pattern' %in% compareTypes) {
  # add all oracle methods
  oracle.methods <- NULL
  other.methods <- NULL
  for (i in methods) {
    if (grepl('.oracle', i)) {
      method.oracle <- i 
    } else {
      other.methods <- c(other.methods, i)
    }
    if (grepl('.gsn', i)) method.oracle <- gsub('.gsn', '.oracle', i)
    if (!method.oracle %in% oracle.methods) {
      oracle.methods <- c(oracle.methods, method.oracle)
    }
  }
  methods <- c(oracle.methods, other.methods)
}

# keep stats for each seed
stats <- list()
for (i in methods) {
  stats[[i]] <- list()
  stats[[i]]$biased <- list()
  stats[[i]]$unbiased <- list()
  for (j in names(stats[[i]])) {
    stats[[i]][[j]] <- list()
    stats[[i]][[j]]$tp <- 0
    stats[[i]][[j]]$fp <- 0
    stats[[i]][[j]]$posscores <- NULL
    stats[[i]][[j]]$negscores <- NULL
    stats[[i]][[j]]$posscores.pattern <- NULL
    stats[[i]][[j]]$negscores.pattern <- NULL
  }
}

for (seed in c(1:nseed)) {
  cat('\n*** seed', seed, '***\n')
  set.seed(seed)
  if (expName == 'fixed') {
    p <- 7
    sysSel <- c(1,2)
    A <- matrix(c(
      0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,0,0,1,0,0,
      0,0,0,0,0,0,0,
      0,0,0,1,0,0,0,
      1,0,0,1,0,0,0 
    ), nrow=p, ncol=p, byrow=TRUE)
  } else if (expName == 'small') {
    p <- 9
    pSel <- 1
    eps <- 0.15
    requireNumColliders <- 3
    requireSimplyConnected <- FALSE
    sampledGraph <- randomGraph(p=p,pSel=pSel,eps=eps,requireNumColliders=requireNumColliders,requireSimplyConnected=requireSimplyConnected,acyclic=TRUE,seed=seed)
    A <- sampledGraph$A
    sysSel <- sampledGraph$sysSel
  } else if (expName == 'large') {
    p <- 17
    pSel <- 3
    eps <- 0.09
    requireNumColliders <- 5
    requireSimplyConnected <- FALSE
    sampledGraph <- randomGraph(p=p,pSel=pSel,eps=eps,requireNumColliders=requireNumColliders,requireSimplyConnected=requireSimplyConnected,acyclic=TRUE,seed=seed)
    A <- sampledGraph$A
    sysSel <- sampledGraph$sysSel
  } else if (expName == 'small_alt') {
    p <- 9
    pSel <- 1
    eps <- 0.12
    requireNumColliders <- 3
    requireSimplyConnected <- TRUE
    sampledGraph <- randomGraph(p=p,pSel=pSel,eps=eps,requireNumColliders=requireNumColliders,requireSimplyConnected=requireSimplyConnected,acyclic=TRUE,seed=seed)
    A <- sampledGraph$A
    sysSel <- sampledGraph$sysSel
  } else if (expName == 'large_alt') {
    p <- 17
    pSel <- 3
    eps <- 0.11
    requireNumColliders <- 5
    requireSimplyConnected <- FALSE
    sampledGraph <- randomGraph(p=p,pSel=pSel,eps=eps,requireNumColliders=requireNumColliders,requireSimplyConnected=requireSimplyConnected,acyclic=TRUE,seed=seed)
    A <- sampledGraph$A
    sysSel <- sampledGraph$sysSel
  } else {
    cat('\n********************\nNo existing experiment!\n********************\n')
    break
  }
  
  pCon <- 1
  pSys <- p - pCon
  sysVars <- c(1:(p-1))
  conVars <- p

  # ancestor matrix without convars
  anc<-sign(expm(A[1:(pSys+pCon),1:(pSys+pCon)]))
  anc<-anc-diag(diag(anc))
  anc[p,] <- 0
  anc[,p] <- 0

  # sample and rescale weights
  weights <- sign(matrix(rnorm(p*p),p,p)) * (matrix(runif(p*p),p,p) + 0.5)
  B <- A * weights
  scales<-matrix(0,p,1)
  for( i in (1):p ) {
    scales[i] <- (sum(B[sysVars,i]^2) + 1)
    B[,i] <- B[,i] / scales[i]^0.5
  }
  # unbiased data
  E <- matrix(data=rnorm(N*p),nrow=N,ncol=p)
  XUnbiased <- t(solve(diag(p)-t(B), t(E)))

  # biased
  max_iters <- 1
  XBiased <- c()
  while(max_iters <= 1000) {
    E <- matrix(data=rnorm(100*N*p),nrow=100*N,ncol=p)
    Xtmp <- t(solve(diag(p)-t(B), t(E)))
    sbVar <- apply(Xtmp[,sysSel, drop=FALSE], 1, sum)
    selectThese <- which(sbVar > 2 & sbVar < 2.5)

    # create or add to XBiased
    if (is.null(XBiased)) {
      XBiased <- Xtmp[selectThese,]
    } else {
      XBiased <- rbind(XBiased, Xtmp[selectThese,])
    }
    if (dim(XBiased)[1] >= N) {
      XBiased <- XBiased[1:N,] # at most N sampled
      break
    }
    max_iters <- max_iters + 1
  }

  ancIndices <- which(anc > 0, arr.ind=TRUE)
  ancIndices <- ancIndices[order(ancIndices[,2]),]
  ancIndices <- ancIndices[order(ancIndices[,1]),]
  totalNonContextEdges <- 0
  cat('True ancestral relations:\n')
  for (i in 1:nrow(ancIndices)) {
    ca <- ancIndices[i,1]
    ef <- ancIndices[i,2]
    if (ca != ef) {
      cat(ca,'->',ef, '\n')
      if (!(ca %in% conVars | ef %in% conVars)) {
        totalNonContextEdges <- totalNonContextEdges + 1
      }
    }
  }
  cat('\n')

  # run methods for biased and unbiased
  for (bias in c(FALSE,TRUE)) {
    if (bias)
      data <- XBiased
    else
      data <- XUnbiased
    b <- ifelse(bias, 'biased', 'unbiased')
    oracle.patterns <- list()
    for (method in methods) {
      cat('- ', method, b, '\n')

      # Y-Structures
      if (startsWith(method, 'yst')) {
        # test
        if (grepl('.gsn', method)) yst.test <- 'gaussCItest'
        if (grepl('.oracle', method)) {
          yst.test <- 'oracle'
          suffStatOracle <- createSuffStatOracle(A=A, sysSel=sysSel, bias=(bias=='bias'))
        } else {
          suffStatOracle <- NULL
        }
        # full or extended Y struc
        yst.fullYStruc <- TRUE
        if (startsWith(method, 'ystext')) yst.fullYStruc <- FALSE
        suppressWarnings(
          resultlist <- ystructures(
            data=data, 
            systemVars=sysVars,
            contextVars=conVars,
            alpha=1e-2,
            beta=NULL,
            test=yst.test,
            fullYStruc=yst.fullYStruc,
            verbose=verbose.method,
            suffStatOracle=suffStatOracle
          )
        )
        result <- resultlist$arel
        pat.size <- 4
      }

      # LCD
      if (startsWith(method, 'lcd.')) {
        if (grepl('.gsn', method)) lcd.test <- 'gaussCItest'
        if (grepl('.oracle', method)) {
          lcd.test <- 'oracle' # todo
          suffStatOracle <- createSuffStatOracle(A=A, sysSel=sysSel, bias=(bias=='bias'))
        } else {
          suffStatOracle <- NULL
        }
        suppressWarnings(
          resultlist <- lcd(
            data=data,
            alpha=1e-2,
            beta=NULL,
            test=lcd.test,
            systemVars=sysVars,
            contextVars=conVars,
            verbose=verbose.method,
            suffStatOracle=suffStatOracle
          )
        )
        result <- resultlist$arel
        pat.size <- 3
      }

      # ICP
      if (method == 'icp') {
        icpdata <- data 
        # threshold continous context to binary for ICP environment
        icpdata[,conVars] <- (icpdata[,conVars] > mean(icpdata[,conVars]))
        result <- simpleicp(data=icpdata,
          contextVars=conVars,
          alpha=1e-2,
          verbose=0
        )
      }

      # compile results for each oracle method
      if (grepl('.oracle', method)) {
        oracle.pat <- resultlist$patterns
        if (length(oracle.pat) > 0) {
          for (i in 1:length(oracle.pat)) {
            oracle.pat[[i]] <- oracle.pat[[i]][1:pat.size]
          }
        }
        oracle.patterns[[method]] <- oracle.pat
      }

      # compile stats
      result <- result-diag(diag(result))
      tp <- sum((result[-conVars,-conVars] != 0) & anc[-conVars,-conVars])
      fp <- sum((result[-conVars,-conVars] != 0) & !anc[-conVars,-conVars])
      stats[[method]][[b]]$tp <- stats[[method]][[b]]$tp + tp
      stats[[method]][[b]]$fp <- stats[[method]][[b]]$fp + fp
      extr <- extractPosNeg(result, anc, conVars)
      stats[[method]][[b]]$posscores <- c(stats[[method]][[b]]$posscores, extr$pos)
      stats[[method]][[b]]$negscores <- c(stats[[method]][[b]]$negscores, extr$neg)
      
      # if comparing to oracle pattern is required, we need to add the relevent method to the method list and run it.
      if ('pattern' %in% compareTypes) {
        if ((!grepl('.oracle', method)) & (startsWith(method, 'yst') | startsWith(method, 'lcd') | startsWith(method, 'weirdpat'))) {
          # extract oracle pat and throw away the oracle score
          if (grepl('.gsn', method)) method.oracle <- gsub('.gsn', '.oracle', method)
          if (grepl('.spr', method)) method.oracle <- gsub('.spr', '.oracle', method)

          # get score prediction for each pattern and compute posscores negscores
          oracle.pat <- oracle.patterns[[method.oracle]]
          pred.pat <- resultlist$patterns
          negscores.pattern.new <- c()
          posscores.pattern.new <- c()
          for (pp in pred.pat) {
            # check for each pattern if it is contained in oracle pattern list
            if (any(sapply(oracle.pat, function(x) all(x == pp[1:pat.size])))) {
              posscores.pattern.new <- c(posscores.pattern.new, pp[length(pp)])
            } else {
              negscores.pattern.new <- c(negscores.pattern.new, pp[length(pp)])
            }
          }
          stats[[method]][[b]]$negscores.pattern <- c(stats[[method]][[b]]$negscores.pattern, negscores.pattern.new)
          stats[[method]][[b]]$posscores.pattern <- c(stats[[method]][[b]]$posscores.pattern, posscores.pattern.new)
        }        
      }
      cat('  ', sum(extr$pos != 0) + sum(extr$neg != 0), 'edges found\n')
    }
  }
}

cat('\n*** STATS ***\n')
for (method in methods) {
  cat('-', method,'TP FP\n')
  cat('unbiased', stats[[method]]$unbias$tp, stats[[method]]$unbias$fp,'\n')
  cat('  biased', stats[[method]]$bias$tp, stats[[method]]$bias$fp, '\n')
}

cat('\n*** PLOTS ***\n')
expID <- paste0(figBase, '.', N, '.', nseed)

# setup folders
figFolder <- paste0('results')
dir.create(figFolder, showWarnings=FALSE)

# setup plot
colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728')
zoomXMax <- c(1, 0.145, 0.105)
zoomNames <- c('', '_view1', '_view2')
addMarkerToCurve <- TRUE
scoreValueForMarker <- -log(0.01) 
markerLwd <- 2.5
markerPch <- 4

if ('anc' %in% compareTypes) { 
  # test if each method has at least some positives in biased and unbiased
  removeThese <- c(NULL)
  for (method in methods) {
    if (sum(stats[[method]][['biased']]$posscores != 0) == 0 | sum(stats[[method]][['unbiased']]$posscores != 0) == 0) {
      removeThese <- c(removeThese,method)
    }
  }
  if (length(removeThese) > 0) {
    cat('** Remove these methods due to zero pos in either biased or unbiased **', '\n')
    cat('\t',removeThese, '\n')
    theseMethods <- setdiff(methods, removeThese)
    cat('   Remaining methods:\n\t', theseMethods, '\n')
  } else {
    theseMethods <- methods
  }

  # pretty name for legend
  methodsPretty <- sapply(theseMethods, function(x) methodsPrettyName[[x]])
  if (length(theseMethods)) {
    for (z in 1:length(zoomNames)) {
      pdf(paste0(figFolder, '/', expName, '_anc_', N, '_', nseed, zoomNames[z], '.pdf'))
      iter <- 0
      for (method in theseMethods) {
        iter <- iter + 1
        b <- 'biased'
        pr.bias <- pr.curve(scores.class0=stats[[method]][[b]]$posscores, scores.class1=stats[[method]][[b]]$negscores, rand.compute=TRUE, curve=TRUE)
        if (removePlotTail) pr.bias$curve <- pr.bias$curve[-which(pr.bias$curve[,3] == 0),]
        if (iter == 1)
          plot(pr.bias, main=paste0('PR'), lwd=2, xaxs="i", yaxs="i", color=colors[iter], pch=0, xlim=c(0,zoomXMax[z]), rand.plot=TRUE, auc.main=FALSE)
        else 
          plot(pr.bias, add=TRUE, main=paste0('PR'), lwd=2, xaxs="i", yaxs="i", color=colors[iter], pch=0, xlim=c(0,zoomXMax[z]), rand.plot=TRUE, auc.main=FALSE)
        if (addMarkerToCurve) plotPRCurveMarker(pr.bias, value=scoreValueForMarker, lwd=markerLwd, pch=markerPch, col=colors[iter])
        b <- 'unbiased'
        pr.unbias <- pr.curve(scores.class0=stats[[method]][[b]]$posscores, scores.class1=stats[[method]][[b]]$negscores, rand.compute=TRUE, curve=TRUE)
        if (removePlotTail) pr.unbias$curve <- pr.unbias$curve[-which(pr.unbias$curve[,3] == 0),]
        plot(pr.unbias, add=TRUE, main=paste0('PR'), lwd=2, xaxs="i", yaxs="i", color=colors[iter], lty='dashed')
      }
      legend("bottomright",legend=methodsPretty,col=colors[1:length(methodsPretty)],pch=15,bty="n",pt.cex=2,cex=1.2,text.col="black",horiz=F,inset=c(0.01, 0.01))
      dev.off()
    }
  }
}

# output: oracle pattern plot
if ('pattern' %in% compareTypes) {
  cat('*** PATTERN PLOTS ***\n')
  theseMethods <- NULL
  for (method in methods) {
    method.oracle <- NULL
    if (grepl('.gsn', method)) method.oracle <- gsub('.gsn', '.oracle', method) 
    if (grepl('.spr', method)) method.oracle <- gsub('.spr', '.oracle', method)
    if (!is.null(method.oracle)) {
      # check no empty pos scores in both pos and neg:
      cat(' - ', method, ' #pos patterns:', sum(stats[[method]][['biased']]$posscores.pattern != 0), 
      '; #neg patterns:', sum(stats[[method]][['biased']]$negscores.pattern != 0), '\n')
      if (!(sum(stats[[method]][['biased']]$posscores.pattern != 0) == 0 | sum(stats[[method]][['unbiased']]$posscores.pattern != 0) == 0)) {
        if (! method %in% theseMethods) {
          theseMethods <- c(theseMethods, method)
        }
      } 
    }
    # safeguard for if negscores is empty
    for (b in c('biased', 'unbiased')) {
      if (is.null(stats[[method]][[b]]$negscores.pattern)) {
        stats[[method]][[b]]$negscores.pattern <- numeric()
      }
    }
  }
  
  # pattern comparisons
  if (length(theseMethods)) {
    methodsPretty <- sapply(theseMethods, function(x) methodsPrettyName[[x]])
    for (z in 1:length(zoomNames)) {
      pdf(paste0(figFolder, '/', expName, '_pattern_', N, '_', nseed, zoomNames[z], '.pdf'))
      iter <- 0
      for (method in theseMethods) {
        iter <- iter + 1
        b <- 'biased'
        pr.bias <- pr.curve(scores.class0=stats[[method]][[b]]$posscores.pattern, scores.class1=stats[[method]][[b]]$negscores.pattern, rand.compute=TRUE, curve=TRUE)
        if (iter == 1)
          plot(pr.bias, main=paste0('PR'), lwd=2, xaxs="i", yaxs="i", color=colors[iter], pch=0, xlim=c(0,zoomXMax[z]), rand.plot=TRUE, auc.main=FALSE)
        else 
          plot(pr.bias, add=TRUE, main=paste0('PR'), lwd=2, xaxs="i", yaxs="i", color=colors[iter], pch=0, xlim=c(0,zoomXMax[z]), rand.plot=TRUE, auc.main=FALSE)
        if (addMarkerToCurve) plotPRCurveMarker(pr.bias, value=scoreValueForMarker, lwd=markerLwd, pch=markerPch, col=colors[iter])
        b <- 'unbiased'
        pr.unbias <- pr.curve(scores.class0=stats[[method]][[b]]$posscores.pattern, scores.class1=stats[[method]][[b]]$negscores.pattern, rand.compute=TRUE, curve=TRUE)
        plot(pr.unbias, add=TRUE, main=paste0('PR'), lwd=2, xaxs="i", yaxs="i", color=colors[iter], lty='dashed')
      }
      if (length(theseMethods) > 0) {
        legend("bottomright",legend=methodsPretty,col=colors[1:length(methodsPretty)],pch=15,bty="n",pt.cex=2,cex=1.2,text.col="black",horiz=F,inset=c(0.01, 0.01))
      }
      dev.off()
    }
  }
}

cat('*** DONE ***\n')