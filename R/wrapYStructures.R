# Copyright (c) 2021-2022, Philip Versteeg (p.j.j.p.versteeg@gmail.com). All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
suppressMessages(library(pcalg))

ystructures <- function(
	data,
	systemVars,
	contextVars,
	alpha=1e-2,
	beta=NULL,
	test='gaussCItest',
	fullYStruc=FALSE,
	verbose=0,
	suffStatOracle=NULL) {

  # some checks
  if (test == 'oracle') 
    stopifnot(!is.null(suffStatOracle))
  if (is.null(beta)) 
    beta <- alpha
  stopifnot(alpha <= beta)

  pSys <- length(systemVars)
  pCon <- length(contextVars)
  arel <- matrix(0,p,p)
  colnames(arel) <- colnames(data)
  patterns <- list()

  # setup independence test
  S <- NULL
  if( test == 'gaussCItest' ) {
    indepTest<-gaussCItest
    suffStat<-list(C=cor(data),n=N,removeNAs=FALSE)
  } else if( test == 'oracle' ) {
    indepTest<-dsepTest
    suffStat<-suffStatOracle
    S <- suffStatOracle$S
  } else {
    stop('unknown test')
  }

  # compute some marg p-values
  p_ij <- matrix(0,p,p)
  for( i in 1:p ) {
    for( j in 1:p ) {
      if( i != j ) {
        p_ij[i,j] <- indepTest(i,j,c(S),suffStat)
      }
    }
  }
  # find ystructures
  for( i in 1:p ) {
    for( j in 1:p ) { 
      if( i != j ) {
        if ( !(p_ij[i,j] < alpha) ) {next}
        for( u in 1:p ) {
          if (u != i & u != j) {
            if( p_ij[u,j] < alpha ) {
              pval <- indepTest(u,j,c(S,i),suffStat)
              if (verbose >= 2) 
                cat('Testing indepTest(c,i,j): ',u,i,j,', pval=',pval,'\n', sep='')
              if( pval >= beta ) {
                if( verbose >= 1 )
                  cat('LCD triple <', u,',',i,',',j,'> found, pval=', p_ij[u,j],', searching for Y-structures..\n', sep='')
                bestScore <- 0
                numYstruct <- 0
                for ( v in 1:p) {
                  if ( v != i & v != j & v != u) {
                    pval_uv_i <- indepTest(u,v,c(S,i),suffStat)
                    if (verbose >= 2) 
                        cat('Testing indepTest(u,v,i): ',u,v,i,', pval=',pval_uv_i,'\n', sep='')
                    if ( pval_uv_i < alpha) {
                      if ( p_ij[u,v] >= beta) {
                        # full y-structures
                        if ( fullYStruc ) {
                          if( p_ij[v,j] < alpha ) { 
                            pval_vj_i <- indepTest(v,j,c(S,i),suffStat)
                            if ( !pval_vj_i >= beta ) {next}
                            if ( verbose >= 1 ) 
                              cat(' -- YStruc <', u, ',', i,',',j,',',v,'> found, pval=', p_ij[u,j], ' & pval=', p_ij[v,j], '\n', sep='')
                            # compute score
                            newScore <- min(-log(p_ij[u,j]), -log(p_ij[v,j]))
                            patterns[[length(patterns)+1]] <- c(i,j,u,v,newScore) 
                            bestScore <- max(bestScore,newScore)
                            numYstruct <- numYstruct + 1
                          }
                        # extended y-structures
                        } else {
                          if ( verbose >= 1 ) 
                            cat(' -- YStruc <', u, ',', i,',',j,',',v,'> found, pval=', p_ij[u,j], ' & pval=', p_ij[v,j], '\n', sep='')
                          # compute score
                          newScore <- min(-log(p_ij[u,j]), -log(p_ij[v,j]))
                          patterns[[length(patterns)+1]] <- c(i,j,u,v,newScore) 
                          bestScore <- max(bestScore,newScore)
                          numYstruct <- numYstruct + 1
                        }
                      }
                    }
                  }
                } 
                if (bestScore > 0) {
                  arel[i,j] <- max(arel[i,j], bestScore)
                  if( verbose >= 2 )
                    cat('Best YStruc found for <',i,',',j,'> out of ',  numYstruct, ' total, score=', bestScore, '\n', sep='')
                }
              }
            }
          }
        }
      }
    }
  }
  result<-list(p=p,systemVars=systemVars,contextVars=contextVars,labels=colnames(data),arel=arel,patterns=patterns)
  return(result)
}
