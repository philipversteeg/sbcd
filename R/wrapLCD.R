# Copyright (c) 2021-2022, Philip Versteeg (p.j.j.p.versteeg@gmail.com). All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
suppressMessages(library(pcalg))

lcd <- function(
  data,
	systemVars,
	contextVars,
	alpha=1e-2,
	beta=NULL,
	verbose=0,
	test='gaussCItest',
	suffStatOracle=NULL,
	conservative=FALSE) {

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
  if( pCon > 0 ) {
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
    p_ij <- matrix(0,pSys,pSys)
    p_ci <- matrix(0,pCon,pSys)
    for( i in 1:pSys ) {
      for( j in 1:pSys ) {
        if( i != j ) {
          p_ij[i,j] <- indepTest(i,j,c(S),suffStat)
        }
      }
      for( c in 1:pCon ) {
        p_ci[c,i] <- indepTest(pSys+c,i,c(S),suffStat)
      }
    }

    # find lcd triples
    for( i in 1:pSys ) {
      for( j in 1:pSys ) { 
        if( i != j ) {
          if ( !(p_ij[i,j] < alpha) ) {next}
          for( c in 1:pCon ) {
            if( p_ci[c,j] < alpha ) {
              pval <- indepTest(pSys+c,j,c(S,i),suffStat)
              if (verbose >= 2) 
                cat('Testing indepTest(pSys+c,i,j): ',pSys+c,i,j,', pval=',pval,'\n', sep='')
              if( pval >= beta ) {
                newp <- -log(p_ci[c,j])
                arel[i,j] <- max(arel[i,j],newp)
                patterns[[length(patterns)+1]] <- c(i,j,pSys+c,newp)
                if( verbose >= 1 ) cat('LCD triple <', pSys+c,',',i,',',j,'> found, pval=',p_ci[c,j],'\n',sep='')
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
