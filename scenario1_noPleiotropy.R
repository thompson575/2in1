#---------------------------------------------------
# NO PLEIOTROPY SCENARIO 
#

library(MASS)
library(tidyverse)
library(AER)
library(metafor)
library(MendelianRandomization)

#-----------------------------------------------------
# Quick MR simulation with no pleiotropy
#
quickSimMR <- function(N,
                       af,
                       sigma,
                       alpha,
                       beta) {
  nSNP <- length(af)
  nX   <- 1
  nY   <- 1
  nV   <- 2
  # Sample 1
  E <- mvrnorm(n = N, mu=c(0, 0), Sigma=sigma)
  X <- E[, 1]
  X2 <- rnorm(n = N, mean=0, sd=sqrt(sigma[1, 1]))
  Y <- E[, 2]
  G <- G2 <- matrix(0, nrow=N, ncol=nSNP)
  
  varG <- 2*af*(1-af)  # var of each G
  
  varGX <- 0  # var of X due to G
  for( j in seq_len(nSNP) ) {
    G[, j]  <- rbinom(N, 2, af[j])
    G2[, j] <- rbinom(N, 2, af[j])
    X  <- X + alpha[j]*G[, j]
    X2 <- X2 + alpha[j]*G2[, j]
    varGX <- varGX + alpha[j]*alpha[j]*varG[j]
  }
  
  varX <- sigma[1,1] + varGX   # total variance of X
  
  Y <- Y + beta*X
  varXY <-  beta*beta*varX  # var of Y due to X
  
  varY <- sigma[2,2] + varXY    # total variance of X
  
  list(varG=varG, varX=varX, varY=varY, varGX=varGX, varXY=varXY,
       G=G, X=X, Y=Y, G2=G2, X2=X2, af=af, sigma=sigma, alpha=alpha, beta=beta)
}

#--------------------------------------------
# Function to control a set of simulations
#
runSim <- function(
  nReplicates,       # number of identical simulations
  seed,              # random number seed
  simulateData,      # function for simulating data
  analyseData,       # function for analysing the data
  simulateParams,    # parameters for simulateData
  analysisParams,    # parameters for analyseData
  resultsHeader,     # names for the evaluation
  resumeSim=1,       # simulation at which to resume
  resumeRep=1,       # replicate at which to resume
  resumeAna=1,       # analysis at which to resume
  resumeDumpFile='', # resume from a previous dump
  verbose=FALSE,     # whether to print results of each iteration
  dump=FALSE) {      # whether to dump results after each iteration
  startTime <- proc.time()
  # split simulation parameters
  sp <- splitParams(simulateParams)
  nSimParams <- length(simulateParams)
  # split analysis parameters
  ap <- splitParams(analysisParams)
  nAnaParams <- length(analysisParams)
  # expand to give all combinations of the varying parameters (may be empty)
  VPs <- expand.grid(c(sp$variableList))
  names(VPs) <- c(sp$namesVaryParams)
  VPa <- expand.grid(c(ap$variableList))
  names(VPa) <- c(ap$namesVaryParams)
  #---------------------------------------------
  # perform simulations
  #---------------------------------------------
  header <- c( sp$namesVaryParams, ap$namesVaryParams, resultsHeader )
  nsCombin <- max(1, nrow(VPs))
  naCombin <- max(1, nrow(VPa))
  if( resumeDumpFile != '' ) {
    load(file=resumeDumpFile)
    nR <- nrow(R)
    resumeAna <- 1 + (nR %% naCombin)
    resumeRep <- 1 + ((nR %/% naCombin) %% nReplicates)
    resumeSim <- 1 + ((nR %/% (naCombin*nReplicates)) %% nsCombin)
  }
  set.seed(seed)
  seeds <- floor( runif(nReplicates*nsCombin,1,1000000))
  s <- 0
  for( simrow in 1:nsCombin ) {
    if( simrow < resumeSim ) {
      s <- s + nReplicates
    } else {
      thisSimulateParams <- list()
      j <- 0
      for( i in seq_along(simulateParams) ) {
        if( sp$variableParam[i] == 0 ) {
          thisSimulateParams <- c(thisSimulateParams, simulateParams[i])
        } else {
          j <- j + 1
          thisSimulateParams <- c(thisSimulateParams, VPs[simrow,j])
        }
      }
      names(thisSimulateParams) <- names(simulateParams)
      for( rep in 1:nReplicates ) {
        s <- s + 1
        if( rep >= resumeRep ) {
          set.seed(seeds[s])
          D <- simulateData(thisSimulateParams)
          for( anarow in 1:naCombin ) {
            if( anarow >= resumeAna ) {
              thisAnalysisParams <- list()
              j <- 0
              for( i in seq_along(analysisParams) ) {
                if( ap$variableParam[i] == 0 ) {
                  thisAnalysisParams <- c(thisAnalysisParams, analysisParams[i])
                } else {
                  j <- j + 1
                  thisAnalysisParams <- c(thisAnalysisParams, VPa[anarow,j])
                }
              }
              names(thisAnalysisParams) <- names(analysisParams)
              analyseData(D,thisAnalysisParams) -> results
              newrow <- c(list(sim=rep,seed=seeds[s]),VPs[simrow,],VPa[anarow,],results)
              vrow <- newrow$sim
              for( j in 2:length(newrow) ) {
                vrow <- c(vrow,newrow[[j]])
              }
              if ( rep == resumeRep & simrow == resumeSim & anarow == resumeAna  ) {
                if( resumeDumpFile == '') {
                  R <- data.frame(newrow)
                  names(R) <- c("sim","seed",header)
                } else R <- rbind(R,vrow)
                resumeRep <- 1
                resumeAna <- 1
              } else {
                R <- rbind(R,vrow)
              }
              if( verbose ) {
                cat( sprintf("Sim %4.0f (%4.0f)     Replicate %6.0f (%6.0f)     Analysis %4.0f (%4.0f)\n",simrow,nsCombin,rep,nReplicates,anarow,naCombin))
                print( vrow )
              }
              if( dump ) {
                save(R, file="RunSimulationDump.RData")
              }
            }
          }
          
        }
      }
    }
  }
  list(nReplicates=nReplicates,
       globalSeed=seed,
       simulateData=simulateData,
       analyseData=analyseData,
       simulateParams=simulateParams,
       analysisParams=analysisParams,
       resume=c(resumeSim,resumeRep,resumeAna),
       time=proc.time()-startTime,
       results=R)
}

#---------------------------------------------
# splitParams: used by runSim
# distinguish between parameters
# that take a fixed value and those that are not fixed
# varying parameters have their values in a list
# For example,
#    list(s=10, list(m=c(2,3,4)), n=c(100,200))
# second parameter varies so 3 simulations run
# with m=2, m=3 and m=4
#---------------------------------------------
splitParams <- function( params ) {
  namesVaryParams <- NULL
  namesFixedParams <- NULL
  variableList <- list()
  fixedList <- list()
  nParams <- length(params)
  vj <- 0
  if ( nParams == 0 ) {
    variableParam <- 0
  } else {
    variableParam <- rep(0,nParams)
    for( i in 1:nParams ) {
      parameter <- params[[i]]
      if( class(parameter) == "list" ) {
        variableParam[i] <- 1
        namesVaryParams <- c(namesVaryParams, names(params)[i])
        vj <- vj + 1
        variableList[[vj]] <- as.vector(unlist(params[i]))
      } else {
        namesFixedParams <- c(namesFixedParams, names(params)[i])
        fixedList <- c(fixedList, params[[i]])
      }
    }
  }
  list( variableParam=variableParam,       # vector of 0=fixed, 1=variable
        fixedList=fixedList,               # list of the values of the fixed parameters
        namesFixedParams=namesFixedParams, # names of the fixed parameters
        variableList=variableList,         # list of vectors containing the variable params
        namesVaryParams=namesVaryParams )  # names of the variable parameters
}

#--------------------------------------------------------
# Function to Analyse the data
#
mrEstimators <- function(R, ap) {
  
  nSNP <- ncol(R$G)
  
  #-------------------------------------
  # OLS analysis
  #
  fit <- lm( R$Y ~ R$X)
  ols <- summary(fit)$coef[2,1]
  olsse <- summary(fit)$coef[2,2]
  
  #-------------------------------------
  # 2SLS analysis - Gold Standard
  #
  tsmodel <- ivreg(R$Y ~ R$X | R$G)
  b   <- summary(tsmodel)$coef[2,1]
  bse <- summary(tsmodel)$coef[2,2]
  
  # diagnostic p-values: weak instruments, Wu-Hausman, Sargan
  # p1 signif implies instruments are strong
  # p2 signif implies should not use OLS
  # p3 signif implies some instruments are not valid
  st <- summary(tsmodel, diagnostics=TRUE)
  p1 <- st$diagnostics[1,4]
  p2 <- st$diagnostics[2,4]
  p3 <- st$diagnostics[3,4]
  
  #-------------------------------------
  # Estimates for each SNP
  #
  bx <- bxse <- bx2 <- bxse2 <- by <- byse <- F <- F2 <- rep(0, nSNP)
  for( i in 1:nSNP ) {
    # GX sample 1
    fit <- lm( R$X ~ R$G[,i])
    bx[i] <- summary(fit)$coef[2,1]
    bxse[i] <- summary(fit)$coef[2,2]
    F[i] <- summary.aov(fit)[[1]][1,4]
    # GX sample 2
    fit <- lm( R$X2 ~ R$G2[,i])
    bx2[i] <- summary(fit)$coef[2,1]
    bxse2[i] <- summary(fit)$coef[2,2]
    F2[i] <- summary.aov(fit)[[1]][1,4]
    # GY sample 1
    fit <- lm( R$Y ~ R$G[,i])
    by[i] <- summary(fit)$coef[2,1]
    byse[i] <- summary(fit)$coef[2,2]
  }
  
  #-------------------------------------
  # I2 GX in sample 1
  #
  rm <- rma(yi=bx, sei=bxse, method='FE')
  i2gx <- rm$I2
  
  #-------------------------------------
  # one-sample: 
  #
  B <- mr_input(bx=bx, bxse=bxse, by=by, byse=byse)
  # FE model
  MR <- mr_ivw(B, model="fixed")
  bfe <- MR$Estimate
  bfese <- MR$StdError
  # RE model
  MR <- mr_ivw(B, model="random")
  bre <- MR$Estimate
  brese <- MR$StdError
  # summarise F statistics
  mF <- mean(F)
  nF <- sum(F<10)
  # median estimator
  MR <- mr_median(B)
  bmd <- MR$Estimate
  bmdse <- MR$StdError
  # modal estimator
  MR <- mr_mbe(B)
  mbe <- MR$Estimate
  mbese <- MR$StdError
  # egger estimator
  MR <- mr_egger(B)
  beg <- MR$Estimate
  begse <- MR$StdError.Est
  # evidence for pleiotropy: I2 and Q is Wald estimates
  bi <- by/bx
  bsei <- byse/abs(bx)
  # Code to omit extreme SE from metafor
  j <- 1:nSNP
  rm <- tryCatch( rma(yi=bi, sei=bsei, method='FE'),
                  error = function(e) {
                    bm <- bsei/median(bsei)
                    j <<- which(bm < 10 )
                    tryCatch(rma(yi=bi[j], sei=bsei[j], method='FE'),
                             error = function(e) {
                               return(NULL)})
                  })
  if( !is.null(rm) ) {
    i2 <- rm$I2
    qe <- rm$QE
    n3 <- length(j)
  } else {
    i2 <- qe <- n3 <- NA
  }
  
  #-------------------------------------
  # two-sample: 
  #
  B <- mr_input(bx=bx2, bxse=bxse2, by=by, byse=byse)
  # FE model
  MR <- mr_ivw(B, model="fixed")
  bfe2 <- MR$Estimate
  bfese2 <- MR$StdError
  # RE model
  MR <- mr_ivw(B, model="random")
  bre2 <- MR$Estimate
  brese2 <- MR$StdError
  # Summary of F statistics
  mF2 <- mean(F2)
  nF2 <- sum(F2<10)
  # median estimator
  MR <- mr_median(B)
  bmd2 <- MR$Estimate
  bmdse2 <- MR$StdError
  # modal estimator
  MR <- mr_mbe(B)
  mbe2 <- MR$Estimate
  mbese2 <- MR$StdError
  # egger estimator
  MR <- mr_egger(B)
  beg2 <- MR$Estimate
  begse2 <- MR$StdError.Est
  # evidence for pleiotropy: I2 and Q is Wald estimates
  bi <- by/bx2
  bsei <- byse/abs(bx2)
  # Code to omit extreme two-sample SEs from metafor
  j <- 1:nSNP
  rm <- tryCatch(rma(yi=bi, sei=bsei, method='FE'),
                 error = function(e) {
                   bm <- bsei/median(bsei)
                   j <<- which(bm < 10 )
                   tryCatch(rma(yi=bi[j], sei=bsei[j], method='FE'),
                            error = function(e) {
                              return(NULL)})
                 })
  if( !is.null(rm) ) {
    i22 <- rm$I2
    qe2 <- rm$QE
    n32 <- length(j)
  } else {
    i22 <- qe2 <- n32 <- NA
  }
  
  # qe <- sum( (by/bx - rep(rm$beta,50))^2/(byse*byse/(bx*bx)))
  # qe <- sum( ((bi - rep(rm$beta,50))/bsei)^2)
  
  #-------------------------------------------
  # R2 statistics
  #
  # GX
  fit <- lm( R$X ~ R$G[,1:nSNP])
  fx <- fitted(fit)
  r2x <- summary(fit)$r.squared
  # GY
  fit <- lm( R$Y ~ R$G[,1:nSNP])
  r2y <- summary(fit)$r.squared
  # G on E(X) & G  r2mg = r2y
  # fit <- lm( R$Y ~ fx + R$G[,1:nSNP])
  # r2mg <- summary(fit)$r.squared 
  # G on E(X)  - pleiotropy measure (r2y-r2m)/r2y
  fit <- lm(R$Y ~ fx)
  r2m <- summary(fit)$r.squared
  
  #Correlation
  rx = R$X - fx
  r <- cor(rx, R$Y-fitted(fit)-summary(fit)$coef[2,1]*rx)
  # R2 GX in second sample
  fit <- lm( R$X2 ~ R$G2[,1:nSNP])
  fx <- fitted(fit)
  r2x2 <- summary(fit)$r.squared

  list( nSNP, ols, olsse, mF, nF, mF2, nF2, 
        r2x, r2x2, r2y, r2m, i2gx, 
        b, bse, p1, p2, p3, 
        bfe, bfese, bre, brese, bfe2, bfese2, bre2, brese2,
        bmd, bmdse, mbe, mbese, beg, begse, i2, qe, r, 
        bmd2, bmdse2, mbe2, mbese2, beg2, begse2, i22, qe2, 
        n3, n32)
}

#-------------------------------------------------
# Simulation with no pleiotropy 
#
quickSimData01 <- function(sp) {
  AF <- runif(sp$nSNP, sp$af[1], sp$af[2])
  S  <- matrix( c(1, 3*sp$rho, 3*sp$rho, 9), nrow=2, ncol=2)
  A  <- sp$a[1] + 2*abs(AF-0.5)*sp$a[1] + rexp(sp$nSNP, rate=9/(sp$a[2]-sp$a[1]))
  quickSimMR(N=sp$n, af=AF, sigma=S, alpha=A, beta=sp$beta)
}


start_time <- Sys.time()
E <- runSim(nReplicates=100, 
            seed=12858,
            simulateData=quickSimData01, 
            analyseData=mrEstimators,
            simulateParams=list(n=100000, 
                                nSNP=100, 
                                af=c(0.01, 0.99), 
                                rho=list(-0.4, -0.2, 0, 0.2, 0.4), 
                                beta=list(0.0, 1.0), 
                                a = c(0.01, 0.1)),
            analysisParams=NULL,
            resultsHeader=c('nSNP', 'ols', 'olsse', 'mF', 'nF', 'mF2', 'nF2', 
                            'r2x', 'r2x2', 'r2y', 'r2m', 'i2gx', 
                            'b', 'bse', 'p1', 'p2', 'p3',
                            'bfe', 'bfese', 'bre', 'brese',
                            'bfe2', 'bfese2', 'bre2', 'brese2',
                            'bmd', 'bmdse', 'mbe',  'mbese', 'beg', 'begse', 'i2', 'qe', 'r',
                            'bmd2', 'bmdse2', 'mbe2',  'mbese2', 'beg2', 'begse2', 'i22', 'qe2', 
                            'n3', 'n32')
)
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

write_rds(E, "scenario1_results_100000.rds")
