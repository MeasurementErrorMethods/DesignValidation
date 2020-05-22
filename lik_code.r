library(TwoPhaseReg)
library(tidyverse)
library(devtools)
library(splines)
library(MASS)

######################################
set.seed(210818)

Nsim <- 250
n  <- 2000
n2 <- 500
nsieve <- 25
res <- list(0)
print_step <- 2
######################################

beta <- c(1, 1, 1)
e_U <- c(sqrt(.5), sqrt(.5))
#######################################



############################################################
###################  Analysis  #############################
############################################################

for (nsim in 1:Nsim) {
  
  simZ   <- rbinom(n, 1, .5)
  simX   <- (1-simZ)*rnorm(n, 0, 1) + simZ*rnorm(n, 0.5, 1)
  epsilon <- rnorm(n, 0, 1)
  simY    <- beta[1] + beta[2]*simX + beta[3]*simZ + epsilon
  simX_tilde <- simX + rnorm(n, 0, e_U[1]*(simZ==0) + e_U[2]*(simZ==1))
  data <- data.frame(Y_tilde=simY, X_tilde=simX_tilde, Y=simY, X=simX, Z=simZ)
  
  indX2_0 <- which(data$Z == 0)
  nX2_0   <- length(indX2_0)
  indX2_1 <- which(data$Z == 1)
  nX2_1   <- length(indX2_1)
  
  Bspline_0 <- bs(data$X_tilde[indX2_0], df=nsieve, degree=3, Boundary.knots=range(data$X_tilde[indX2_0]), intercept=TRUE)
  Bspline_1 <- bs(data$X_tilde[indX2_1], df=nsieve, degree=3, Boundary.knots=range(data$X_tilde[indX2_1]), intercept=TRUE)
  Bspline <- matrix(0, n, 2*nsieve)
  Bspline[indX2_0,1:nsieve]              <- Bspline_0
  Bspline[indX2_1,(nsieve+1):(2*nsieve)] <- Bspline_1
  colnames(Bspline) <- paste("bs", 1:(2*nsieve), sep="")
  
  data <- cbind(data, Bspline)
  
  ##### Designs
  ## SRS
  id_phase2 <- c(sample(n, n2))
  dat_srs   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[1]]  <- smle_MEXY(Y="Y", X="X", Z="Z", Y_tilde="Y_tilde", X_tilde="X_tilde",
                         Bspline=colnames(Bspline), data=dat_srs, noSE=TRUE)
  
  ## ODS
  order_Y   <- order(data$Y_tilde)
  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  dat_ods   <- data %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res[[2]]  <- smle_MEXY(Y="Y", X="X", Z="Z", Y_tilde="Y_tilde", X_tilde="X_tilde",
                         Bspline=colnames(Bspline), data=dat_ods, noSE=TRUE)

  if (nsim == 1){
    results_est  <- list()
    results_conv <- matrix(NA, nrow=Nsim, ncol=length(res))
    for (i in 1:length(res)) {
      results_est[[i]]  <- matrix(NA, nrow=Nsim, ncol=length(beta))
      colnames(results_est[[i]]) <- c("intercept", "X", "Z")
    }
  }
  
  ## Save and print
  for (i in 1:length(res)){
    results_est[[i]][nsim,] <- t(res[[i]]$coefficients[,1])
    results_conv[nsim,i] <- res[[i]]$converge
  }
  
  if (nsim %% print_step==0) {
    print(paste(nsim, "replicates done."))
    for (i in 1:length(res)) {
      bias <- sapply(results_est, FUN = function(x) colMeans(x[1:nsim,], na.rm=TRUE)) -
        matrix(rep(beta), length(results_est), nrow=length(beta))
      se   <- sapply(results_est, FUN = function(x) apply(x[1:nsim,], 2, var, na.rm=TRUE))
    }
    mse <- t(bias)^2 + t(se)
    cat('bias'); print(round(bias*1000,3))
    cat('se');   print(round(se*1000,3))
    cat('\n')
    cat('conv'); print(colMeans(results_conv[1:nsim,]))
    cat('\n')
  }
}

