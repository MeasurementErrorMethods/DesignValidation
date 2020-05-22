sim_function <- function(n=500, n2=250, beta, ods_input, vary_probs=FALSE, ...) {

  nODS  <- NeymanAllocation(data=data, q=ods_input$q1, beta=beta, type='Yds', n2=n2)$res2[,1]
  nODS  <- nODS[order(names(nODS))]
  nRS  <- NeymanAllocation(data=data, q=ods_input$q1, beta=beta, type='RS', n2=n2)$res2[,1]
  nRS  <- nRS[order(names(nRS))] 
#  nWRS <- NeymanAllocation(data=data, q=ods_input$q1, beta=beta, type='WRS', n2=n2)$res2[,1]
#  nWRS <- nWRS[order(names(nWRS))] 
  nIFS  <- NeymanAllocation(data=data, q=ods_input$q1, beta=beta, type='IF', n2=n2)$res2[,1]
  nIFS  <- nIFS[order(names(nIFS))]

  #sdX1 <- sd(residuals(lm(X ~ X_tilde + Z, data=data))[data$Z==0])
  #sdX2 <- sd(residuals(lm(X ~ X_tilde + Z, data=data))[data$Z==1])
  #sdX <- sdX[data$Z+1]
  
  #a <- func_cut(ifd, qs=c(.2,.85))
  #c(sum(a==1)*sd(ifd[a==1]), sum(a==2)*sd(ifd[a==2]), sum(a==3)*sd(ifd[a==3]))
  
  ######### Functions for sampling
  func_cut <- function(x, qs) 
    cut(x, breaks=c(-Inf, quantile(x, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
  
  func_samp <- function(x,y,z=nrow(data)) {
    samps <- c()
    for (jj in 1:length(table(x))){
      samps <- c(samps, sample((1:z)[x==jj], y[jj]))
    }
    samps
  }
  
  ##########
  
  #library(stratifyR)
  #ifs <- (data$X_tilde - mean(data$X_tilde))*residuals(lm(Y ~ X_tilde + Z, data=data))
  #opt_strata <- strata.data(ifs + 500, h = 3, n=n2)

  methods <- rep(c('srs', 'ods', 'rs', 'ifs'), each=2)
  
  id_phase2 <- strata <- Design <- w. <- Methods <- list()
  qods <- qrs <- qwrs <- qifs <- ods_input$q1
  
  ############## SRS
  #set.seed(1+kk)
  psrs <- which(methods=='srs')[1]
  id_phase2[[psrs]] <- id_phase2[[psrs+1]] <- c(sample(n, n2))
  strata[[psrs]]    <- strata[[psrs+1]]    <- rep(1, n)
  Design[[psrs]]    <- Design[[psrs+1]]    <- 'SRS'
  w.[[psrs]]        <- w.[[psrs+1]]        <- rep(1, n)
  Methods[[psrs]]   <- 'Raking'
  Methods[[psrs+1]] <- 'HT'

  ############## ODS sampling
  #set.seed(1+kk)
  pods <- which(methods=='ods')[1]
  ODS_strata <- func_cut(data$Y, qods)
  samps <- func_samp(ODS_strata, nODS)
  id_phase2[[pods]] <- id_phase2[[pods+1]] <- samps
  strata[[pods]]    <- strata[[pods+1]]    <- ODS_strata
  Design[[pods]]    <- Design[[pods+1]]    <- 'ODS'
  w.[[pods]]        <- w.[[pods+1]]        <- nODS/table(ODS_strata)
  Methods[[pods]]   <- 'Raking'
  Methods[[pods+1]] <- 'HT'

  ############## RS sampling
  #set.seed(1+kk)
  prs <- which(methods=='rs')[1]
  RS        <- resid(lm(Y ~ Z, data=data))
  RS_strata <- func_cut(RS, qrs)
  samps <- func_samp(RS_strata, nRS)
  id_phase2[[prs]] <- id_phase2[[prs+1]] <- samps
  strata[[prs]]    <- strata[[prs+1]]    <- RS_strata
  Design[[prs]]    <- Design[[prs+1]]    <- 'RS'
  w.[[prs]]        <- w.[[prs+1]]        <- nRS/table(RS_strata)
  Methods[[prs]]   <- 'Raking'
  Methods[[prs+1]] <- 'HT'

  ############## IF sampling
  #set.seed(1+kk)
  pifs <- which(methods=='ifs')[1]
  IFS  <- data$X_tilde*RS
  IFS  <- dfbeta(lm(Y ~ X_tilde + Z, data=data))[,2]
  IFS_strata <- func_cut(IFS, qifs)
  samps <- func_samp(IFS_strata, nIFS)
  id_phase2[[pifs]] <- id_phase2[[pifs+1]] <- samps
  strata[[pifs]]    <- strata[[pifs+1]]    <- IFS_strata
  Design[[pifs]]    <- Design[[pifs+1]]    <-'IFS'
  w.[[pifs]]        <- w.[[pifs+1]]        <- nIFS/table(IFS_strata)
  Methods[[pifs]]   <- 'Raking'
  Methods[[pifs+1]] <- 'HT'

  ############# Run calibration
  
  res <- matrix(NA, 3, length(methods))
  for (jj in 1:length(strata)){
  #  browser()
    res[,jj] <- coef(cal_function(data=data, id_phase2=id_phase2, i=jj, strata=strata, Methods=Methods, w.=w.)$res)[,1]   
  }

  ############# Combine results and return
  colnames(res) <- paste(methods, rep(c('raking', 'ht'), length(methods)), sep='-')[1:length(methods)]
  res
}



res_function <- function(res) {
  
  #### Getting the results
  est  <- lapply(1:ncol(res[[1]]), FUN=function(i) t(sapply(res, FUN = function(x) x[,i])))
  bias <- sapply(est, colMeans) - matrix(rep(beta, ncol(res[[1]])), ncol=ncol(res[[1]]))
  se   <- sapply(est, FUN=function(x) apply(x, 2, var))
  mse  <- bias^2 + se
  
  res_list <- list(bias=bias, se=se, mse=mse)
  res_list <- lapply(res_list, function(x) {colnames(x) <- colnames(res[[1]]); round(x*100,5)})
  
  cat('Results multiplied by 100')
  print(res_list) 
}

cal_function <- function(data, id_phase2, i, strata, Methods, w.){
  
  ###### Get specific variables: phase-2 indicator; variable for stratification
  id_phase2 <- id_phase2[[i]]
  data$strata <- strata[[i]]
  Method <- Methods[[i]]
  
  ###### Not selected for phase-2
  data$R <- 0
  data$R[id_phase2] <- 1
  data$X[data$R==0] <-  NA
  
  ###### Start calibration: impute missing variables
  data$id <- 1:nrow(data)
  
  ###### Get influence functions evaluated at imputed values
  inffun  <- dfbeta(lm(Y ~ X_tilde + Z, data=data))
  colnames(inffun) <- paste("if", 1:ncol(inffun), sep="")
  data_if <- cbind(data, inffun)
  
  wgt <- 1/(w.[[i]][data$strata])
  if (all(is.na(data$strata)))
    wgt <- rep(1, nrow(data))
  
  data_if$wgt <- wgt
  
  ######
  if_design <- twophase(id = list(~id, ~id), subset = ~(R==1), weights = list(NULL, ~ wgt), data = data_if,  method='approx')
  
  if (Method=='HT')     if_cal <- if_design
  if (Method=='Raking') {
      if_cal <- calibrate(if_design, phase=2, calfun="raking", formula=~if1+if2+if3 + Z*X_tilde, data=data_if)
  }
  
  res <- summary(svyglm(Y ~ X+Z, design=if_cal))
  res_lm <- NULL
  res_lm <- summary(glm(Y ~ X+Z, data=data_if[!is.na(data$X),], weights=wgt))
  return(list(res=res, res_lm=res_lm))
}





NeymanAllocation <- function(data,qs,beta=c(1,1,1),type,binary=FALSE,n=1e6,missclas_prop=0.05,n2=500) {
  
  #data <- gen_data(n, beta=beta)
  
  if (!binary) {
    
    if (type=='Yds') {
      Y1 <- cut(data$Y, breaks=c(-Inf, quantile(data$Y, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='Xds') {
      Y1 <- cut(data$X_tilde, breaks=c(-Inf, quantile(data$X_tilde, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='RS') {
      if (is.null(data$Z)) {
        rs <- resid(lm(Y ~ X_tilde, data=data))
      } else {
        rs <- resid(lm(Y ~ X_tilde + Z, data=data))
      }
      Y1 <- cut(rs, breaks=c(-Inf, quantile(rs, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='WRS') {
      #browser()
      sdX <- sd(residuals(lm(X ~ X_tilde + Z, data=data)))
      if (is.null(data$Z)) {
        wrs <- resid(lm(Y ~ X_tilde, data=data))*sdX
      } else {
        wrs <- resid(lm(Y ~ X_tilde + Z, data=data))*sdX
      }
      Y1 <- cut(wrs, breaks=c(-Inf, quantile(wrs, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='IF') {
      if (is.null(data$Z)) {
        #ifd <- data$X_tilde*resid(lm(Y ~ X_tilde, data=data))
        ifd <- dfbeta(lm(Y ~ X_tilde, data=data))[,2]
      } else {
        #ifd <- data$X_tilde*resid(lm(Y ~ X_tilde + Z, data=data))
        ifd <- dfbeta(lm(Y ~ X_tilde + Z, data=data))[,2]
      }
      Y1 <- cut(ifd, breaks=c(-Inf, quantile(ifd, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    if (type=='SR') {
      S1   <- data$X_tilde*resid(lm(Y ~ X_tilde, data=data))
      S2   <- resid(lm(Y ~ X_tilde, data=data))
      S3   <- resid(lm(Y ~ X_tilde, data=data))^2
      if (is.null(data$Z)) {
        S1S2 <- resid(lm(S1 ~ S2+S3))
      } else {
        S4   <- data$Z*resid(lm(Y ~ X_tilde, data=data))
        S1S2 <- resid(lm(S1 ~ S2+S3+S4))
      }
      Y1 <- cut(S1S2, breaks=c(-Inf, quantile(S1S2, probs=qs), Inf), labels=paste(1:(length(qs)+1), sep=','))
    }
    
    mu1 <- residuals(lm(Y ~ X_tilde, data=data))
    if (!is.null(data$Z))
      mu1 <- residuals(lm(Y ~ X_tilde + Z, data=data))
    
    #s1 <- sapply(1:(length(qs)+1), FUN=function(x) sd((data$X_tilde*mu1)[Y1==x]))
    s1 <- sapply(1:(length(qs)+1), FUN=function(x) sd(dfbeta(lm(Y ~ X_tilde + Z, data=data))[,2][Y1==x]))
    
    # Get/retun proportions
    res <- table(Y1)*s1/sum(table(Y1)*s1)
    
    res2 <- cbind(round(res*n2), table(Y1))
    res2 <- res2[order(res2[,1], decreasing=TRUE),]
    while(any(res2[,1] > res2[,2]))
    {
      pos <- which(res2[,1] > res2[,2]); pos <- pos[1]
      res2. <- matrix(c(res2[res2[,1] < res2[,2],]), ncol=2)
      res2[1,1] <- res2.[1,1] + (res2[pos,1] - res2[pos,2])
      res2[pos,1] <- res2[pos,2]
      res2[res2[,1] < res2[,2],] <- res2.
    }
  }
  if (binary){
    
    # Generate data and calculate Influence function
    data <- gen_data(N, beta=beta, binary=binary, missclas_prop=missclas_prop)
    mu <- sum(data$Y==1)/length(data$Y)
    U <- data$X_tilde*(data$Y - mu)
    
    # Get SD for each stratum
    Y1 <- data$Y==1
    Y0 <- data$Y==0
    
    # Get SD for each stratum
    s1 <- sd(U[Y1])
    s2 <- sd(U[Y0])
    
    # Get influence function
    mu1 <- sum(data$Y==1)/length(data$Y)
    mu0 <- sum(data$Y==1)/length(data$Y)
    U1 <- data$X_tilde[Y1]*(data$Y[Y1] - mu1)
    U0 <- data$X_tilde[Y0]*(data$Y[Y0] - mu0)
    
    # Get SD for each stratum
    s1 <- sd(U[Y1])
    s0 <- sd(U[Y0])
    
    # Get/return proportions
    p1 <- mu / (mu + (1-mu)*s2/s1) 
    res <- c(p1, (1 - p1))
  }
  return(list(res=res, res2=res2))
}




