setwd('C:/Users/amorigg1/Documents/VCCC')
library(data.table)
library(tidyverse)
library(mvtnorm)
library(survey)
library(MASS)

source('C:/Users/amorigg1/Documents/VCCC/raking_ipw_vccc.R')


art.v = read.csv("shepherd_valid_art_13feb2014.csv", as.is=TRUE)
art.nv = read.csv("shepherd_non_valid_art_13feb2014.csv", as.is=TRUE)
lab = read.csv("shepherd_labs_13feb2014.csv", as.is=TRUE)
demo = read.csv("shepherd_demo_13feb2014.csv", as.is=TRUE)

#### START: create ART #################################################
# exclude people with missing AGE_AT_MED_START, becuase those people might start ART before enrollment
art.v.ex.id = unique(art.v$CFAR_PID[which(is.na(art.v$AGE_AT_MED_START))])
art.v1 = art.v[-which(art.v$CFAR_PID %in% art.v.ex.id),]
art.v2 <- art.v1 %>% left_join(demo, by='CFAR_PID') %>% group_by(CFAR_PID) %>% 
  mutate(YEARS_IN_STUDY = YEAR_OF_LAST_VISIT - YEAR_OF_ENROLLMENT,
         MORE_ONE_YEAR = max(YEARS_IN_STUDY)) %>% arrange(AGE_AT_MED_START) %>% slice(1)

art.nv.ex.id = unique(art.nv$CFAR_PID[which(is.na(art.nv$AGE_AT_MED_START))])
art.nv1 = art.nv[-which(art.nv$CFAR_PID %in% art.nv.ex.id),]
art.nv2 <- art.nv1 %>% left_join(demo, by='CFAR_PID') %>% group_by(CFAR_PID) %>% 
  mutate(YEARS_IN_STUDY = YEAR_OF_LAST_VISIT - YEAR_OF_ENROLLMENT,
         MORE_ONE_YEAR = max(YEARS_IN_STUDY)) %>% arrange(AGE_AT_MED_START) %>% slice(1)

art.cb = merge(art.v2, art.nv2, by="CFAR_PID", suffixes=c("_v","_nv"))
#### END: create ART ###################################################

# retain lab values for CD4 and VL
lab_new <- lab %>% group_by(CFAR_PID) %>%  mutate(AGE_DIFF = AGE_AT_RESULT_DATE - AGE_AT_RESULT_DATE[1],
                                                  MAX_DIFF = max(AGE_DIFF)) #%>% filter(MAX_DIFF > .75)

#### START: map CD4 ####################################################
art.cb$CD4_COUNT_v = NA
ids_both <- intersect(art.cb$CFAR_PID, lab_new$CFAR_PID)

art.cb_temp  <- art.cb  %>% filter(CFAR_PID %in% ids_both)
lab_new_temp <- lab_new %>% filter(CFAR_PID %in% ids_both) %>% filter(testName == 'CD4 COUNT')
lab_new_temp_VL <- lab_new %>% filter(CFAR_PID %in% ids_both) %>% filter(testName != 'CD4 COUNT' & testName != 'CD4 PERCENT')

art.cb$VL_COUNT_BSL_v = NA
art.cb$VL_COUNT_BSL_nv = NA
art.cb$CD4_COUNT_BSL_v = NA
art.cb$CD4_COUNT_1Y_v = NA

for (i in 1:nrow(art.cb)) {
  lab.ids = which(lab_new_temp$CFAR_PID == art.cb_temp$CFAR_PID[i])
  ages.cd4 = lab_new_temp$AGE_AT_RESULT_DATE[lab.ids]
  cd4 = lab_new_temp$RESULT_NUMERIC[lab.ids]

  age.s.v = art.cb_temp$AGE_AT_MED_START_v[i]
  age.s.nv = art.cb_temp$AGE_AT_MED_START_nv[i]
  if (length(lab.ids) > 0) {
    close.ids.v = which(ages.cd4 >= age.s.v-180/365.25 & ages.cd4 <= age.s.v+30/365.25)
    if (length(close.ids.v) > 0) {
      id.v = close.ids.v[which.min(abs(ages.cd4[close.ids.v]-age.s.v))]
      art.cb$CD4_COUNT_BSL_v[i] = cd4[id.v]
    }
    
    close.ids.v = which(ages.cd4 >= age.s.v + 1 - 180/365.25 & ages.cd4 <= age.s.v + 1 + 30/365.25)
    if (length(close.ids.v) > 0) {
      id.v = close.ids.v[which.min(abs(ages.cd4[close.ids.v]-age.s.v))]
      art.cb$CD4_COUNT_1Y_v[i] = cd4[id.v]
    }
  }
}

for (i in 1:nrow(art.cb)) {
  lab.ids = which(lab_new_temp_VL$CFAR_PID == art.cb_temp$CFAR_PID[i])
  ages.vl = lab_new_temp_VL$AGE_AT_RESULT_DATE[lab.ids]
  vl = lab_new_temp_VL$RESULT_NUMERIC[lab.ids]
  
  age.s.v = art.cb_temp[art.cb_temp$CFAR_PID == art.cb_temp$CFAR_PID[i],]$AGE_AT_MED_START_v
  age.s.nv = art.cb_temp[art.cb_temp$CFAR_PID == art.cb_temp$CFAR_PID[i],]$AGE_AT_MED_START_nv
  if (length(lab.ids) > 0) {
    close.ids.v = which(ages.vl >= age.s.v-180/365.25 & ages.vl <= age.s.v+30/365.25)
    if (length(close.ids.v) > 0) {
      id.v = close.ids.v[which.min(abs(ages.vl[close.ids.v]-age.s.v))]
      art.cb$VL_COUNT_BSL_v[i] = vl[id.v]
    }
    
    close.ids.nv = which(ages.vl >= age.s.nv-180/365.25 & ages.vl <= age.s.nv+30/365.25)
    if (length(close.ids.nv) > 0) {
      id.nv = close.ids.nv[which.min(abs(ages.vl[close.ids.nv]-age.s.nv))]
      art.cb$VL_COUNT_BSL_nv[i] = vl[id.nv]
    }
  }
}

ex.ids = which(is.na(art.cb$CD4_COUNT_BSL_v) | is.na(art.cb$CD4_COUNT_1Y_v) | is.na(art.cb$VL_COUNT_BSL_v) | is.na(art.cb$VL_COUNT_BSL_nv))
art.cb1 = art.cb[-ex.ids,]

art.cb1$CD4_COUNT_BSL_sqrt_v = sqrt(art.cb1$CD4_COUNT_BSL_v)
art.cb1$CD4_COUNT_1Y_sqrt_v = sqrt(art.cb1$CD4_COUNT_1Y_v)

art.cb1$VL_COUNT_BSL_LOG_v = log(art.cb1$VL_COUNT_BSL_v)
art.cb1$VL_COUNT_BSL_LOG_nv = log(art.cb1$VL_COUNT_BSL_nv)

mean(art.cb1$VL_COUNT_BSL_LOG_v != art.cb1$VL_COUNT_BSL_LOG_nv)
cor(art.cb1$VL_COUNT_BSL_LOG_v, art.cb1$CD4_COUNT_BSL_sqrt_v)
cor(art.cb1$VL_COUNT_BSL_LOG_v, art.cb1$CD4_COUNT_1Y_sqrt_v)
#### END: map CD4 ######################################################

#### full-cohort analysis #########################################
data <- art.cb1 %>% mutate(Y = CD4_COUNT_1Y_sqrt_v, X = VL_COUNT_BSL_LOG_v, X_tilde = VL_COUNT_BSL_LOG_nv, Z = CD4_COUNT_BSL_sqrt_v, W = SEX_v) %>%
  ungroup() %>% dplyr::select(Y, X, X_tilde, Z, W)
ods_input <- list(q1 = c(.2, .85), prop_ods = c(.3, .4, .3),
                  nstrata = 3, mx = 0, zrange = 1, zprob =.5,
                  e_U = c(sqrt(.5),sqrt(.5)), calculate_if_only = TRUE)
run <- sim_function(data, n=nrow(data), n2=250, beta=1, ods_input=ods_input, kk=1)

## Fit true model
fit <- lm(Y ~ X + Z + W, data=data)
res_full  <- c(summary(fit)$coef[2,1], summary(fit)$coef[2,2]^2)

## Fit error-prone model
fit_error <- lm(Y ~ X_tilde + Z + W, data=data)
res_error <- c(summary(fit_error)$coef[2,1], summary(fit_error)$coef[2,2]^2)





#### spmle ##################################################

library(TwoPhaseReg)
library(splines)
sdXT <- residuals(lm(X_tilde ~ Y + Z + W, data=data))
sdX  <- c(sd(sdXT[data$W==1]), sd(sdXT[data$W==2]))

spmlefunc <- function(data, nsieve, nsieve2, nsieve3, dg, n=nrow(data), n2=250, sdX){
  
  indX2_0 <- which(data$W == 1)
  nX2_0   <- length(indX2_0)
  indX2_1 <- which(data$W == 2)
  nX2_1   <- length(indX2_1)
  
  Bspline_0 <- bs(data$X_tilde[indX2_0], df=nsieve, degree=dg, Boundary.knots=range(data$X_tilde[indX2_0]), intercept=TRUE)
  Bspline_1 <- bs(data$X_tilde[indX2_1], df=nsieve, degree=dg, Boundary.knots=range(data$X_tilde[indX2_1]), intercept=TRUE)
  Bspline <- matrix(0, nrow(data), 2*nsieve)
  Bspline[indX2_0,1:nsieve]              <- Bspline_0
  Bspline[indX2_1,(nsieve+1):(2*nsieve)] <- Bspline_1
  colnames(Bspline) <- paste("bs", 1:(2*nsieve), sep="")
  
  
  data2 <- data.frame(data, Bspline)
  id_phase2 <- c(sample(n, n2))
  dat_srs <- data2 %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res1    <- smle_MEXY(Y="Y", X="X", Z=c("Z", "W"), Y_tilde="Y", X_tilde="X_tilde",
                       Bspline=colnames(Bspline), data=dat_srs)
  
  id_phase2  <- c(sample((1:n)[data$W==1], n2/2), sample((1:n)[data$W==2], n2/2))
  dat_ssrs <- data2 %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res2     <- smle_MEXY(Y="Y", X="X", Z=c("Z","W"), Y_tilde="Y", X_tilde="X_tilde",
                        Bspline=colnames(Bspline), data=dat_ssrs)
  ## RS
  #  RS <- resid(lm(Y ~ Z, data=data))
  RS <- resid(lm(Y ~ X_tilde + Z, data=data2))
  order_RS  <- order(RS)
  id_phase2 <- c(order_RS[1:(n2/2)], order_RS[(n-n2/2+1):n])
  dat_rs <- data2 %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res4   <- smle_MEXY(Y="Y", X="X", Z=c("Z","W"), Y_tilde="Y", X_tilde="X_tilde",
                      Bspline=colnames(Bspline), data=dat_rs)
  
  ## WRS
  order_WRS <- order(RS*ifelse(data$Z==0, sdX[1], sdX[2]))
  id_phase2 <- c(order_WRS[1:(n2/2)], order_WRS[(n-n2/2+1):n])
  dat_wrs <- data2 %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res5    <- smle_MEXY(Y="Y", X="X", Z=c("Z","W"), Y_tilde="Y", X_tilde="X_tilde",
                       Bspline=colnames(Bspline), data=dat_wrs)
  
  
  Bspline_0 <- bs(data$X_tilde[indX2_0], df=nsieve2, degree=dg, Boundary.knots=range(data$X_tilde[indX2_0]), intercept=TRUE)
  Bspline_1 <- bs(data$X_tilde[indX2_1], df=nsieve3, degree=dg, Boundary.knots=range(data$X_tilde[indX2_1]), intercept=TRUE)
  Bspline <- matrix(0, n, nsieve2+nsieve3)
  Bspline[indX2_0,1:nsieve2]                     <- Bspline_0
  Bspline[indX2_1,(nsieve2+1):(nsieve2+nsieve3)] <- Bspline_1
  colnames(Bspline) <- paste("bs", 1:(nsieve2+nsieve3), sep="")
  data3 <- cbind(data, Bspline)
  
  order_Y   <- order(data2$Y)
  id_phase2 <- c(order_Y[1:(n2/2)], order_Y[(n-n2/2+1):n])
  dat_ods <- data3 %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res3    <- smle_MEXY(Y="Y", X="X", Z=c("Z", "W"), Y_tilde="Y", X_tilde="X_tilde",
                       Bspline=colnames(Bspline), data=dat_ods)
  
  ## SFS
  order_SFS <- order(RS*data3$X_tilde)
  id_phase2 <- c(order_SFS[1:(n2/2)], order_SFS[(n-n2/2+1):n])
  dat_sfs <- data3 %>% mutate(R = ifelse(c(1:n) %in% id_phase2, 1, 0), X = ifelse(R==0, NA, X))
  res6    <- smle_MEXY(Y="Y", X="X", Z=c("Z","W"), Y_tilde="Y", X_tilde="X_tilde",
                       Bspline=colnames(Bspline), data=dat_sfs)
  ests <- rbind(res1$coefficients[2,1], res2$coefficients[2,1], res3$coefficients[2,1],
                res4$coefficients[2,1], res5$coefficients[2,1], res6$coefficients[2,1])
  vars <- rbind(res1$coefficients[2,2], res2$coefficients[2,2], res3$coefficients[2,2],
                res4$coefficients[2,2], res5$coefficients[2,2], res6$coefficients[2,2])
  
  return(list(ests=ests, vars=vars))
}

dg <- 3
nsieve  <- 15
nsieve2 <- 15
nsieve3 <- 15
a1 <- spmlefunc(data, nsieve, nsieve2, nsieve3, dg, n=nrow(data), n2=250, sdX=sdX)

dg <- 3
nsieve  <- 20
nsieve2 <- 10
nsieve3 <- 20
a2 <- spmlefunc(data, nsieve, nsieve2, nsieve3, dg, n=nrow(data), n2=500, sdX=sdX)

2*1.96*a1$vars

res_VL_SPMLE_250 <- cbind(a1$ests, a1$vars^2)
res_VL_SPMLE_500 <- cbind(a2$ests, a2$vars^2)






#### ipw and raking #########################################

Nsim <- 100
res_VL1 <- var_VL1 <- matrix(NA, Nsim, ncol(run$res))
res_VL2 <- var_VL2 <- matrix(NA, Nsim, ncol(run$res))
for (i in 1:Nsim) {
  res_VL_temp <- sim_function(data, n=nrow(data), n2=250, beta=1, ods_input=ods_input, kk=1)
  res_VL1[i,] <- res_VL_temp$res[2,]
  var_VL1[i,] <- res_VL_temp$var[2,]
  
  res_VL_temp <- sim_function(data, n=nrow(data), n2=750, beta=1, ods_input=ods_input, kk=1)
  res_VL2[i,] <- res_VL_temp$res[2,]
  var_VL2[i,] <- res_VL_temp$var[2,]
  
  if (i %% 25 == 0) {
    cat('iter = ',i,'\n')
    cat('var_wgt_avg (250 sim) = ', round(apply(var_VL1[1:i,], 2, mean), 5), '\n')
    cat('var_wgt_avg (500 sim) = ', round(apply(var_VL2[1:i,], 2, mean), 5), '\n')
  }
}
res_ipw_250 <- cbind(apply(res_VL1, 2, mean)[1:ncol(res_VL1) %% 2 == 0], apply(var_VL1, 2, mean)[1:ncol(var_VL1) %% 2 == 0])
res_rak_250 <- cbind(apply(res_VL1, 2, mean)[1:ncol(res_VL1) %% 2 != 0], apply(var_VL1, 2, mean)[1:ncol(var_VL1) %% 2 != 0])

res_ipw_500 <- cbind(apply(res_VL2, 2, mean)[1:ncol(res_VL2) %% 2 == 0], apply(var_VL2, 2, mean)[1:ncol(var_VL2) %% 2 == 0])
res_rak_500 <- cbind(apply(res_VL2, 2, mean)[1:ncol(res_VL2) %% 2 != 0], apply(var_VL2, 2, mean)[1:ncol(var_VL2) %% 2 != 0])

#### mi #########################################

set.seed(2108)
res_VL_MI_250 <- fun_mi(data, n=nrow(data), n2=1250, k=20)
res_VL_MI_500 <- fun_mi(data, n=nrow(data), n2=750, k=20)

#### combine #########################################

res_all_250 <- rbind(res_full, res_VL_MI_250, res_ipw_250, res_rak_250, res_VL_SPMLE_250)
res_all_250 <- cbind(res_all_250[,1], 2*1.96*sqrt(res_all_250[,2]))
rownames(res_all_250) <- c('full', c('srs-mi', 'ods-mi', 'rs-mi', 'wrs-mi', 'sfs-mi', 'sfs-mi-ipw'),
                           colnames(run$res)[1:ncol(run$res) %% 2 == 0],
                           colnames(run$res)[1:ncol(run$res) %% 2 != 0],
                           c('srs-spmle', 'ssrs-spmle', 'ods-spmle', 'rs-spmle', 'wrs-spmle', 'sfs-spmle'))

res_all_500 <- rbind(res_full, res_VL_MI_500, res_ipw_500, res_rak_500, res_VL_SPMLE_500)
res_all_500 <- cbind(res_all_500[,1], 2*1.96*sqrt(res_all_500[,2]))
rownames(res_all_500) <- c('full', c('srs-mi', 'ods-mi', 'rs-mi', 'wrs-mi', 'sfs-mi', 'sfs-mi-ipw'),
                           colnames(run$res)[1:ncol(run$res) %% 2 == 0],
                           colnames(run$res)[1:ncol(run$res) %% 2 != 0],
                           c('srs-spmle', 'ssrs-spmle', 'ods-spmle', 'rs-spmle', 'wrs-spmle', 'sfs-spmle'))


cbind(res_all_250, res_all_500)

### n = 250
res_250_mi  <- rbind(res_VL_MI_250[,1] - res_full[1], res_VL_MI_250[,2]/res_full[2])
res_250_ipw <- rbind(res_ipw_250[,1] - res_full[1], res_ipw_250[,2]/res_full[2])
res_250_rak <- rbind(res_rak_250[,1] - res_full[1], res_rak_250[,2]/res_full[2])

colnames(res_250_mi)  <- c('srs-mi', 'ods-mi', 'rs-mi', 'wrs-mi', 'sfs-mi', 'sfs-mi-ipw')
colnames(res_250_ipw) <- colnames(run$res)[1:ncol(run$res) %% 2 == 0]
colnames(res_250_rak) <- colnames(run$res)[1:ncol(run$res) %% 2 != 0]
rownames(res_250_mi) <- rownames(res_250_ipw) <- rownames(res_250_rak) <- c('bias', '95% CI width')

res_250_mi
res_250_ipw
res_250_rak

### n = 500

res_500_mi  <- rbind(res_VL_MI_500[,1] - res_full[1], res_VL_MI_500[,2]/res_full[2])
res_500_ipw <- rbind(res_ipw_500[,1] - res_full[1], res_ipw_500[,2]/res_full[2])
res_500_rak <- rbind(res_rak_500[,1] - res_full[1], res_rak_500[,2]/res_full[2])


colnames(res_500_mi)  <- c('srs-mi', 'ods-mi', 'rs-mi', 'wrs-mi', 'sfs-mi', 'sfs-mi-ipw')
colnames(res_500_ipw) <- colnames(run$res)[1:ncol(run$res) %% 2 == 0]
colnames(res_500_rak) <- colnames(run$res)[1:ncol(run$res) %% 2 != 0]
rownames(res_500_mi) <- rownames(res_500_ipw) <- rownames(res_500_rak) <- c('bias', '95% CI width')

res_500_mi
res_500_ipw
res_500_rak

###

tab_final <- round(rbind(
    res_250_mi[,-6],
    res_500_mi[,-6],
    res_250_ipw,
    res_500_ipw,
    res_250_rak,
    res_500_rak
  ),3)
library(xtable)
xtable(tab_final)





a <- matrix(rep(c(res_full[1], 0, res_full[1], 0), 10), ncol=10)
a0 <- round(t(cbind(res_all_250, res_all_500))[,c(8:17)] - a, 3)
a00 <- rbind(a0[,1:5], a0[,6:10])
xtable(a00, digits = 3)

a <- matrix(rep(c(res_full[1], 0, res_full[1], 0), 5), ncol=5)
a1 <- round(t(cbind(res_all_250, res_all_500))[,2:6] - a, 3)

a <- matrix(rep(c(res_full[1], 0, res_full[1], 0), 5), ncol=5)
a2 <- round(t(cbind(res_all_250, res_all_500))[,c(18,20:23)] - a, 3)

a3 <- rbind(a1, a2)
xtable(a3, digits = 3)



