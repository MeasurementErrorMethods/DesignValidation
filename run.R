library(MASS)
library(survey)
library(ggplot2)
library(stratification)

run_func <- function(Nsim, N, n, ods_input, beta){

  Ntemp <- 1e6
  simZ   <- rbinom(Ntemp, ods_input$zrange, ods_input$zprob)
  simX   <- (1-simZ)*rnorm(Ntemp, 0, 1) + simZ*rnorm(Ntemp, 0.5, 1)
  epsilon <- rnorm(Ntemp, 0, 1)
  simY    <- beta[1] + beta[2]*simX + beta[3]*simZ + epsilon
  simX_tilde <- simX + rnorm(n, 0, ods_input$e_U[1]*(simZ==0) + ods_input$e_U[2]*(simZ==1))
  data <- data.frame(Y_tilde=simY, X_tilde=simX_tilde, Y=simY, X=simX, Z=simZ)
  
  sdX1 <- sd(residuals(lm(X ~ X_tilde + Z, data=data))[data$Z==0])
  sdX2 <- sd(residuals(lm(X ~ X_tilde + Z, data=data))[data$Z==1])
  sdX <<- round(c(sdX1, sdX2),3)
  
  single_run <- sim_function(n=N, n2=n, beta=beta, ods_input=ods_input, kk=1)[2,]
  
  res <- matrix(NA, Nsim, length(single_run))
  col_names <- names(single_run)
  set.seed(210818)
  for (kk in 1:Nsim){
    res[kk,] <- sim_function(n=N, n2=n, beta=beta, ods_input=ods_input, kk=kk)[2,]
    if (kk > 1){
      ratioCal  <- apply(res[1:kk,], 2, var)/apply(res[1:kk,], 2, var)[1]
      ratioCal[(1:ncol(res)) %% 2 == 0] <- NA
      ratioHT <- apply(res[1:kk,], 2, var)/apply(res[1:kk,], 2, var)[2]
      ratioHT[(1:ncol(res)) %% 2 == 1] <- NA
      res_bias_var  <- data.frame(variance_times = round(apply(res[1:kk,], 2, var)*1000,3),
                                  ratioHT_SRS = round(ratioHT,3), ratioCal_SRS = round(ratioCal,3),
                                  bias_times = round((apply(res[1:kk,], 2, mean) - beta[2])*1000,3))
      rownames(res_bias_var) <- col_names
      cat('\nk = ', kk)
      print(res_bias_var)
      
      if (kk %% 100 == 0){
        res2 <- as.vector(res[1:kk,])
        methods <- rep(col_names, each=kk)
        res2 <- data.frame(values = res2, methods = methods)
        plot_print <- ggplot(res2, aes(x=methods, y=values)) + geom_boxplot() + geom_point(size=.5)
        print(plot_print)
      }
    }
  }
  res
}


setwd('C:/Users/amorigg1/Desktop/Final_paper/IPW')

Nsim <- 2000
N <- 2000
n <- 500
beta  <- c(1,0,1)
ods_input <- list(q1 = c(.19, .81), prop_ods = c(.3, .4, .3),
                  nstrata = 3, mx = 0, zrange = 1, zprob =.5,
                  e_U = c(sqrt(.5),sqrt(.5)), calculate_if_only = TRUE)

file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,.5,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,1,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)



Nsim <- 2000
N <- 2000
n <- 500
beta  <- c(1,0,1)
ods_input <- list(q1 = c(.19, .81), prop_ods = c(.3, .4, .3),
                  nstrata = 3, mx = 0, zrange = 1, zprob =.5,
                  e_U = c(sqrt(1),sqrt(1)), calculate_if_only = TRUE)

file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,.5,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,1,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

Nsim <- 2000
N <- 2000
n <- 500
beta  <- c(1,0,1)
ods_input <- list(q1 = c(.19, .81), prop_ods = c(.3, .4, .3),
                  nstrata = 3, mx = 0, zrange = 1, zprob =.5,
                  e_U = c(sqrt(3),sqrt(3)), calculate_if_only = TRUE)

file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,.5,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,1,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)


Nsim <- 2000
N <- 2000
n <- 500
beta  <- c(1,0,1)
ods_input <- list(q1 = c(.19, .81), prop_ods = c(.3, .4, .3),
                  nstrata = 3, mx = 0, zrange = 1, zprob =.5,
                  e_U = c(sqrt(.5),sqrt(1)), calculate_if_only = TRUE)

file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,.5,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)

beta  <- c(1,1,1)
file_est  <- paste0('Design_Errors_X_3strata_Zprob',ods_input$zprob,'e_U1',ods_input$e_U[1],'e_U2',ods_input$e_U[2],"_beta",beta[2],"_betaZ",beta[3],"_N",N,"_n",n,".txt")
save2 <- paste('N = ', N, 'n = ', n, 'beta = ', paste(beta, collapse=" "))
res <- run_func(Nsim, N, n, ods_input, beta)
file.create(file_est)
write.table(save2, file = file_est,  sep = "\t", row.names = TRUE, col.names = TRUE, append=TRUE)
write.table(res, file = file_est, sep = "\t", row.names = TRUE, col.names = FALSE, append=TRUE)
