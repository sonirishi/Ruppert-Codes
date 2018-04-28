rm(list=ls(all=T))
setwd("E:/Documents/Practice")

library(dplyr)
library(Ecdat)
library(data.table)
library(fGarch)
library(MASS)

data("CRSPday")

print(cov(CRSPday[,4:7]))

covar_matrix <- matrix(c(2,1.47,0,1.47,3,0,0,0,5),3,3)

weight_vector <- c(1,1,0.5)

var_linear <- t(weight_vector) %*% covar_matrix %*% weight_vector

print(var_linear)

################## Covariance between benchmark and portfolio

W1 <- c(1/3,1/3,1/3)

W2 <- c(1/2,1/2,0)

covar_bet_port <- t(W1) %*% covar_matrix %*% W2

print(covar_bet_port)

####################################

W <- matrix(c(1/3,1/2,1/3,1/2,1/3,0),2,3)

print(W %*% covar_matrix %*% t(W))  ### Covariance matrix of the port and benchmark

#######################

library(readxl)

data <- readxl::read_xlsx("sp_monthly.xlsx")

var_montly <- var(data$returns)

risk_annualized_mont <- sqrt(12*var_montly)

data_1 <- readxl::read_xlsx("sp_yearly.xlsx")

var_yrly <- var(data_1$returns)

risk_annual <- sqrt(var_yrly)

print(sd(data_1$returns[2:6]))

#print(sd(data_1$returns))

data_2 <- readxl::read_xlsx("sp_daily.xlsx")

data_2 <- data_2 %>% mutate(lag_price = lag(price), return = (price - lag_price)/lag_price)

var_daily <- var(data_2$return[7141:nrow(data_2)])

risk_annualized_daily <- sqrt(252*var_daily)

#################

library(mvtnorm)

pairs(CRSPday[,4:7])

multnorm <- as.data.frame(mvrnorm(n=1000, mu = c(0,0), Sigma = matrix(c(1,0.5,0.5,1),2,2)))

colnames(multnorm) <- c("x","y")

multnorm_1 <- as.data.frame(mvrnorm(n=1000, mu = c(0,0), Sigma = matrix(c(1,-0.5,-0.5,1),2,2)))

plot(multnorm_1)

multi_t <- mvtnorm::rmvt(2500,df=3,delta=c(0,0),sigma = matrix(c(1,0,0,1),2,2))

plot(multi_t)
cor(multi_t)

library(fGarch)

t_1 <- fGarch::rstd(2500,0,1,3)
t_2 <- fGarch::rstd(2500,0,1,3)

t_full <- cbind(t_1,t_2)

plot(t_full)

midcap <- read.csv("midcapD.csv")

pairs(midcap[,c(2:7)])

#################################### Multivariate T Fit

library(MASS)

mfit_t <- cov.trob(CRSPday[,c(4:7)])

mtfit <- function(series){
  
  start <- append(append(as.vector(apply(series,2,mean)), as.vector(cov(series))),c(2.1))
  
  loglik <- function(par) {
    
    mean <- par[1:4]
    covar <- matrix(par[5:20],4,4)
    df <- par[21]
    
    f <- -sum(log(dmvt(series, mean, covar, df)))
    f
  }
  
  cov1 <- c(1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,
            1e-10,1e-10,1e-10,1e-10,1e-10,1e-10)
  cov2 <- c(1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000)
  
  lower <- append(append(c(-1000,-1000,-1000,-1000),cov1),2.1)
  upper <- append(append(c(1000,1000,1000,1000),cov2),60)
  results <- optim(start, loglik, method = "L-BFGS-B",lower=lower,upper=upper,hessian=T)
  return(results)
}

fit_mine <- mtfit(CRSPday[,c(4:7)])
