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
  
  loglik <- function(par) {
    
    mean <- par[1:4]
    rho12 <- par[5]; rho13 <- par[6]; rho14 <- par[7]
    rho23 <- par[8]; rho24 <- par[9]; rho34 <- par[10]
    var1 <- par[11]; var2 <- par[12]; var3 <- par[13]; var4 <- par[14]
    nu <- par[15]
    
    cov12 <- rho12*sqrt(var1*var2); cov13 <- rho13*sqrt(var1*var3); cov14 <- rho14*sqrt(var1*var4)
    cov23 <- rho23*sqrt(var2*var3); cov24 <- rho24*sqrt(var2*var4); cov34 <- rho34*sqrt(var3*var4)
    
    covar <- matrix(c(var1,cov12,cov13,cov14,cov12,var2,cov23,cov24,cov13,cov23,var3,cov34,
                      cov14,cov24,cov34,var4),4,4)
    
    f <- -sum(log(dmvt(x=series, delta=mean, sigma=covar, df=nu, log=FALSE)))
    f
  }
  
  cov1 <- c(1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,-0.1,-0.1,-0.1,-0.1)
  
  cov2 <- c(0.99,0.99,0.99,0.99,0.99,0.99,-0.1,0.1,0.1,0.1)
  
  lower <- append(append(c(-.02,-.02,-.02,-.02),cov1),2.1)
  upper <- append(append(c(0.02,0.02,0.02,0.02),cov2),15)
  
  start <- lower
  
  results <- optim(start, loglik, method = "L-BFGS-B",lower=lower,upper=upper,hessian=T)
  return(results)
}

fit_mine <- mtfit(CRSPday[,c(4:7)])

#####################

multi_t <- function(t,u){
  dmvt(cbind(t,u),df=4,sigma=matrix(c(2,1.1,1.1,2),2,2), log=F)
}

rang <- seq(-10,10,0.04)

contour(x=rang,y=rang,z=outer(rang,rang,multi_t))

covar <- matrix(c(2,1.1,1.1,2),2,2)

eigen_covar <- eigen(covar)

covt <- diag(c(1,1,3),3,3)

eigen(covt)

########
library(fMultivar)

multi_st <- function(t,u){
  fMultivar::dmvst(cbind(t,u),Omega=matrix(c(2,1.1,1.1,2),2,2), alpha = c(-1,0.25))
}

rang <- seq(-2,2,0.04)

contour(x=rang,y=rang,z=outer(rang,rang,multi_st))

library(sn)

fit_skew <- mstFit(CRSPday[,4:7])

fit_skew@fit$dp$beta

fit_skew@fit$dp$alpha

library(nlme)

boot <- 200

alpha_matrix <- matrix(rep(0,nrow(CRSPday)*4),nrow(CRSPday),4)

for(i in 1:boot){
  
  index <- sample(1:nrow(CRSPday),size = nrow(CRSPday), replace = T) ## bootstapping index
  
  new_data <- CRSPday[index,4:7]
  
  fit_skew <- mstFit(new_data)
  
  alpha_matrix[i,] <- fit_skew@fit$dp$alpha
  
  print(alpha_matrix[i,])
  
  print("Bootrapping the shit out stocks...")

}

plot(density(alpha_matrix[,1]))  # all are centered around 0 which means little skew

plot(density(alpha_matrix[,2]))

plot(density(alpha_matrix[,3]))

plot(density(alpha_matrix[,4]))

########################################

berndtInvest <- read.csv("berndtInvest.csv")

datalab <- berndtInvest[,c(1:5,18)]

print(cov(datalab[,2:5]))

pairs(datalab[,2:5])

weights <- as.matrix(c(0.5,0.3,0.2))

covar_3_assets <- cov(datalab[,2:4])

cov_port <- t(weights) %*% covar_3_assets %*% weights

###############

library(MASS)
library(mnormt)

df <- seq(2.5,8,.01)

n <- length(df)

loglik_max <- rep(0,n)

for(i in 1:n){
  fit <- cov.trob(datalab[,2:5],nu=df[i])
  mu <- as.vector(fit$center)
  sigma <- matrix(fit$cov,nrow=4)
  loglik_max[i] <- sum(log(dmt(datalab[,2:5],mean=fit$center,S=fit$cov,df=df[i])))
}

print(df[which.max(loglik_max)])

plot(y = loglik_max, x = df)

fit_skew <- mstFit(x=as.matrix(datalab[,2:5]))

lik_ci <- max(loglik_max) - 1/2*qchisq(.9,1)  ## 90% profile likelihood interval

index <- which(loglik_max > lik_ci)

rangedf <- df[index]

######################################
par(mfrow=c(1,4))
N = 2500
nu = 3

set.seed(5640)
cov=matrix(c(1,.8,.8,1),nrow=2)
x= mvrnorm(N, mu = c(0,0), Sigma=cov)
w = sqrt(nu/rchisq(N, df=nu))
x1 = x * cbind(w,w)
plot(x1,main="(a)")

set.seed(5640)
cov=matrix(c(1,.8,.8,1),nrow=2)
x= mvrnorm(N, mu = c(0,0), Sigma=cov)
w1 = sqrt(nu/rchisq(N, df=nu))
w2 = sqrt(nu/rchisq(N, df=nu))
x2 = x * cbind(w1,w2)
plot(x2,main="(b)")

set.seed(5640)
cov=matrix(c(1,0,0,1),nrow=2)
x= mvrnorm(N, mu = c(0,0), Sigma=cov)
w1 = sqrt(nu/rchisq(N, df=nu))
w2 = sqrt(nu/rchisq(N, df=nu))
x3 = x * cbind(w1,w2)
plot(x3,main="(c)")

set.seed(5640)
cov=matrix(c(1,0,0,1),nrow=2)
x= mvrnorm(N, mu = c(0,0), Sigma=cov)
w = sqrt(nu/rchisq(N, df=nu))
x4 = x * cbind(w,w)
plot(x4,main="(d)")

############

mtfit <- function(series){
  
  loglik <- function(par) {
    
    mean <- par[1:4]
    A <- matrix(c(par[5],par[6],par[7],par[8],0,par[9],par[10],par[11],0,0,par[12],
                  par[13],0,0,0,par[14]),nrow=4,byrow=T)
    
    covar <- t(A)%*%A
    
    nu <- par[15]
    
    f <- -sum(log(dmvt(x=series, delta=mean, sigma=covar, df=nu, log=FALSE)))
    f
  }
  
  A <- chol(cov(series))
  
  cov1 <- c(-.1,-.1,-.1,-.1,-.1,-.1,-.1,-.1,-.1,-.1)
  cov2 <- c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1)
  
  lower <- append(append(c(-0.02,-0.02,-0.02,-0.02),cov1),3)
  upper <- append(append(c(0.02,0.02,0.02,0.02),cov2),15)
  
  start <- as.vector(c(apply(series,2,mean),A[1,1],A[1,2],A[1,3],A[1,4],A[2,2],
                       A[2,3],A[2,4],A[3,3],A[3,4],A[4,4],4))
  
  results <- optim(start, loglik, method = "L-BFGS-B",lower=lower,upper=upper,hessian=T)
  return(results)
}

fit_mine <- mtfit(CRSPday[,c(4:7)])

hessian_MLE <- fit_mine$hessian

fisher_info <- -hessian_MLE

print(dim(hessian_MLE))

SE_par <- solve(fisher_info)

par_mle <- fit_mine$par

mean_mle <- par_mle[1:4]

MLE_A <- matrix(c(par_mle[5],par_mle[6],par_mle[7],par_mle[8],0,par_mle[9],par_mle[10],par_mle[11],
                  0,0,par_mle[12],par_mle[13],0,0,0,par_mle[14]),nrow=4,byrow=T)

covariance_MLE <- t(MLE_A) %*% MLE_A

rho_MLE <- cov2cor(covariance_MLE)
  
##########################
library(zoo)

library(readxl)

data <- readxl::read_xlsx("sp_monthly_1.xlsx")

data <- data %>% mutate(lagprice = lag(price), return = (price-lagprice)/lagprice)

data <- data[-1,]

cumreturn <- function(vec){
  val <- prod(1+vec) - 1
  return(val)
}

cumulative_return <- rollapply(data$return,12,cumreturn)

cum_ret_2010_fwd <- rollapply(data$return[272:375],12,cumreturn)

mth_returns <- data$return[283:375]

plot(density(mth_returns))
plot(density(cum_ret_2010_fwd))

qqnorm(mth_returns,datax = T)
qqline(mth_returns,datax = T)

qqnorm(data$return,datax = T)
qqline(data$return,datax = T)

qqnorm(data$return[200:375],datax = T)
qqline(data$return[200:375],datax = T)

qqnorm(cumulative_return,datax = T)
qqline(cumulative_return,datax = T)

data_2 <- readxl::read_xlsx("sp_daily.xlsx")

data_2 <- data_2 %>% mutate(lag_price = lag(price), return = (price - lag_price)/lag_price)

qqnorm(data_2$return,datax = T)
qqline(data_2$return,datax = T)

qqnorm(data_2$price,datax = T)
qqline(data_2$price,datax = T)

plot(density(data_2$price))

plot(density(data$price))

qqnorm(data$price,datax = T)
qqline(data$price,datax = T)

############################

library(mvtnorm)
df <- 5
mean <- c(0.001,0.002)
covar <- matrix(c(0.1,0.03,0.03,0.15),nrow=2,byrow = T)
weight <- as.matrix(c(0.5,0.5))

var_port <- t(weight) %*% covar %*% weight

mean_port <- 0.5*(sum(mean))

## port will be t distributed, don't know the proof

n <- 10000

var_change <- sqrt(var_port * (df-2)/df)[1]  # student t distribtion has a variance of df/(df-2)

return_port <- rt(n, df=df)*var_change + mean_port  # we use sqrt as variance is weight

quant_99 <- quantile(return_port,0.99)

return_above <- return_port[which(return_port > quant_99)]

print(mean(return_above))

################################### Exercises

wt_mat <- matrix(c(0.2,0.8),2,1)

mean_port <- 0.2*1 + 0.8*1.5

cov_mat <- matrix(c(2,0.8,0.8,2.7),2,2,byrow = T)

var_port <- t(wt_mat) %*% cov_mat %*% wt_mat

wt_optim_var <- function(cov_mat){
  var_def <- function(w){
    weight_mat <- matrix(c(w,1-w),2,1)
    var_port <- t(weight_mat) %*% cov_mat %*% weight_mat
    return(var_port)
  }
  results <- optim(0,var_def,lower=0,upper=1,method ="L-BFGS",hessian=T,control=list(trace=T))
  return(results)
}

optim_par_var <- wt_optim_var(cov_mat)

wt_optim_var <- optim_par_var$par

wt_optim_exp <- function(exp_asset){
  exp_def <- function(w){
    weight_mat <- c(w,1-w)
    exp_port <- -sum(weight_mat * exp_asset)
    return(exp_port)
  }
  results <- optim(0,exp_def,lower=0,upper=1,method ="L-BFGS",hessian=T,control=list(trace=T))
  return(results)
}

optim_par_exp <- wt_optim_exp(c(1,1.5))

wt_optim_exp <- optim_par_exp$par

