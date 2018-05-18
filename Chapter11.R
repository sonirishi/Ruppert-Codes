rm(list=ls(all=T))

setwd("E:/Documents/Practice")

library(dplyr)
library(data.table)

capital <- 1000000

loss_cover <- 150000

per_loss <- loss_cover/capital

w <- -0.21/(0.25*qnorm(0.01) + 0.09)  ## This is basically to get to > -15% returns

mu1 <- 0.14
mu2 <- 0.08
sigma1 <- 0.2
sigma2 <- 0.15
rho12 <- 0

portfolio <- function(mu,sigma,rho){
  
  ######### Variance not given, convert by squaring the SD
  
  var_matrix <- matrix(c(sigma[1]^2,rho*sigma[1]*sigma[2],rho*sigma[1]*sigma[2],sigma[2]^2),2,2)
  
  risk_return <- function(par){
    
    w <- matrix(c(par,1-par),2,1)
  
    returns <- w %*% mu
  
    risk <- t(w) %*% var_matrix %*% w
    
    return(risk)
  }
  
  #### One dimensional optimization by Brent is better, Nelder mead not reliable
  
  value <- optim(0,risk_return,method = "L-BFGS",lower=0,upper=1)  
  
  return(value)
}

optim_par <- portfolio(mu=c(mu1,mu2),sigma=c(sigma1,sigma2),rho=rho12)

tangency_portfolio <- function(mu1,mu2,muf,sigma1,sigma2,rho12){
  num <- (mu1-muf)*sigma2^2 - (mu2-muf)*rho12*sigma1*sigma2
  denom <- (mu1-muf)*sigma2^2 + (mu2-muf)*sigma1^2 - (mu1+mu2-2*muf)*rho12*sigma1*sigma2
  return(num/denom)
}

w1 <- tangency_portfolio(0.14,0.08,0.06,0.2,0.15,0)

port_return <- w1*0.14 + (1-w1)*0.08
std <- sqrt(w1^2*0.2^2 + (1-w1)^2*0.15^2)

new_std <- 0.05

weight_risky <- new_std/std

new_return <- 0.06*(1-weight_risky) + port_return*weight_risky

target_return <- 0.1

new_wt <- uniroot(function(x){x*0.14 + (1-x)*0.08 - 0.1},lower=0,upper=1)$root

new_sd <- sqrt(new_wt^2*0.2^2 + (1-new_wt)^2*0.15^2)

new_wt1 <- uniroot(function(x){x*port_return + (1-x)*0.06 - 0.1},lower=0,upper=1)$root

std_fin <- new_wt1*std

print((new_sd-std_fin)/std_fin)  # reduction in risk by 27%

rho12 <- c(0.7,0.3,0,-0.7)

weight_range <- seq(0,1,0.01)

expected_return <- list()
expected_std <- list()

for(i in 1:length(weight_range)){
  expected_return[[i]] <- weight_range[i]*0.14 + (1-weight_range[i])*0.09
  expected_std[[i]] <- sqrt(weight_range[i]^2*0.2^2 + (1-weight_range[i])^2*0.15^2)
}

plot(y=unlist(expected_return),x=unlist(expected_std),type='l',col="green")

expected_return <- list()
expected_std <- list()

for(i in 1:length(weight_range)){
  expected_return[[i]] <- weight_range[i]*0.14 + (1-weight_range[i])*0.09
  expected_std[[i]] <- sqrt(weight_range[i]^2*0.2^2 + (1-weight_range[i])^2*0.15^2 + 
                              2*weight_range[i]*(1-weight_range[i])*0.7*0.2*0.15)
}

plot(y=unlist(expected_return),x=unlist(expected_std),type='l',col="orange")

expected_return <- list()
expected_std <- list()

for(i in 1:length(weight_range)){
  expected_return[[i]] <- weight_range[i]*0.14 + (1-weight_range[i])*0.09
  expected_std[[i]] <- sqrt(weight_range[i]^2*0.2^2 + (1-weight_range[i])^2*0.15^2 + 
                              2*weight_range[i]*(1-weight_range[i])*(-0.7)*0.2*0.15)
}

plot(y=unlist(expected_return),x=unlist(expected_std),type='l',col="orange")

################3

library(quadprog)
library(Ecdat)

data("CRSPday")

R <- 100*CRSPday[,4:6]

########## inequality contraints are converted to equality contraints

mean_vect <- as.vector(apply(R,2,mean))  ## mean of individual stock returns

cov_mat <- cov(R)  # covariance between assets

sd_vect <- sqrt(diag(cov_mat))  # get std of stocks

Amat <- cbind(rep(1,3),mean_vect)

muP <- seq(0.05,0.14,length=300)  ## range of target return of portfolio

sdP <- muP  # setup storage for std  smart way

weights <- matrix(0,nrow=300,ncol=3)

#### Minimum covariance portfolio with different means

for(i in 1:length(muP)){
  bvec <- c(1,muP[i])
  result <- solve.QP(Dmat = 2*cov_mat,dvec = rep(0,3),Amat = Amat, bvec=bvec,meq=2)
  ### 2*covariance matrix because of the way the code needs inputs
  sdP[i] = sqrt(result$value)
  weights[i,] <- result$solution
}

mufree <- 1.3/252 ## no of trading days

sharpe <- (muP-mufree)/sdP

ind <- which(sharpe == max(sharpe))

weights[ind,]  # tangency portfolio weights

ind2 <- sdP == min(sdP)

sdP[ind2]

ind3 <- (muP > muP[ind2])

############### Put no shorting constraints

Amat <- as.matrix(cbind(rep(1,3),mean_vect,diag(-1,nrow=3),diag(1,nrow=3)))

muP <- seq(0.08,0.1,length=300)  ## range of target return of portfolio, keeping less than max of stocks

sdP <- muP  # setup storage for std  smart way

weights <- matrix(0,nrow=300,ncol=3)

#### Minimum covariance portfolio with different means

for(i in 1:length(muP)){
  bvec <- c(1,muP[i],rep(-1,3),rep(0,3))
  result <- solve.QP(Dmat = 2*cov_mat,dvec = rep(0,3),Amat = Amat, bvec=bvec,meq=2)
  ### 2*covariance matrix because of the way the code needs inputs
  sdP[i] = sqrt(result$value)
  weights[i,] <- result$solution
}

### QC the constraints

min(weights[,1]); max(weights[,1])

min(weights[,2]); max(weights[,2])

min(weights[,3]); max(weights[,3])

apply(weights,1,sum)  # weights are 1 which is required

##############################################################

ret_data <- read.table("countries.txt",header = T)

R = 100*(ret_data[2:nrow(ret_data),]/ret_data[1:(nrow(ret_data)-1),] - 1)

mean_returns <- as.matrix(apply(R[,4:13],2,mean))

cov_mat <- cov(R[,4:13])

muP <- seq(0,2.5,length=300)
sdP <- muP

weights <- matrix(rep(0,length(muP)*10),length(muP),10)

Amat <- cbind(rep(1,10),mean_returns)

for(i in 1:length(muP)){
  bvec <- c(1,muP[i])
  result <- solve.QP(Dmat = 2*cov_mat,dvec = rep(0,10),Amat = Amat, bvec=bvec,meq=2)
  ### 2*covariance matrix because of the way the code needs inputs
  sdP[i] = sqrt(result$value)
  weights[i,] <- result$solution
}

mufree = 1/24 

sharpe_ratio_act <- (muP - mufree)/sdP

ind_max_sharpe <- which(sharpe_ratio_act == max(sharpe_ratio_act))

tang_sharpe <- sharpe_ratio_act[ind_max_sharpe]

print(sharpe_ratio_act[ind_max_sharpe])  ### max sharpe ratio

############ Boostrap Model ####
boot <- 250

estimated_sharpe <- matrix(rep(0,boot),boot,1)

actual_sharpe <- matrix(rep(0,boot),boot,1)

for(j in 1:boot){
  
  index_boot <- sample(1:nrow(R),nrow(R),replace = T)
  
  R_boot <- R[index_boot,]

  mean_return_boot <- as.matrix(apply(R_boot[,4:13],2,mean))

  cov_mat_boot <- cov(R_boot[,4:13])

  muP <- seq(0,2.5,length=300)
  sdP <- muP

  weights <- matrix(rep(0,length(muP)*10),length(muP),10)

  Amat <- cbind(rep(1,10),mean_return_boot)

  for(i in 1:length(muP)){
    bvec <- c(1,muP[i])
    result <- solve.QP(Dmat = 2*cov_mat_boot,dvec = rep(0,10),Amat = Amat, bvec=bvec,meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] <- result$solution
  }

  mufree = 1/24 

  sharpe_ratio <- (muP - mufree)/sdP

  ind_max_sharpe <- which(sharpe_ratio == max(sharpe_ratio))
  
  estimated_sharpe[j] <- sharpe_ratio[ind_max_sharpe]  ## This is the one using boot

  ########## This is the sharpe ratio using the sample mean/covar which is the actual
  
  actual_sharpe[j] <- (weights[ind_max_sharpe,] %*% mean_returns - mufree)/sqrt(weights[ind_max_sharpe,] %*% cov_mat %*% weights[ind_max_sharpe,])
}

final_data <- cbind(rbind(estimated_sharpe,actual_sharpe),
                    rbind(matrix(rep("estimated",250)),matrix(rep("actual",250))))

final_data <- as.data.frame(final_data)

colnames(final_data) <- c("sharpe","type")

final_data$sharpe <- as.numeric(as.character(final_data$sharpe))

boxplot(sharpe~type,data=final_data)
abline(h=tang_sharpe)

###############################

ret_data <- read.table("countries.txt",header = T)

R = 100*(ret_data[2:nrow(ret_data),]/ret_data[1:(nrow(ret_data)-1),] - 1)

mean_returns <- as.matrix(apply(R[,4:13],2,mean))

cov_mat <- cov(R[,4:13])

boot <- 250

estimated_sharpe <- matrix(rep(0,boot),boot,1)

actual_sharpe <- matrix(rep(0,boot),boot,1)

for(j in 1:boot){
  
  index_boot <- sample(1:nrow(R),nrow(R),replace = T)
  
  R_boot <- R[index_boot,]
  
  mean_return_boot <- as.matrix(apply(R_boot[,4:13],2,mean))
  
  total_mean_boot <- mean(mean_return_boot)
  
  mean_return_boot <- (mean_return_boot + total_mean_boot)/2   # with shrinkage parameter 1/2
  
  cov_mat_boot <- cov(R_boot[,4:13])
  
  muP <- seq(min(mean_returns) + 0.001,max(mean_returns) - 0.001,length=300)
  sdP <- muP
  
  weights <- matrix(rep(0,length(muP)*10),length(muP),10)
  
  Amat <- cbind(rep(1,10),mean_return_boot,diag(rep(1,10)))
  
  for(i in 1:length(muP)){
    bvec <- c(1,muP[i],rep(-1,10))
    result <- solve.QP(Dmat = 2*cov_mat_boot,dvec = rep(0,10),Amat = Amat, bvec=bvec,meq=2)
    sdP[i] = sqrt(result$value)
    weights[i,] <- result$solution
  }
  
  mufree = 1/24 
  
  sharpe_ratio <- (muP - mufree)/sdP
  
  ind_max_sharpe <- which(sharpe_ratio == max(sharpe_ratio))
  
  estimated_sharpe[j] <- sharpe_ratio[ind_max_sharpe]  ## This is the one using boot
  
  ########## This is the sharpe ratio using the sample mean/covar which is the actual
  
  actual_sharpe[j] <- (weights[ind_max_sharpe,] %*% mean_returns - mufree)/sqrt(weights[ind_max_sharpe,] %*% cov_mat %*% weights[ind_max_sharpe,])
}

final_data <- cbind(rbind(estimated_sharpe,actual_sharpe),
                    rbind(matrix(rep("estimated",250)),matrix(rep("actual",250))))

final_data <- as.data.frame(final_data)

colnames(final_data) <- c("sharpe","type")

final_data$sharpe <- as.numeric(as.character(final_data$sharpe))

boxplot(sharpe~type,data=final_data)

##############################################

dat <- read.csv("Stock_FX_Bond.csv")

data <- dat %>% select(ends_with("AC"))  ## not to be used like this

prices <- data[,c("GM_AC","F_AC","CAT_AC","UTX_AC","MRK_AC","IBM_AC")]

R = 100*(prices[2:nrow(prices),]/prices[1:(nrow(prices)-1),] - 1)

mean_returns <- as.matrix(apply(R,2,mean))

cov_mat <- cov(R)

muP <- seq(0.037,0.073,length=300)
sdP <- muP

weights <- matrix(rep(0,length(muP)*6),length(muP),6)

Amat <- cbind(rep(1,6),mean_returns,diag(rep(1,6)),diag(rep(-1,6)))

for(i in 1:length(muP)){
  bvec <- c(1,muP[i],rep(-0.1,6),rep(-0.5,6))
  result <- solve.QP(Dmat = 2*cov_mat,dvec = rep(0,6),Amat = Amat, bvec=bvec,meq=2)
  ### 2*covariance matrix because of the way the package needs inputs
  sdP[i] = sqrt(result$value)
  weights[i,] <- result$solution
}

mufree = 3/365 ## interest earned for non trading days as well so divide by 365 for daily value

sharpe_ratio_act <- (muP - mufree)/sdP

ind_max_sharpe <- which(sharpe_ratio_act == max(sharpe_ratio_act))

tang_sharpe <- sharpe_ratio_act[ind_max_sharpe]

print(sharpe_ratio_act[ind_max_sharpe])  ### max sharpe ratio

ind_min_var <- which(sdP == min(sdP))

print(sdP[ind_min_var])  ### min var portfolio

###### For 0.07% returns

bvec <- c(1,0.07,rep(-0.1,6),rep(-0.5,6))

result <- solve.QP(Dmat = 2*cov_mat,dvec = rep(0,6),Amat = Amat, bvec=bvec,meq=2)

print(result$solution)  # sums to 1

################## Exercise ################
library(quadprog)

expected_A <- 2.3/100
expected_B <- 4.5/100

std_A <- sqrt(6)/100
std_B <- sqrt(11)/100

rhoAB <- 0.17

expected_port <- 3/100

optimal_wt_mean <- 
  uniroot(function(w){w*expected_A + (1-w)*expected_B - expected_port},lower = 0, upper = 1)$root

std_port <- sqrt(5.5)/100

covar_matrix <- matrix(c(std_A^2,rhoAB*std_A*std_B,rhoAB*std_A*std_B,std_B^2),byrow = T,nrow=2) ## not required

library(rootSolve)

multiroot(function(w){w^2*std_A^2 + (1-w)^2*std_A^2 + 2*std_A*std_B*rhoAB*w*(1-w) - std_port^2},start = 0)

####

weight_tangency <- 0.65

expected_tangency <- 5/100
std_tangency <- 7/100

std_final <- 5/100
mu_free <- 1.5/100

weight_fin <- std_final/std_tangency; weight_rf <- 1-weight_fin

wa <- weight_tangency*weight_fin; wb <- (1-weight_tangency)*weight_fin

expected_final <- weight_fin*expected_tangency + (1-weight_fin)*mu_free

###

price_a <- 100; price_b <- 125

total_a <- price_a*200; total_b <- 100*price_b

weight_a <- total_a/(total_a+total_b)

return_a <- 0.001; return_b <- 0.0015

std_a <- 0.03; std_b <- 0.04

rho <- 0.35

covar <- matrix(c(std_a^2,rho*std_a*std_b,rho*std_a*std_b,std_b^2),byrow = T,nrow=2)

mean_port <- weight_a*return_a + (1-weight_a)*return_b

wt_mat <- matrix(c(weight_a,(1-weight_a)),2,1)

std_port <- sqrt(t(wt_mat) %*% covar %*% (wt_mat))

var_return <- qnorm(0.05,mean,std_port)  ## assumption that returns are normal

### This is the value of returns at probability distribution of returns = 0.05

total_invested <- total_a + total_b

VaR <- var_return*total_invested

print(paste0("First ever VaR at 95% confidence: ",-1*VaR))
