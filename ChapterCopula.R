rm(list=ls(all=T))

setwd("E:/Documents/Practice")

library(dplyr)
library(data.table)

norm_data <- rnorm(100,0,1)

CDF_y <- pnorm(norm_data)

plot(CDF_y)  # This proves that CDF of y is uniform distribution

hist(CDF_y)

library(copula)

fc_obj <- frankCopula(param=500)

fc_random <- rCopula(200,fc_obj)

plot(fc_random)

library(infotheo)

data(USArrests)

dat<-discretize(USArrests)

mutinformation(USArrests[,1],USArrests[,2])

cor(USArrests[,1],USArrests[,2],method = "spearman")

#####################

gasflow <- read.csv("FlowData.csv")

library(ggplot2)
library(fGarch)

ggplot(data=gasflow) + aes(x=Flow1) + geom_histogram()

ggplot(data=gasflow) + aes(x=Flow2) + geom_histogram()

ggplot(data=gasflow) + geom_point(aes(x=Flow1,y=Flow2),color="blue")

skew_t_flow1 <- fGarch::sstdFit(gasflow$Flow1)

skew_t_flow2 <- fGarch::sstdFit(gasflow$Flow2)

cdf_transformed_flow1 <- psstd(gasflow$Flow1,skew_t_flow1$estimate[1],
                               skew_t_flow1$estimate[2],skew_t_flow1$estimate[3],
                               skew_t_flow1$estimate[4])  ### Estimate marginal and then use CDF

cdf_transformed_flow2 <- psstd(gasflow$Flow2,skew_t_flow2$estimate[1],
                               skew_t_flow2$estimate[2],skew_t_flow2$estimate[3],
                               skew_t_flow2$estimate[4])  ### Estimate marginal and then use CDF

hist(cdf_transformed_flow1,col="blue")

hist(cdf_transformed_flow2,col="green")

plot(cdf_transformed_flow1,cdf_transformed_flow2,col="orange")  ## Scatter Plot of this data

cor(gasflow$Flow1,gasflow$Flow2,method = "spearman")

val <- cor(gasflow$Flow1,gasflow$Flow2,method="kendall")

print(sin(pi*val/2))

p_corr <- cor(gasflow$Flow1,gasflow$Flow2)

boot <- 10000
pearson_corr <- matrix(rep(0,1000),1000,1)

###### Right way to do the correlation bootstrap

for(i in 1:boot){
  index <- sample(1:nrow(gasflow),nrow(gasflow),replace=T)  ## boostrap same index values
  sample1 <- gasflow$Flow1[index]  # Get both data for same index values not different
  sample2 <- gasflow$Flow2[index]  ## Getting different values changes correlation which makes sense
  pearson_corr[i] <- cor(sample1,sample2)
}

conf_corr <- quantile(pearson_corr,c(0.025,0.975))

print(2*p_corr - conf_corr[1])  ## dont forget this boostrap interval
print(2*p_corr - conf_corr[2])   # Can use bcanon function as well

cor.test(gasflow$Flow1,gasflow$Flow2)

######### Fir copula to this data

cop_t_dim3 <- tCopula(c(-0.6,0.75,0),dim = 3, dispstr = "un",df = 1)

set.seed(5640)

rand_t_Cop <- copula::rCopula(500,cop_t_dim3)

pairs(rand_t_Cop)

corrs <- cor(rand_t_Cop)  

cor(rand_t_Cop,method = "kendall")

cor.test(rand_t_Cop[,1], rand_t_Cop[,2])

taildep_tCop <- function(mu,rho){
  val <- -sqrt((mu+1)*(1-rho)/(1+rho))
  return(val)
}

tail_dep_val_1_2 <- 2*pt(taildep_tCop(2,-0.54999514),2)  ### correlation between tail

tail_dep_val_1_3 <- 2*pt(taildep_tCop(2,0.70707296),2)  ## correlation values from corrs

tail_dep_val_2_3 <- 2*pt(taildep_tCop(2,-0.06538499),2)  ## Seems like tail dependence


cop_n_dim3 <- normalCopula(c(-0.6,0.75,0),dim = 3, dispstr = "un")

mvdc_normal <- mvdc(cop_n_dim3,c("exp","exp","exp"),list(list(rate=2),list(rate=3),list(rate=4)))

set.seed(5640)

rand_mvdc <- rMvdc(1000,mvdc_normal)

pairs(rand_mvdc)

plot(density(rand_mvdc[,1]))

plot(density(rand_mvdc[,2]))

plot(density(rand_mvdc[,3]))

############################################

library(Ecdat)
library(copula)
library(MASS)
library(fGarch)
library(fCopulae)

data("CRSPday")

ibm <- CRSPday[,5]
crsp <- CRSPday[,7]

est.ibm <- as.numeric(fitdistr(ibm,"t")$estimate)
est.crsp <- as.numeric(fitdistr(crsp,"t")$estimate)

est.ibm[2] <- est.ibm[2]*sqrt(est.ibm[3]/(est.ibm[3]-2)) ## standard deviation

est.crsp[2] <- est.crsp[2]*sqrt(est.crsp[3]/(est.crsp[3]-2))  ## Convert to classical T Dist.

cor_tau <- cor(ibm,crsp,method = "kendall")

omega <-  sin(pi*cor_tau/2)

cop_t_dim2 <- tCopula(omega,dim=2,dispstr = "un",df=4)  ## copula depends on only correlation

## We don't want it to depend on the marginal hence kendall tau or spearman

n <- length(ibm)

data1 <- cbind(pstd(ibm,mean=est.ibm[1],sd=est.ibm[2],nu=est.ibm[3]),
               pstd(crsp,mean=est.crsp[1],sd=est.crsp[2],nu=est.crsp[3]))

data2 <- cbind(rank(ibm)/(n+1),rank(crsp)/(n+1))

#### two parameters: Correlation and Degrees of Freedom

ft1 <- fitCopula(cop_t_dim2,method="mpl",data=data1,start=c(omega,5),lower=c(0,2.5),
                 upper=c(0.5,15),optim.method="L-BFGS-B")

ft2 <- fitCopula(cop_t_dim2,method="mpl",data=data2,start=c(omega,5),lower=c(0,2.5),
                 upper=c(0.5,15),optim.method="L-BFGS-B")

print(ft1@estimate)

mvdc_t_t <- mvdc(cop_t_dim2,c("std","std"),
                 list(list(mean=est.ibm[1],sd=est.ibm[2],nu=est.ibm[3]),
                 list(mean=est.crsp[1],sd=est.crsp[2],nu=est.crsp[3])))

print(ft2@estimate)

start <- c(est.ibm,est.crsp,ft1@estimate)

objFn <- function(param){
  -loglikMvdc(param,cbind(ibm,crsp),mvdc_t_t)
}

fit_cop <- optim(start,objFn,method="L-BFGS-B",
                 lower = start,upper = c(0.1,0.03,15,0.1,0.03,15,8,15))

taildep_tCop(ft1@estimate[2],ft1@estimate[1])

2*pt(taildep_tCop(ft1@estimate[2]+1,-0.54999514),ft1@estimate[2]+1)  # No Tail dependence

###################

fnorm <- fitCopula(data=data1,copula = normalCopula(-0.3,dim = 2),method="ml",start=0.5,
                   optim.method="BFGS")

fgumbel <- fitCopula(data=data1,copula = gumbelCopula(1,dim = 2),method="ml",start=1,
                   optim.method="BFGS")

ffrank <- fitCopula(data=data1,copula = frankCopula(3,dim = 2),method="ml",start=1,
                     optim.method="BFGS")

fclayton <- fitCopula(data=data1,copula = claytonCopula(1,dim = 2),method="ml",start=1,
                    optim.method="BFGS")

#######################

u1 <- data1[,1]
u2 <- data1[,2]

dem <- pempiricalCopula(u1,u2)  ## empirical copula probability, it uses heaviside function

contour(dem$x,dem$y,dem$z)

contour(kde2d(u1,u2))

rn1 <- runif(100,-10,0)
rn2 <- runif(100,1,10)

cor(rn1,rn2,method = "kendall")  # Negative of actual kendall tau

cor(1/rn1,1/rn2,method = "kendall")  ## same as original data
