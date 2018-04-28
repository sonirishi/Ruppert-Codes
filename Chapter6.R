rm(list=ls(all=T))
setwd("E:/Documents/Practice")

library(dplyr)
library(Ecdat)
library(data.table)
library(fGarch)
library(MASS)

data("CRSPday")

B <- 1000

fit_t_dist_orig <- fitdistr(CRSPday[,4],densfun = "t")

## We can use sqrt of covariance matrix to get the variance 11,22,33

mlist <- list()
sdlist <- list()
dflist <- list()

set.seed(2018)

for(i in 1:B){
  y <- sample(CRSPday[,4],length(CRSPday[,4]),replace=T)  # bootstrap sample
  fit_t_dist <- fitdistr(y,densfun = "t")
  mlist[i] <- fit_t_dist$estimate[1]
  sdlist[i] <- fit_t_dist$estimate[2]
  dflist[i] <- fit_t_dist$estimate[3]
}

print(mean(unlist(mlist))); print(sd(unlist(mlist)))
      
print(mean(unlist(sdlist))); print(sd(unlist(sdlist)))

print(mean(unlist(dflist))); print(sd(unlist(dflist)))

##############################

B <- 1000

fit_t_dist_orig <- fitdistr(CRSPday[1:250,4],densfun = "t")

## We can use sqrt of covariance matrix to get the variance 11,22,33

mlist <- list()
sdlist <- list()
dflist <- list()

set.seed(2018)

for(i in 1:B){
  y <- sample(CRSPday[1:250,4],length(CRSPday[1:250,4]),replace=T)  # bootstrap sample
  fit_t_dist <- fitdistr(y,densfun = "t")
  mlist[i] <- fit_t_dist$estimate[1]
  sdlist[i] <- fit_t_dist$estimate[2]
  dflist[i] <- fit_t_dist$estimate[3]
}

print(mean(unlist(mlist))); print(sd(unlist(mlist)))

print(mean(unlist(sdlist))); print(sd(unlist(sdlist)))

print(mean(unlist(dflist))); print(sd(unlist(dflist)))

plot(density(unlist(dflist)))

######## Test CLT on Dice #####

B <- 5000

roll_list <- list()
mean_roll_list <- list()

for(i in 1:B){
  y <- sample(1:6,B,replace=T)
  roll_list <- list(roll_list,y)
  mean_roll_list[i] <- mean(unlist(y))
  gc()
}

roll_list[[1]] <- NULL

plot(density(unlist(mean_roll_list)))

qqnorm(unlist(mean_roll_list),datax = T)
qqline(unlist(mean_roll_list),datax = T)  ## CLT proved for dice

#####################################

set.seed(5678)
B <- matrix(0,1000,2)
for(i in 1:1000)
{
  x = rnorm(4000) 
  y = 1 + 2*x + rt(4000,2.01)
  g = lm(y~x)
  B[i,] = coef(g)  #### coefficient of x
}
qqnorm(B[,2])
qqline(B[,2])

##################

quKurt <- function(y,p1=0.025,p2=0.25){
  Q <- quantile(y,c(p1,p2,1-p2,1-p1))
  return(as.numeric(((Q[4]-Q[1])/(Q[3]-Q[2]))))
}

library(bootstrap)

bmwRet <- read.csv("bmwRet.csv")

set.seed(1988)

kurtosis <- quKurt(bmwRet[,3])

bcakurt <- bcanon(bmwRet[,3],5000,quKurt)

print(bcakurt$confpoints)

bcakurt_2 <- bcanon(bmwRet[,3],5000,moments::kurtosis)

print(bcakurt_2$confpoints)

qkurtosis <- matrix(rep(0,5000),5000,1)

for(i in 1:5000){
  yboot <- sample(bmwRet[,3],nrow(bmwRet),replace = T)
  qkurtosis[i] <- quKurt(yboot)
}

confikurt <- quantile(qkurtosis,c(0.025,0.05,0.1,0.16,0.84,0.9,0.95,0.975))  # this is simple percentile method

plot(density(qkurtosis))  ## decently symmetric

print(2*quKurt(bmwRet[,3])-confikurt[2])  # THis is the basic bootstrap interval
print(2*quKurt(bmwRet[,3])-confikurt[7])

library(zoo)

rollmeankurt <- rollapply(qkurtosis,50,mean)
rolldskurt <- rollapply(qkurtosis,50,sd)

plot(rollmeankurt,type="l")

plot(rolldskurt,type="l")

meanrolldskurt <- rollapply(rolldskurt,50,mean)

plot(rolldskurt,type="l")
lines(meanrolldskurt,type="l",col="yellow")

midcap <- read.csv("midcapD.csv")

qqplot(midcap$LSCC,midcap$CSGS)

comparequkurt <- function(x,p1=0.025,p2=0.25,xdata){
  quKurt(xdata[x,1],p1,p1)/quKurt(xdata[x,2],p1,p2)
}

quKurt(midcap$LSCC)

quKurt(midcap$CSGS)

xdata <- cbind(midcap$LSCC,midcap$CSGS)

comparequkurt(1:nrow(midcap),xdata=xdata)

bca_midcap <- bcanon((1:nrow(midcap)),5000,comparequkurt,xdata=xdata)

print(bca_midcap$confpoints)                     

####################

library(fGarch)

kurt <- kurtosis(bmwRet[,3])
skew <- skewness(bmwRet[,3])

fit_skewt <- sstdFit(bmwRet[,3])

q.grid <- (1:nrow(bmwRet))/(nrow(bmwRet)+1)

qqplot(bmwRet[,3],qsstd(q.grid,fit_skewt$estimate[1],fit_skewt$estimate[2],
                        fit_skewt$estimate[3],fit_skewt$estimate[4]))


######################

nboot <- 5000

ModelFree_kurt <- rep(0,nboot)

ModelBased_kurt <- rep(0,nboot)

set.seed(1988)

for(i in 1:nboot){
  samp_modelfree <- sample(bmwRet[,3],nrow(bmwRet),replace=T)
  samp_modelbased <- rsstd(nrow(bmwRet),fit_skewt$estimate[1],fit_skewt$estimate[2],
                          fit_skewt$estimate[3],fit_skewt$estimate[4])
  ModelFree_kurt[i] <- quKurt(samp_modelfree)
  ModelBased_kurt[i] <- quKurt(samp_modelbased)
}

plot(density(ModelBased_kurt),type="l")
lines(density(ModelFree_kurt),col="red",type="l")

mf_low <- quantile(ModelFree_kurt,0.05)
mf_high <- quantile(ModelFree_kurt,0.95)

mb_low <- quantile(ModelBased_kurt,0.05)
mb_high <- quantile(ModelBased_kurt,0.95)

print(paste0("Model Based Bootstrap: ",2*quKurt(bmwRet[,3])-mb_high," - ",2*quKurt(bmwRet[,3])-mb_low))

print(paste0("Model Free Bootstrap: ",2*quKurt(bmwRet[,3])-mf_high," - ",2*quKurt(bmwRet[,3])-mf_low))

bootinterval <- bcanon(bmwRet[,3],5000,quKurt)

bootinterval$confpoints

###########################

sboot5bys <- 0.71
sboot95bys <- 1.67
s <- 0.31

print(paste0("90percentile CI of SD: ",2*s-s*sboot95bys," - ",2*s-s*sboot5bys))

################

bias <- 0.68431-0.69119
sd <- .11293

MSE <- bias^2 + sd^2

print(MSE)

bias^2/sd^2
