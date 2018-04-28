setwd("E:/Documents/Practice")

library(Ecdat)

data("Mishkin")

library(datasets)

data("AirPassengers")

library(forecast)

Acf(Mishkin[,1])        

Acf(diff(Mishkin[,1]))        

Box.test(Mishkin[,1])

library(evir)

data("bmw")

Acf(bmw)

fit = arima(bmw,order = c(1,0,0))

resid = residuals(fit)

acf(resid)

plot(resid)

qqnorm(resid,datax = T)
qqline(resid,datax = T)
plot(density(resid)$y)

fit_mish = arima(Mishkin[,1],order=c(1,0,1))

fit_mish_resid = residuals(fit_mish); fit_mish_resid = cbind.data.frame(as.data.frame(fit_mish_resid),as.data.frame(Mishkin[,1]))

colnames(fit_mish_resid) = c("resid","act")

plot(fit_mish_resid)

library(ggplot2)

ggplot(data=fit_mish_resid) + geom_point(aes(x=act,y=resid),color="red")

Acf(residuals(fit_mish))

#############################################

library(forecast)

infl_rt = Mishkin[,1]

aic_fit = auto.arima(diff(infl_rt),max.p=15,max.q=0,ic="aic")

bic_fit = auto.arima(diff(infl_rt),max.p=15,max.q=0,ic="bic")

### Inflation Rate AR process

aic_fit = auto.arima(infl_rt,max.p=15,max.q=0,ic="aic")

bic_fit = auto.arima(infl_rt,max.p=15,max.q=0,ic="bic")

#### MA process

aic_fit = auto.arima(infl_rt,max.p=0,max.q=10,ic="aic")

bic_fit = auto.arima(infl_rt,max.p=0,max.q=10,ic="bic")

##################################

library(data.table)
library(ggplot2)

CPI = fread("CPI.csv")

CPI$CPI = as.numeric(CPI$CPI); colnames(CPI)[1] = "date"

CPI[,":="(CPI = log(CPI), diffCPI = diff(CPI), diff2CPI = diff(diff(CPI)))]

ggplot(data=CPI) + geom_point(aes(y=CPI,x=date),color="red") + 
  geom_point(aes(y=diffCPI,x=date),color="black") + geom_point(aes(y=diff2CPI,x=date),color="orange")

Acf(CPI$diff2CPI)

####

ip = fread("IP.csv")

ip$IP = as.numeric(ip$IP)

ip[,":="(diffIP = diff(IP))]; colnames(ip)[1] = "date"

ggplot(data=ip) + geom_point(aes(y=IP,x=1:nrow(ip)),color="red") + 
  geom_point(aes(y=diffIP,x=1:nrow(ip)),color="black") + 
  scale_x_continuous(breaks=c(0,10,50,100,150,500,900))

####### Unit Root

infl_rt = Mishkin[,1]

library(tseries)

adf.test(infl_rt,alternative = "explosive")

arm_infl = forecast::Arima(infl_rt,order=c(2,0,1),include.drift=T)

root = function(x){
  func = 1 - 1.1949*x + 0.2178*x^2
  return(func)
}

polyroot(c(1,-1.1949,0.2178))  ## Root of the polynomial

print(pp.test(infl_rt,alternative = "explosive"))

print(kpss.test(infl_rt))

randwalk = function(x){
  val = 0.05 + 1*x + rnorm(1,0,4)
}

rands = c(0.5)

for(i in 2:1000){
  rands[[i]] = randwalk(rands[[i-1]])
}

list_rands = unlist(rands)

library(zoo)

rollmean = rollapply(list_rands,50,mean)  # Rolling with window = 50, apply function

rollstd = rollapply(list_rands,10,sd)

plot(rollmean,rollstd)

rollmean = list(); rollvar = list()

y = 1

for(i in seq(1,1000,20)){
  rollmean[[y]] = mean(list_rands[seq(i,i+20,1)])
  rollvar[[y]] = var(list_rands[seq(i,i+20,1)])
  y = y + 1
}

plot(unlist(rollmean),unlist(rollvar),type = "l")

plot(unlist(rollmean),type = "l", col="blue")
par(new=TRUE)
plot(list_rands, type="l", col="green")
#par(new=TRUE)
#plot(unlist(rollvar), type="l", col="red")
#par(new=T)
#plot(diff(list_rands), type="l", col="black")
#par(new=T)
#plot(log(list_rands+150), type="l", col="black")

plot(unlist(rollmean),unlist(rollvar))

plot(temp$residuals,type='l')
par(new=T)
plot(diff(diff(x2)),col="red",type='l')

####################

arima_fit = arima(infl_rt,order=c(0,1,3))

ma_fit = arima(diff(infl_rt),order=c(0,1,3))

arima_resid = arima_fit$residuals

ma_resid = ma_fit$residuals

#arima_forecast = forecast(infl_rt,100)

#ma_forecast = forecast(diff(infl_rt),100)

## These 2 models are same but predict different series ###

plot(arima_forecast,col="red")  ### Diverging predictions
par(new=T)
plot(infl_rt,col="blue")

plot(ma_forecast,col="orange")  ### Straight line predictions
par(new=T)
plot(infl_rt,col="blue")

###################

plot(arima_forecast$upper[,1],col="red")
par(new=T)
plot(arima_forecast$lower[,1],col="green")

simulate_infl = list()

for(i in 1:10){
  simulate_infl[[i]] = simulate(arima_fit,30)
}

plot(simulate_infl[[1]],col="red")
lines(simulate_infl[[2]],col="orange")
lines(simulate_infl[[3]],col="blue")
lines(simulate_infl[[4]],col="green")

simulate_infl = matrix(rep(0,10000*30),10000,30)

for(i in 1:10000){
  simulate_infl[i,] = simulate(arima_fit,30)
}

arima_forecast = forecast(arima_fit,30)

###########

library(dplyr)
library(Ecdat)
library(tseries)

pacf(infl_rt)  # slowly dying pacf meaning pure AR is not a good fit

data("Tbrate")

plot(Tbrate)

acf(Tbrate)

adf.test(Tbrate[,1]) # non stationary

adf.test(Tbrate[,2]) # non stationary - Big p value

adf.test(Tbrate[,3]) # non stationary

diff_tb = diff(Tbrate)

acf(diff_tb)

adf.test(diff_tb[,1]) 

adf.test(diff_tb[,2]) 

adf.test(diff_tb[,3]) 

pairs(diff_tb)

rollmean = rollapply(Tbrate[,1],width = 10,mean)
rollsd = rollapply(Tbrate[,1],width = 10,sd)

plot(rollmean,col="red")
lines(rollsd,col="black")

rollmean = rollapply(diff_tb[,1],width = 10,mean)
rollsd = rollapply(diff_tb[,1],width = 10,sd)

plot(rollmean,col="red")
lines(rollsd,col="black")

library(fpp)

data("a10")
plot(a10)

rollmean = rollapply(a10,width = 10,mean)
rollsd = rollapply(a10,width = 10,sd)

plot(rollmean,col="red")
plot(rollsd,col="black")

plot(a10,col="red")
lines(rollmean,col="black")

a11 = log(a10)

rollmean = rollapply(a11,width = 10,mean)
rollsd = rollapply(a11,width = 10,sd)

plot(rollmean,col="red")
plot(rollsd,col="black")

rollmeansq = rollmean^2

model_sd = lm(rollsd ~ rollmean + rollmeansq) ## check relation between mean and SD

diff_series = diff(diff(a10),12)
rollmean = rollapply(diff_series,width = 10,mean)

plot(diff_series,col="orange")
lines(rollmean,col="black")
lines(rollapply(diff_series,width = 10,sd),col="red")   #SD has still some trend, mean is ok

###################
boxplot(diff_tb[,1]~cycle(diff_tb))

auto.arima(Tbrate[,1],max.P = 0,max.Q = 0,ic="aic")

auto.arima(Tbrate[,1],max.P = 0,max.Q = 0,ic="bic")

model = arima(Tbrate[,1],order=c(0,1,1))

acf(model$residuals)

Box.test(model$residuals,lag=10,type="Ljung",fitdf=1) ## null is independence

resid = model$residuals^2

acf(resid)

Box.test(resid,lag=10,type="Ljung",fitdf=1) ## null is independence; autocorrelation present

#########

auto.arima(Tbrate[,3],max.P = 0,max.Q = 0, ic = "aic")

fit = arima(Tbrate[,3],order=c(1,1,1))

forecasts = predict(fit,36)

plot(Tbrate[,3],col="red")

plot(forecasts$pred,col="black")
lines(forecasts$pred + 1.96*forecasts$se,col="blue")
lines(forecasts$pred - 1.96*forecasts$se,col="orange")

#####################################################

library(Ecdat)

data("CRSPday")

crsp = CRSPday[,7]
acf(crsp)

acf(as.numeric(crsp))

fit1 = arima(crsp,order=c(1,0,0))

fit2 = arima(crsp,order=c(2,0,0))

fit1$coef[1] + 1.96*0.0198; fit1$coef[1] - 1.96*0.0198

arima(crsp,c(0,0,1))

tlist = matrix(rep(1,1000),1000,1)

for(i in 1:1000){
  tlist[i+1] = 5 - 0.55*tlist[i] + rnorm(1,0,1.414)
}

acf(tlist)
rollmean = rollapply(tlist,width = 10,mean)

plot(tlist,col="red",type="l")
plot(rollmean,col="black",type="l")

tlist = matrix(rep(0.5,1000),1000,1)

for(i in 1:1000){
  tlist[i+1] = 0.3 + 0.4*tlist[i] + rnorm(1,0,1.414)
}

plot(tlist,col="pink",type='l')

rollmean = rollapply(tlist,width = 10,mean)
rollsd = rollapply(tlist,width = 10,sd)

plot(rollmean,type='l')
plot(rollsd,type='l')

library(Ecdat)
data("Mishkin")

tb1 = as.numeric(log(Mishkin[,3]))

plot(tb1,type='l')

acf(tb1)

plot(diff(tb1),type='l')

acf(diff(tb1))

fit = auto.arima(tb1,ic="aic")

plot(fit$residuals)

acf(fit$residuals)

fit = arima(tb1, order = c(5,1,3))

plot(fit$residuals)

acf(fit$residuals)

fit = arima(tb1, order = c(0,1,1))

plot(fit$residuals)

acf(fit$residuals)

fit = auto.arima(tb1,ic="bic") # So much different between AIC and BIC

plot(fit$residuals)

acf(fit$residuals)
