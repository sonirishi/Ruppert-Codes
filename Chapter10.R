setwd("E:/Documents/Practice")
library(dplyr)
library(Ecdat)
library(zoo)

data("Hstarts")

plot(Hstarts)

acf(as.numeric(Hstarts[,1]))

meaned = rollapply(Hstarts[,1],10,mean)
sdd = rollapply(Hstarts[,1],10,sd)

plot(Hstarts[,1])
lines(meaned,col="red")

plot(sdd)
plot(meaned,col="red")

plot(log(Hstarts[,1]))

acf(log(Hstarts[,1]))

boxplot(Hstarts[,1],~cycle(Hstarts[,1]))

plot(diff(diff(Hstarts[,1]),4))

tss = diff(Hstarts[,1],4)

meaned = rollapply(tss,10,mean)

plot(tss,col="red")
lines(meaned,col="black")

acf(as.numeric(tss)); pacf(as.numeric(tss))

auto.arima(tss)

fit = arima(tss, order=c(1,0,2),(seasonal = list(order=c(0,0,2),period=4)))

plot(fit$residuals)

resid_roll = rollmean(fit$residuals,10,mean)

plot(resid_roll)
lines(tss,col="red")

fit = arima(Hstarts[,1], order=c(1,1,1),(seasonal = list(order=c(0,1,1),period=4)))

pred_hs = predict(fit,12)

plot(pred_hs$pred)

forecast_hs = forecast(fit,12)

plot(forecast_hs)

library(FitAR)

#BoxCox.Arima(fit)

cpi = read.csv("CPI.csv"); 
ip = read.csv("IP.csv")

acf(cpi$CPI)

plot(cpi$CPI)

plot(ip$IP)

ccf(cpi$CPI,ip$IP)

cor(cpi$CPI[769:900],ip$IP[697:828]) # as there is time trend, correlation seems very high

cor(diff(cpi$CPI[769:900]),diff(ip$IP[697:828])) # reality correlation without time trend

ccf(as.numeric(diff(cpi$CPI)),diff(ip$IP))

ts = rbind(diff(log(cpi$CPI))[769:900],diff(log(ip$IP))[697:828])

ts = t(ts)

multifit = ar(ts); multifit$aic  # 5 is the minimum AIC hence ar(5) process

#phi = matrix(c(0.767,-.33,0.0112,.3014),2,2)

multifit$ar; multifit$var.pred

eigen(phi)

multifit_n = ar(ts,order.max = 1)

eigen(multifit_n$ar[,,])

predict_multifit = predict(multifit_n,10)

library(longmemo)

sr1 = simARMA0(2500,0.7-1+0.5) #(count, hurst parameter H);  H = d + 1/2
sr2 = simARMA0(2500,0.35+0.5)
sr3 = simARMA0(2500,-0.35+0.5)

plot(sr3); acf(sr3)

plot(sr2); acf(sr2)

plot(cumsum(sr1),type='l'); acf(cumsum(sr1))   # This is equivalent to integrating the series

library(fracdiff)

acf(diffseries(cumsum(sr1),0.7))

acf(diffseries(cumsum(sr1),1))

################

infl_rt = Mishkin[,1]

fit1 = fracdiff(infl_rt,nar = 0, nma = 0)

diff_infl = diffseries(infl_rt,d = .4)

library(forecast)

auto.arima(diff_infl)

acf(diff_infl)

library(boot)

# tsboot(infl_rt)

################################

library(Ecdat)

data("IncomeUK")

consumption = IncomeUK[,2]

plot(consumption)

acf(consumption) # slowly dying long term process

plot(diff(consumption)) # seems seasonal

plot(diff(diff(consumption),4))

acf(as.numeric(diff(diff(consumption),4)))

pacf(as.numeric(diff(diff(consumption),4)))

arima(consumption,order=c(0,1,0),seasonal = c(c(1,1,1),4)) ## both c(1,1,1) and c(1,1,0) are ok

rollsd = rollapply(consumption,10,sd)

rollmean = rollapply(consumption,10,mean)

plot(consumption,col="red",type="l")  # some increase in variance is visible
lines(rollmean,col="orange",type="l")
par(new=T)
plot(rollsd,col="black",type="l")

rollmsq = rollmean^2

model = lm(rollsd ~ rollmean + rollmsq)  # not the best way but some intuition i guess

auto.arima(log(consumption))

library(lmtest)

gqtest(consumption~1)  # shows that variance is not stable in time series heteroskacidity test

modl = auto.arima(log(consumption),ic="bic")

forcs = forecast(modl,8)

plot(forcs)

forcs_pred = predict(modl,8)

plot(forcs,col="red")
lines(forcs_pred$pred,col="orange") # predict only gives forecast,forecast gives confidence interval

data(Tbrate,package = "Ecdat")

del_dat = diff(Tbrate)

var1 = ar(del_dat,order.max=4,aic=T)

acf(var1$resid[-1,])

data(Mishkin,package = "Ecdat")

plot(Mishkin[,5])

acf(Mishkin[,5])

plot(diff(sqrt(Mishkin[,5])))

acf(as.numeric(diff(sqrt(Mishkin[,5]))))

acf(as.numeric(diff(diff(Mishkin[,5]))))

pacf(as.numeric(diff(diff(Mishkin[,5]))))

gqtest(Mishkin[,5]~1)

rollmean = rollapply(Mishkin[,5],20,mean)
rollsd = rollapply(Mishkin[,5],20,sd)
rollmsq=rollmean^2

summary(lm(rollsd~rollmean+rollmsq)) #unsure of this now

library(fracdiff)

diffsq_mish = diff(sqrt(Mishkin[,5]))

fit.frac = fracdiff(diffsq_mish,nar=0,nma=0)

summary(fit.frac)

fit.frac$d

fdiff = diffseries(diffsq_mish,fit.frac$d)

acf(as.numeric(fdiff))

auto.arima(fdiff,ic="aic")
auto.arima(fdiff,ic="bic")

library(AER)

data("FrozenJuice")

price = FrozenJuice[,1]

plot(price)

rollmean = rollapply(price,30,mean)
rollsd = rollapply(price,30,sd)

plot(rollmean,col="red")
par(new=T)
plot(rollsd)

plot(log(price))

plot(diff(price))

plot(diff(log(price)))

auto.arima(price,ic="aic")

pacf(as.numeric(diff(price)))

acf(as.numeric(diff(price)))

auto.arima(price,ic="bic")

n = length(price)
set.seed(1988)
for(iter in 1:10){
  eps = rnorm(n+20)
  y = rep(0,n+20)
  for(t in 3:n+20){
    y[t] = 0.2825*y[t-1] + 0.057*y[t-2] + eps[t]
  }
  y = y[101:n+20]
  y = cumsum(y)  ## Integrated series of order 1 means difference will make it stationary
  y = ts(y,frequency = 12)
  fit=auto.arima(y,d=1,D=0,ic="bic")
  print(fit)
}

set.seed(1988)
niter=250
estimates=matrix(0,nrow=niter,ncol=2)
for (iter in 1:niter)
{
  eps = rnorm(n+20)
  y = rep(0,n+20)
  for (t in 3:(n+20))
  {
    y[t] = .2825 *y[t-1] + 0.0570*y[t-2] + eps[t] }
  y = y[101:n+20]
  y = cumsum(y)
  y = ts(y,frequency=12)
  fit=arima(y,order=c(2,1,0))
  estimates[iter,]=fit$coef
}

apply(estimates,2,mean)  #0.28203325 0.0522196

data(IncomeUK, package = "Ecdat")

income = IncomeUK[,1]

acf(as.numeric(income))
plot(income,type="l")
plot(diff(income))

pacf(as.numeric(diff(income)))

acf(as.numeric(diff(income)))

library(forecast)

auto.arima(income)

library(fracdiff)

fit = fracdiff(diff(income), nar=0,nma=0)

library(AER)

data("USMacroG")

plot(USMacroG[,"unemp"])

plot(diff(USMacroG[,"unemp"]))

acf(as.numeric(diff(USMacroG[,"unemp"])))

pacf(as.numeric(diff(USMacroG[,"unemp"])))

plot(diff(diff(USMacroG[,"unemp"])))

library(tseries)

adf.test(diff(USMacroG[,"unemp"]))

adf.test(diff(diff(USMacroG[,"unemp"]))) # ADF says this is stationary

adf.test(USMacroG[,"unemp"])  # not stationary

acf(as.numeric(diff(diff(USMacroG[,"unemp"]))))

pacf(as.numeric(diff(diff(USMacroG[,"unemp"]))))

arima(USMacroG[,"unemp"], order=c(0,2,0), seasonal = list(order=c(0,0,2),period=4))

arima(USMacroG[,"unemp"], order=c(1,1,0), seasonal = list(order=c(0,0,2),period=4))

auto.arima(USMacroG[,"unemp"])

set.seed(1988)

n = length(USMacroG[,"unemp"])

for(iter in 1:8){
  y = rep(0,n+50)
  error = rnorm(n + 50)
  for(t in 9:(n+50)){
    y[t] = 0.6611*y[t-1] - 0.4199*error[t-4] - 0.2623*error[t-8] + error[t]
    #print(y[t]); print(t)
  }
  y_new = y[40:length(y)]
  y_new = cumsum(y_new)
  y_new = ts(y_new,frequency = 4)
  print(auto.arima(y_new,d=1,D=0,ic="bic"))
}

library(Ecdat)

data("Tbrate")

auto.arima(Tbrate[,1])

auto.arima(Tbrate[,2])

auto.arima(Tbrate[,3])

plot(Tbrate[,2])

plot(diff(Tbrate[,2]))

adf.test(diff(Tbrate[,2]))

#kpss.test(diff(Tbrate[,2]))

acf(as.numeric(diff(Tbrate[,2])))

pacf(as.numeric(diff(Tbrate[,2])))

fit1 = arima(Tbrate[,2],order=c(1,1,2))

fit2 = arima(Tbrate[,2],order=c(1,1,1))

plot(fit1$residuals)

plot(fit2$residuals)

Box.test(fit1$residuals,lag=10,type="Ljung-Box",fitdf=1) ## null is independence

Box.test(fit2$residuals,lag=10,type="Ljung-Box",fitdf=1) ## null is independence

acf(as.numeric(fit1$residuals)) ## Seems ok to me

acf(as.numeric(fit2$residuals)) ## Seems autocorrelation a lag 1
