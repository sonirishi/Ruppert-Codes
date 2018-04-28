setwd("E:/Documents/Practice")

library(Ecdat)
library(ggplot2)
library(dplyr)

data("Garch")

Garch_1 = Garch %>% mutate(rate_1 = lag(dm)) %>% mutate(rate = dm - rate_1)

plot(Garch_1$rate,type="l",col="blue")

data("Capm")

capm_1 = Capm %>% mutate(rate_1 = lag(rf)) %>% mutate(rate = rf - rate_1)

plot(capm_1$rate,type="l",col="blue")

data("SP500")

ndata = as.data.frame(rnorm(2783,mean(SP500$r500),sd(SP500$r500)));
colnames(ndata) = "normal"

ggplot() + geom_density(data=SP500,aes(x=r500),color="red") + 
  geom_density(data=ndata,aes(x=normal),color="blue")   ### Different dataframes in one plot

tdata = as.data.frame(rt(1000,df=10))

colnames(tdata) = "student"

f = uniroot(function(f) {1.5 - sqrt(f/(f-2))}, c(2.1,100))$root # Formula for t's SD

ndata = as.data.frame(rnorm(100,0,f));colnames(ndata) = "normal"

ggplot() + geom_density(data=tdata,aes(x=student),color="red") + 
  geom_density(data=ndata,aes(x=normal),color="blue")   ### Different dataframes in one plot

library(lmenssp)

qqplot.t(SP500$r500,dof=4)  ## qqplot with t distribution and df = 4

plot(SP500$r500,type="l")  ## high variation during 1807 index which as black monday

nw_data = as.data.frame(rnorm(100000,2,5) + jitter(rep(0,100000))) 

colnames(nw_data) = "col"

ggplot() + geom_density(data=nw_data,aes(x=col),color = "blue")

library(nortest)

nortest::ad.test(nw_data$col)

shapiro.test(SP500$r500)

Capm %>% mutate(lagrate = lag(rf)) %>% mutate(chg_rate = rf - lagrate) %>% 
  ggplot() + geom_point(aes(x=lagrate,y=chg_rate),color="blue")

Capm %>% mutate(loglag = lag(log(rf)), lograte = log(rf), lagrate = lag(rf)) %>% 
  mutate(chg_rate = lograte - loglag) %>% 
  ggplot() + geom_point(aes(x=lagrate,y=chg_rate),color="blue")

Capm %>% mutate(loglag = lag(log(rf)), lograte = log(rf), lagrate = lag(rf)) %>% 
  mutate(chg_rate = lograte - loglag) %>% 
  ggplot() + geom_point(aes(x=seq(1,516),y=chg_rate),color="blue")

capmdata = Capm %>% mutate(lag = lag(rf), chgrate = abs(rf - lag),chgrate_sq = (chgrate)^2)

boxcox_func = function(x,alpha){
  return((x^alpha-1)/alpha)
}

for(i in seq(0.01,1,0.01)){
  capmdata[,paste0("boxcoxchgrate",i)] = sapply(capmdata$chgrate, boxcox_func, alpha = i)
  capmdata[,paste0("boxcoxchgrate_sq",i)] = sapply(capmdata$chgrate_sq, boxcox_func, alpha = i)
}

corr_capm = as.data.frame(matrix(rep(0,300),nrow=100,ncol=3))

colnames(corr_capm) = c("corr_chng","corr_sqchng","alpha")

k = 1

for(i in seq(0.01,1,0.01)){
 corr_capm[k,1] = cor(capmdata[-1,paste0("boxcoxchgrate",i)],capmdata[-1,"lag"])
 corr_capm[k,2] = cor(capmdata[-1,paste0("boxcoxchgrate_sq",i)],capmdata[-1,"lag"])
 corr_capm[k,3] = i
 k = k + 1
}  

ggplot() + geom_line(data=corr_capm,aes(y=corr_chng, x=alpha,colour="blue")) + 
  geom_line(data=corr_capm,aes(y=corr_sqchng, x=alpha,colour="red")) + 
  scale_color_identity("Line.Color", labels=c("Sq Change","Abs Change"), guide="legend")

## aesthetics inside aes is linked to a variable, outside its static

data("Earnings")

earn = Earnings %>% mutate(sqrt_y = sqrt(y))

ggplot() + geom_density(data=earn,aes(x=y),color="blue")

ggplot() + geom_density(data=earn,aes(x=sqrt_y),color="red")

######## Lab exercise

data("EuStockMarkets")

stock_data = as.data.frame(matrix(rep(0,1860*4),nrow=1860,ncol=4));
colnames(stock_data) = c("ts1","ts2","ts3","ts4")

for(i in 1:4){
  stock_data[,paste0("ts",i)] = as.numeric(EuStockMarkets[,i])
}

ggplot(data=stock_data) + geom_line(aes(y=ts1,x=seq(1,1860)),color="black")

ggplot(data=stock_data) + geom_line(aes(y=ts2,x=seq(1,1860)),color="blue")

ggplot(data=stock_data) + geom_line(aes(y=ts3,x=seq(1,1860)),color="orange")

ggplot(data=stock_data) + geom_line(aes(y=ts4,x=seq(1,1860)),color="green")

##

logR = diff(log(EuStockMarkets))

plot(logR)

plot(as.data.frame(logR))

## Q2

index.names = dimnames(logR)[[2]]
par(mfrow=c(2,2))

for(i in 1:4){
  qqnorm(logR[,i],datax=T,main=index.names[i])
  qqline(logR[,i],datax=T)  ## Theoretical first and 3 quantile line through it which would be N
  print(shapiro.test(logR[,i]))
}

## Q3

n = dim(logR)[1]
q.grid = (1:n)/(n+1)  # vector of probab
df = c(1,4,6,10,20,30)
i=1 # Just study the DAX index 
par(mfrow=c(3,2))
for(j in 1:6){
  qqplot(logR[,i], qt(q.grid,df=df[j]), main=paste(index.names[i], ", df=", df[j]))
  abline(lm( qt(c(0.25,0.75),df=df[j]) ~ quantile(logR[,i],c(0.25,0.75)) ))
  ## first and 3 quantile line through them
}

## Q4

library("fGarch")
x = seq(-0.1,+0.1,by=0.001)
par(mfrow=c(1,1))
#plot( density(logR[,1]), lwd=2, ylim=c(0,60) )
plot( density(logR[,1]), lwd=2, ylim=c(0,20), xlim=c(-0.05,-0.0) )
lines(x,dstd(x,mean=median(logR[,1]),sd=mad(logR[,1]),nu=5),lty=5,lwd=2)
lines(x,dnorm(x,mean=median(logR[,1]),sd=mad(logR[,1])),lty=3,lwd=4)
legend( "topleft", c("KDE","t: df=5","normal"),lwd=c(2,2,4),lty=c(1,5,3))

ggplot(data=stock_data) + geom_density(aes(x=stock_data$ts1),color="blue")

## Exercises

ford = read.csv("ford.csv")

plot(ford$FORD,type="l")

mean(ford$FORD); median(ford$FORD); sd(ford$FORD)

qqnorm(ford$FORD,datax = TRUE)
qqline(ford$FORD,datax = TRUE)

shapiro.test(ford$FORD)

par(mfrow=c(3,2))

for (i in c(1,4,6,10,20,30)){
  qqplot(ford$FORD,rt(length(ford$FORD),df=i))
  qqline(ford$FORD,datax=TRUE)
}

## Remove Black monday 
par(mfrow=c(3,2))

for (i in c(1,4,6,10,20,30)){
  qqplot(ford[-which.min(ford$FORD),"FORD"],rt(length(ford[-which.min(ford$FORD),"FORD"]),df=i))
  qqline(ford[-which.min(ford$FORD),"FORD"],datax=TRUE)
}

F_inv_q = median(ford$FORD)
d = density(ford$FORD)
a = approx(d$x, d$y, F_inv_q, method = "linear")
q = 0.5

sqrt( ( q*(1-q) ) / ( length(ford$FORD) * a$y^2 ) ) 

sd(ford$FORD)/sqrt(length(ford$FORD))

## Q2

library(Ecdat)
library(dplyr)
library(ggplot2)
data("Garch")

garch1 = Garch %>% mutate(dy_lag = lag(dy), diff_lag = dy - dy_lag)

norm = as.data.frame(rnorm(nrow(garch1)-1,mean(garch1$diff_lag,na.rm=T),sd(garch1$diff_lag,na.rm=T)))
colnames(norm) = "norm"
norm$norm_1 = rnorm(nrow(garch1)-1,median(garch1$diff_lag,na.rm=T),mad(garch1$diff_lag,na.rm=T))

par(mfrow=c(1,1))
ggplot() + geom_density(data=garch1[-1,],aes(x=diff_lag,colour="diff_lag")) + 
  geom_density(data=norm,aes(x=norm,colour="norm")) +
  geom_density(data=norm,aes(x=norm_1,colour="norm_1")) + 
  xlab("Lag Difference") + ylab("Density") +
  ggtitle("Density Chart")

## Q4

garch1 = garch1 %>% mutate(lagbp = lag(bp), diffbp = bp - lagbp)

par(mfrow = c(3,2))

p = c(0.25,0.1,0.05,0.025,0.01,0.0025)

for(i in 1:6){
  qqnorm(garch1$diffbp,datax = T)
  qqline(garch1$diffbp,probs=c(p[i],1-p[i]),datax = T)
}


par(mfrow = c(3,2))

p = c(0.25,0.1,0.05,0.025,0.01,0.0025)

for(i in 1:6){
  rand = rnorm(length(garch1$diffbp))
  qqnorm(rand,datax = T)
  qqline(rand,probs=c(p[i],1-p[i]),datax = T)
}

garch1 = garch1 %>% mutate(loglagbp = log(lag(bp)), logdiffbp = log(bp) - loglagbp)

par(mfrow = c(3,2))

p = c(0.25,0.1,0.05,0.025,0.01,0.0025)

for(i in 1:6){
  qqnorm(garch1$logdiffbp,datax = T)
  qqline(garch1$logdiffbp,probs=c(p[i],1-p[i]),datax = T)
}
