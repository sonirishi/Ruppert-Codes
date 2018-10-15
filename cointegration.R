source("E:/Documents/Practice/standard_code.R")

library(tseries)

yieldDat = read.table("treasury_yields.txt",header=T)

plot(yieldDat$X3mo)

dat <- yieldDat[,3:7]

res <- residuals(lm(dat[,3]~dat[,1]+dat[,2]+dat[,4]+dat[,4]))

acf(res)  ## Autocorrelated

plot(res)

library(urca)

po.test(dat[,c(3,1,2,4,5)])

############## Simulate the series

error1 <- rnorm(n = 1000,0,1)
error2 <- rnorm(n = 1000,0,1)

y1 <- rep(0,1000)
y2 <- rep(0,1000)

for(i in 2:1000){
  y1[i] <- y1[i-1] + 0.5*(y1[i-1] - y2[i-1]) + error1[i]
  y2[i] <- y2[i-1] + 0.55*(y1[i-1] - y2[i-1]) + error2[i]
}

plot(y1,type='l')
lines(y2,col='blue',type = 'l')

plot(y1-y2,type = 'l',col="green")
