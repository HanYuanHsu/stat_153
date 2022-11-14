#There is an option to get the pacf in the R function ARMAacf. 
#ACF of AR(1)
L = 15
corrs = ARMAacf(ar = c(0.7), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

par.corrs = ARMAacf(ar = c(0.7), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Partial Autocorrelation")  
abline(h=0)

#AR(1) with negative phi:
L = 15
corrs = ARMAacf(ar = c(-0.7), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

par.corrs = ARMAacf(ar = c(-0.7), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Partial Autocorrelation")  
abline(h=0)

#Sample Partial Autocorrelation Function
#For AR(1)
ph = 0.8
n = 200
dt = arima.sim(n = 200, list(ar = ph))
plot(1:n, dt, type = "o", xlab = "Time", ylab = "Time Series", main = "AR(1)") 

par(mfrow = c(1, 2))
acf(dt, lag.max = L, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = L,  plot = T, main = "Sample Partial Autocorrelation") 

#Negative Phi
ph = -0.8
n = 200
dt = arima.sim(n = 200, list(ar = ph))
plot(1:n, dt, type = "o", xlab = "Time", ylab = "Time Series", main = "AR(1)") 

par(mfrow = c(1, 2))
acf(dt, lag.max = L, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = L,  plot = T, main = "Sample Partial Autocorrelation") 

#For AR(2)
L = 25
ph1 = 0.5
ph2 = 0.25

par(mfrow = c(1, 2))
corrs = ARMAacf(ar = c(ph1, ph2), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

par.corrs = ARMAacf(ar = c(ph1, ph2), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

#From Data
dt = arima.sim(n = 200, list(ar = c(ph1, ph2)))
plot(1:n, dt, type = "o", xlab = "Time", ylab = "Time Series", main = "AR(2)") 

par(mfrow = c(1, 2))
acf(dt, lag.max = L, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = L,  plot = T, main = "Sample Partial Autocorrelation") 
par(mfrow = c(1, 1))

ph1 = 1
ph2 = -0.25
par(mfrow = c(1, 2))
corrs = ARMAacf(ar = c(ph1, ph2), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

par.corrs = ARMAacf(ar = c(ph1, ph2), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

#From Data
dt = arima.sim(n = 200, list(ar = c(ph1, ph2)))
plot(1:n, dt, type = "o", xlab = "Time", ylab = "Time Series", main = "AR(2)") 

par(mfrow = c(1, 2))
acf(dt, lag.max = L, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = L,  plot = T, main = "Sample Partial Autocorrelation") 
par(mfrow = c(1, 1))

ph1 = 1.5
ph2 = -0.75
par(mfrow = c(1, 2))
corrs = ARMAacf(ar = c(ph1, ph2), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

par.corrs = ARMAacf(ar = c(ph1, ph2), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

#From Data
dt = arima.sim(n = 200, list(ar = c(ph1, ph2)))
plot(1:n, dt, type = "o", xlab = "Time", ylab = "Time Series", main = "AR(2)") 

par(mfrow = c(1, 2))
acf(dt, lag.max = L, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = L,  plot = T, main = "Sample Partial Autocorrelation") 

#AR(3)
ph1 = 1
ph2 = -0.6
ph3 = 0.9
par(mfrow = c(1, 2))
corrs = ARMAacf(ar = c(ph1, ph2, ph3), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

par.corrs = ARMAacf(ar = c(ph1, ph2, ph3), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)
#These correlations seem to become larger and larger (also larger than 1). This means that we do not have a causal stationary solution here. 
polyroot(c(1, -ph1, -ph2, -ph3))
#Check that the first two root have magnitude smaller than one.
Mod(polyroot(c(1, -ph1, -ph2, -ph3)))

#Valid AR(3)
ph1 = 3/2
ph2 = -3/4
ph3 = 1/8
polyroot(c(1, -ph1, -ph2, -ph3))
Mod(polyroot(c(1, -ph1, -ph2, -ph3)))
par(mfrow = c(1, 2))
corrs = ARMAacf(ar = c(ph1, ph2, ph3), lag.max = L)
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)
par.corrs = ARMAacf(ar = c(ph1, ph2, ph3), lag.max = L, pacf = T)
plot(x = 1:L, y = par.corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation") 
abline(h = 0)

#From Data
dt = arima.sim(n = 200, list(ar = c(ph1, ph2, ph3)))
plot(1:n, dt, type = "o", xlab = "Time", ylab = "Time Series", main = "AR(2)") 

par(mfrow = c(1, 2))
acf(dt, lag.max = L, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = L,  plot = T, main = "Sample Partial Autocorrelation") 







