#Example for fitting MA(1)
library(astsa)
help(varve)
plot(varve, main = "varve data", type = "l")
logvarve = log(varve)
par(mfrow = c(2, 1))
plot(varve, main = "varve data", type = "l")
plot(logvarve, main = "logarithm of varve data", type = "l")
par(mfrow = c(1, 1))
plot(logvarve, main = "logarithm of varve data", type = "l")
#Because the level of the data seems to be changing with time, one cannot directly fit an MA model to this dataset: 
#A common strategy is to difference the data:
dt = diff(logvarve)
plot(dt, type = "l")
acf(dt) #this suggests fitting the MA(1) model
mamod = arima(dt, order = c(0, 0, 1))
mamod
mamod.orig = arima(logvarve, order = c(0, 1, 1))
mamod.orig
#there are different "Methods" in the arima function
mamod = arima(dt, order = c(0, 0, 1), method = "CSS")
mamod

#We can fit the MA(1) model directly (without using the custom function arima) as follows: 
Sfunc = function(alpha) #alpha consists of mu and the theta parameters (theta is indexed from 1; theta_0 is always 1)
{
    mu = alpha[1]
    thet = alpha[-1]
    n = length(dt)
    q = length(thet)
    zvals = c(rep(0, q), rep(9999, n))
    for(t in 1:n)
    {
       zvals[q+t] = dt[t] - mu - sum(thet * zvals[(q+t-1):t])
    }
    ans = sum(zvals^2)
    return(ans)
}
q = 1
alphaest = optim(rep(0, (q+1)), Sfunc)$par #We have chosen 0 as the initialization of both mu and theta
#Below uses an alternative initializer. Recall the correlation between y_t and y_{t+1} is thet/(1+thet^2).
# Since the sample acf at lag 1 approximates that correlation, we can get an approximate thet from the
# sample acf at lag 1 and use it as an initializer.
muinit = mean(dt)
rhohat = acf(dt, plot = F)$acf[2]
thetainit = (1 - (sqrt(1-4*(rhohat^2))))/(2*rhohat)
alphaest1 = optim(c(muinit, thetainit), Sfunc)$par #Both estimates of alphaest are almost the same

#Standard Error Calculation
#First compute Hessian
library(numDeriv)
H = hessian(Sfunc, alphaest)
n = length(dt)
sighat = sqrt(Sfunc(alphaest)/(n-length(alphaest)))
#c(sighat, sighat^2)
covmat = (sighat^2)*(solve(0.5*H))
#covmat
stderrs = sqrt(diag(covmat))
cbind(alphaest, stderrs)
sighat^2

#Prediction with MA(1)
L = 100 #number of future predictions desired
muest = alphaest[1]
thetaest = alphaest[2]
q = 1
n = length(dt)
zvals = c(rep(0, q), rep(9999, n))
for(t in 1:n)
{
   zvals[q+t] = dt[t] - muest - sum(thetaest * zvals[(q+t-1):t])
}
predest = muest + (thetaest*(zvals[(q+n)]))
predvec = c(predest, rep(muest, L - 1))
predse = c(sighat, rep(sighat*(sqrt(1+(thetaest^2))), L-1))
cbind(predvec, predse)
#Compare to the builtin function:
pma = predict(mamod, n.ahead = L)
cbind(pma$pred, pma$se)
#Plot predictions:
plot(c(dt, pma$pred), ylab = "Differenced Log Varve", xlab = "Years" , type = "l")
points((n+1) : (n+L), predvec, type = "l", col = "blue")
points((n+1) : (n+L), pma$pred, type = "l", col = "red")
points((n+1):(n+L), predvec + 2*predse, type = "l", col = "green")
points((n+1):(n+L), predvec - 2*predse, type = "l", col = "green")

#This is the prediction for the difference logvarve data. Predictions for the actual data can be obtained in the following way: 
predactualdata = logvarve[(n+1)] + cumsum(predvec)
plot(c(logvarve, predactualdata), type = "l", xlab = "Time (Years)", ylab = "Data")
points(1:length(logvarve), logvarve, type = "l")
points((1 + length(logvarve)): (L + length(logvarve)), predactualdata, type = "l", col = "blue")

#GDP Growth Rate (in percent)
gdpgr.raw = read.csv("A191RP1Q027SBEA.csv") #Data taken from https://fred.stlouisfed.org/series/A191RP1Q027SBEA 
#This is yearly data on the percent change in GDP from preceding period (year) with some seasonal adjustment
dt = gdpgr.raw[,2]
plot(dt, type = "l", main = "GDP Growth Rate")
acf(dt, lag.max = 50)  #Two spikes at lags 1 and 2. It might make sense to fit an MA(2) model:
mamod = arima(dt, order = c(0, 0, 2))
mamod
#We can fit the MA(2) model directly (without using the custom function arima) as follows: 
Sfunc = function(alpha) #alpha consists of mu and the theta parameters (theta is indexed from 1; theta_0 is always 1)
{
    mu = alpha[1]
    thet = alpha[-1]
    n = length(dt)
    q = length(thet)
    zvals = c(rep(0, q), rep(9999, n))
    for(t in 1:n)
    {
       zvals[q+t] = dt[t] - mu - sum(thet * zvals[(q+t-1):t])
    }
    ans = sum(zvals^2)
    return(ans)
}
q = 2
alphaest = optim(rep(0, (q+1)), Sfunc)$par
alphaest
#Hessian calculation
library(numDeriv)
H = hessian(Sfunc, alphaest)
n = length(dt)
sighat = sqrt(Sfunc(alphaest)/(n-length(alphaest)))
#c(sighat, sighat^2)
covmat = (sighat^2)*(solve(0.5*H))
#covmat
stderrs = sqrt(diag(covmat))
cbind(alphaest, stderrs)
sighat^2
#Both give similar results.
pacf(dt, lag.max = 50) 
armod = arima(dt, order = c(2, 0, 0))
armod
#The aic of mamod and armod are comparable. Actually the two models are similar. The intercept and the estimate of sigma are similar. The coefficients are also similar because: 
ARMAtoMA(armod$coef[c(1,2)], lag.max = 10)
#compare the first two coefficients to the fitted MA model:
cbind(mamod$coef[1:2], ARMAtoMA(armod$coef[c(1,2)], lag.max = 10)[1:2])
#Comparing the AutoCorrelation Functions of these models 
par(mfrow = c(2, 1))
plot(ARMAacf(ar = armod$coef[c(1, 2)], lag.max = 50), type = "h")
plot(ARMAacf(ma = mamod$coef[c(1, 2)], lag.max = 50), type = "h")
par(mfrow = c(1, 1))
#These plots have many similarities.


