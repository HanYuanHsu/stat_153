#The ARMAacf function:
#AR(1)
L = 20
ph = -0.8
corrs = ARMAacf(ar = c(ph), lag.max = L)
pacfcorrs = ARMAacf(ar = c(ph), lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation", main = "AR Process") 
abline(h=0)
plot(x = 0:L, y = c(0, pacfcorrs), type = "h",  xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))
c(corrs["1"], corrs["2"]/corrs["1"], corrs["3"]/corrs["2"], corrs["4"]/corrs["3"]) #all equal to phi

#MA(1)
L = 20
th = -0.8
corrs = ARMAacf(ma = c(th), lag.max = L)
pacfcorrs = ARMAacf(ma = c(th), lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation", main = "MA Process") 
abline(h=0)
plot(x = 0:L, y = c(0, pacfcorrs), type = "h",  xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#ARMA(1, 1)
L = 30
corrs = ARMAacf(ma = c(0.5), ar = c(-0.8), lag.max = L)
pacfcorrs = ARMAacf(ma = c(0.5), ar = c(-0.8), lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h",  xlab = "Lag k", ylab = "Autocorrelation", main = "ARMA(1, 1) Process") 
abline(h = 0)
plot(x = 0:L, y = c(0, pacfcorrs), type = "h",  xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))
c(corrs["1"], corrs["2"]/corrs["1"], corrs["3"]/corrs["2"], corrs["4"]/corrs["3"])


#Lake Huron Dataset
plot(LakeHuron, main = "Annual Measurements of the level of Lake Huron 1875-1972", xlab = "Year", ylab = "Level in feet", type = "o")
#Looks like there is a decreasing trend.
#Let us fit a line to it.

t = 1: length(LakeHuron)
lin.mod = lm(LakeHuron ~ 1 + t)
plot(1875:1972, LakeHuron, type = "o", xlab = "Year", ylab = "Level in Feet", main = "Lake Huron level measurements")
points(1875:1972, lin.mod$fitted, type = "l", col = "red")

#Residuals:
plot(1875:1972, lin.mod$residuals, type = "o", xlab = "Year", ylab = "Residuals after fitting the line")
sd(lin.mod$residuals)

#Do the residuals look purely random?
par(mfrow = c(2, 1))
acf(lin.mod$residuals, lag.max = 20, type = "correlation", plot = T, main = "Sample Autocorrelation of Residuals")
pacf(lin.mod$residuals, lag.max = 20, plot = T, main = "Sample Partial Autocorrelation of Residuals")
par(mfrow = c(1, 1))

#AR(2) should be a good model for residuals
armod2 = arima(lin.mod$residuals, order = c(2, 0, 0))
armod2
n = length(lin.mod$residuals)
bic2 = (-2)*(armod2$loglik) + (log(n))*(4) #there are four parameters in AR(2) including sigma

armod1 = arima(lin.mod$residuals, order = c(1, 0, 0))
armod1
bic1 = (-2)*(armod1$loglik) + (log(n))*(3) #there are three parameters in AR(1) including sigma

c(bic1, bic2)
c(armod1$aic, armod2$aic)

armod3 = arima(lin.mod$residuals, order = c(3, 0, 0))
armod3
bic3 = (-2)*(armod3$loglik) + (log(n))*(5) #there are five parameters in AR(3) including sigma

c(bic1, bic2, bic3)
c(armod1$aic, armod2$aic, armod3$aic)

#Predictions for the next 20 years:
k = 20
preds = predict(armod2, n.ahead = k)
plot(c(lin.mod$residuals, preds$pred), type = "l", xlab = "Time", ylab = "Data")
points((n+1):(n+k), preds$pred, col = "blue", type = "l")

#Predictions for original data:
n = length(LakeHuron)
predsorig = lin.mod$coefficients[1] + lin.mod$coefficients[2]*((n+1):(n+k)) + preds$pred
plot(c(LakeHuron, predsorig), type = "l", xlab = "Time", ylab = "LakeHuron Data")
points((n+1):(n+k), predsorig, col = "blue", type = "l")

#Differencing for this dataset:
dt = diff(LakeHuron)
plot(dt, type = "l")
#Let us fit a stationary model to this dataset:
par(mfrow = c(2, 1))
acf(dt, lag.max = 20, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dt, lag.max = 20, plot = T, main = "Sample Partial Autocorrelation")
par(mfrow = c(1, 1))

pmax = 5
qmax = 5
aicmat = matrix(-999, (1+pmax), (1+qmax))
bicmat = matrix(-999, (1+pmax), (1+qmax))
for(i in 0:pmax)
{
    for(j in 0:qmax)
    {
        md = arima(dt, order = c(i, 0, j))
        aicmat[(i+1), (j+1)] = ((-2)*md$loglik)+2*(i+j+2)
        bicmat[(i+1), (j+1)] = ((-2)*md$loglik)+(log(length(dt)))*(i+j+2)
    }
}

#Best AR model:
which.min(aicmat[,1]) #AR(3)
which.min(bicmat[,1]) #AR(0) or just white noise

#Best MA model:
which.min(aicmat[1,]) #MA(4) 
which.min(bicmat[1,])

#Best ARMA model:
which(aicmat == min(aicmat), arr.ind = T)
which(bicmat == min(bicmat), arr.ind = T)
#ARMA(1, 2) seems to be the best
armabest = arima(dt, order = c(1, 0, 2))
armabest
#Compare to 
arbest = arima(dt, order = c(3, 0, 0))
arbest
#and
mabest = arima(dt, order = c(0, 0, 4))
mabest

#Fit the ARMA(1, 2) model directly to the original data without differencing:
armaorig = arima(LakeHuron, order = c(1, 1, 2))
armaorig
predsarma = predict(armaorig, n.ahead = k)
plot(c(LakeHuron, predsorig), type = "l", xlab = "Time", ylab = "LakeHuron Data")
points((n+1):(n+k), predsorig, col = "blue", type = "l")
points(c(LakeHuron, predsarma$pred), type = "l", xlab = "Time", ylab = "LakeHuron Data")
points((n+1):(n+k), predsarma$pred, col = "red", type = "l")


#Sales for Beer, Wine and Liquor Stores
sales.raw = read.csv("MRTSSM4453USN.csv")
dt = sales.raw[,2]
plot(dt, type = "l", ylab = "Millions of Dollars", main = "Retail Sales: Beer, Wine and Liquor Stores")
n = length(dt)
#Differencing
dt12 = diff(dt, lag = 12)
plot(dt12, type = "l")
plot(diff(dt12), type = "l")
dtnew = diff(dt12) #this is the transformed data

par(mfrow = c(2, 1))
acf(dtnew, lag.max = 60, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dtnew, lag.max = 60, plot = T, main = "Sample Partial Autocorrelation")
par(mfrow = c(1, 1))

pmax = 5
qmax = 5
aicmat = matrix(-999, (1+pmax), (1+qmax))
bicmat = matrix(-999, (1+pmax), (1+qmax))
for(i in 0:pmax)
{
    for(j in 0:qmax)
    {
        md = arima(dtnew, order = c(i, 0, j))
        aicmat[(i+1), (j+1)] = ((-2)*md$loglik)+2*(i+j+2)
        bicmat[(i+1), (j+1)] = ((-2)*md$loglik)+(log(length(dtnew)))*(i+j+2)
    }
}

mamod = arima(dtnew, order = c(0, 0, 1))
mamod

#Best AR model:
which.min(aicmat[,1]) 
which.min(bicmat[,1]) 

#Best MA model:
which.min(aicmat[1,])
which.min(bicmat[1,])

#Best ARMA model:
which(aicmat == min(aicmat), arr.ind = T)
which(bicmat == min(bicmat), arr.ind = T)

armabest = arima(dtnew, order = c(4, 0, 3))
armabest

armaorig = arima(dt, order = c(4, 1, 3), seasonal = list(order = c(0, 1, 0), period = 12))
armaorig

L = 36

preds = predict(armaorig, n.ahead = L)
plot(c(dt, preds$pred), type = "l", xlab = "Time", ylab = "Data")
points((n+1):(n+L), preds$pred, type = "l", col = "blue")

#Compare to predictions by MA(4)
mabestorig = arima(dt, order = c(0, 1, 4), seasonal = list(order = c(0, 1, 0), period = 12))
mabestorig
predsma = predict(mabestorig, n.ahead = L)
#plot(c(dt, predsma$pred), type = "l", xlab = "Time", ylab = "Data")
points((n+1):(n+L), predsma$pred, type = "l", col = "red")

#Seasonal ARMA and ARIMA models:
#Motivation: co2 dataset. Taken from the book by Cryer and Chan.
library(TSA)
data(co2, package = "TSA")
help(co2, package = "TSA")
#There is a different dataset with the name co2 in the package datasets.

plot(co2, xlab = "year", ylab = "CO2 Levels", main = "monthly Carbon Dioxide Levels at Alert, NWT, Canada", type = "o")

#Both Trend and Seasonality:
#Seasonal difference:
t1 = diff(co2, 12)
plot(t1, main = "Seasonal Difference of co2")
#First difference of t1:
t2 = diff(t1)
plot(t2, main = "Seasonal and First Difference of co2") 
#No trend or seasonality now. We can now fit a stationary model. 

par(mfrow = c(2, 1))
acf(t2, lag.max = 40)
#Large negative autocorrelation at lag 1 and then a bunch of small autocorrelations followed by large autocorrelations at lags 11, 12, 13 (the one at lag 12 is quite large)
pacf(t2, lag.max = 40)
par(mfrow = c(1, 1))
#Not very interpretable.

#To fit a good model for this dataset, we need to learn about seasonal and multiplicative seasonal ARMA models.

#Seasonal MA(1) model
th = c(rep(0, 11), 0.7)#theta = 0.7
L = 50
corrs = ARMAacf(ma = th, lag.max = L)
par.corrs = ARMAacf(ma = th, lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#Seasonal MA(2) model
th = c(rep(0, 11), 0.7, rep(0, 11), 0.8)
L = 100
corrs = ARMAacf(ma = th, lag.max = L)
par.corrs = ARMAacf(ma = th, lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#Seasonal AR(1) model
ph = c(rep(0, 11), 0.7)#phi = 0.7
L = 150
corrs = ARMAacf(ar = ph, lag.max = L)
par.corrs = ARMAacf(ar = ph, lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#Seasonal AR(2) model
ph = c(rep(0, 11), 1.5, rep(0, 11), -0.75)
L = 150
corrs = ARMAacf(ar = ph, lag.max = L)
par.corrs = ARMAacf(ar = ph, lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#Seasonal ARMA(1,1) model
ph = c(rep(0, 11), 0.5)#phi = 0.4
th = c(rep(0, 11), 0.8) #th = 0.8
L = 150
corrs = ARMAacf(ar = ph, ma = th, lag.max = L)
par.corrs = ARMAacf(ar = ph, ma = th, lag.max = L, pacf = T)
par(mfrow = c(2, 1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))
#Therefore the acf and pacf of seasonal ARMA models look just like those of usual ARMA models at the seasonal lags.

#Multiplicative seasonal models
#ARMA(0, 1) X (1, 0)_12 model: regular MA(1) multiplied by seasonal AR(1)
ph = c(rep(0, 11), 0.8) #Phi = 0.8
L = 50
corrs = ARMAacf(ar = ph, ma = -0.5, lag.max = L)
par.corrs = ARMAacf(ar = ph, ma = -0.5, lag.max = L, pacf = T)
par(mfrow = c(2,1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))


#ARMA(0, 2) X (1, 0)_12 model: regular MA(2) multiplied  by seasonal AR(1)
ph = c(rep(0, 11), 0.8) #Phi = 0.8
L = 50
corrs = ARMAacf(ar = ph, ma = c(-0.5, 1), lag.max = L)
par.corrs = ARMAacf(ar = ph, ma = c(-0.5, 1), lag.max = L, pacf = T)
par(mfrow = c(2,1))
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#ARMA(0, 1) X (0, 1)_12 model: regular MA(1) multiplied by seasonal MA(1)
th = c(0.6, rep(0, 10), 0.7, 0.6*0.7)
corrs = ARMAacf(ma = th, lag.max = L)
par(mfrow = c(2, 1))
par.corrs = ARMAacf(ma = th, lag.max = L, pacf = T)
plot(x = 0:L, y = corrs, type = "h", xlab = "Lag k", ylab = "Autocorrelation")
abline(h = 0)
plot(x = 1:L, y = par.corrs, type = "h", xlab = "Lag k", ylab = "Partial Autocorrelation")
abline(h = 0)
par(mfrow = c(1, 1))

#Back to the co2 dataset:
#co2 dataset
data(co2, package = "TSA")
plot(co2, xlab = "year", ylab = "CO2 Levels", main = "monthly Carbon Dioxide Levels at Alert, NWT, Canada", type = "o")

#Both Trend and Seasonality:
#Seasonal difference:
t1 = diff(co2, 12)
plot(t1, main = "Seasonal Difference of co2")
#First difference of t1:
t2 = diff(t1)
plot(t2, main = "Seasonal and First Difference of co2")
#No trend or seasonality now. We can now fit a stationary model.

par(mfrow = c(2,1))
acf(t2, lag.max = 40)
#Large negative autocorrelation at lag 1 and then a bunch of small autocorrelations followed by large autocorrelations at lags 11, 12, 13 (the one at lag 12 is quite large)
pacf(t2, lag.max = 40)
#Not very interpretable.

#Multiplicative Seasonal ARIMA for the co2 dataset:
m1.co2 = arima(co2, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
m1.co2

#Prediction
m = 24
fcast = predict(m1.co2, n.ahead = m)

newx = 1:(length(co2) + m)
newy = c(co2, fcast$pred)

plot(newx, newy, type = "l")
points(newx[((length(co2)+1):(length(co2) + m))], newy[((length(co2)+1):(length(co2) + m))], col = "blue", type = "l" )

#Back to the Sales Data for Beer, Wine and Liquor Stores
sales.raw = read.csv("MRTSSM4453USN.csv")
dt = sales.raw[,2]
plot(dt, type = "l", ylab = "Millions of Dollars", main = "Retail Sales: Beer, Wine and Liquor Stores")
n = length(dt)
#Differencing twice:
dt12 = diff(dt, lag = 12)
dtnew = diff(dt12) #this is the transformed data

par(mfrow = c(2, 1))
acf(dtnew, lag.max = 60, type = "correlation", plot = T, main = "Sample Autocorrelation")
pacf(dtnew, lag.max = 60, plot = T, main = "Sample Partial Autocorrelation")
par(mfrow = c(1, 1))

#Previously we fit an ARMA(4, 3) model to this dataset. Now let us multiply it with a seasonal ARMA model
armaorig = arima(dt, order = c(4, 1, 3), seasonal = list(order = c(0, 1, 0), period = 12))
armaorig
L = 36
preds = predict(armaorig, n.ahead = L)
plot(c(dt, preds$pred), type = "l", xlab = "Time", ylab = "Data")
points((n+1):(n+L), preds$pred, type = "l", col = "blue")

#Multiplicative Seasonal ARIMA
sarimamod = arima(dt, order = c(4, 1, 3), seasonal = list(order = c(0, 1, 1), period = 12))
sarimamod #much improved likelihood and aic
preds.sarima = predict(sarimamod, n.ahead = L)
#plot(c(dt, predsma$pred), type = "l", xlab = "Time", ylab = "Data")
points((n+1):(n+L), preds.sarima$pred, type = "l", col = "green")

#Another Multiplicative Seasonal ARIMA model
sarimamod2 = arima(dt, order = c(4, 1, 3), seasonal = list(order = c(1, 1, 3), period = 12))
sarimamod2 #much improved likelihood and aic
preds.sarima = predict(sarimamod2, n.ahead = L)
#plot(c(dt, predsma$pred), type = "l", xlab = "Time", ylab = "Data")
points((n+1):(n+L), preds.sarima$pred, type = "l", col = "red")













