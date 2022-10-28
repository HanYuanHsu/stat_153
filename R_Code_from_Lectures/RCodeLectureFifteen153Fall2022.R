#Fitting AR models to data:
#Sunspots Data:
sunspots.data = read.delim("SN_y_tot_V2.0_25Aug2022.txt", header = F, sep = "")
head(sunspots.data)
sunspots = sunspots.data[,1:2]
plot(sunspots[,1], sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")

dt = sunspots[,2]
n = length(dt)
p = 2
yt = dt[(p+1):n]

Xmat = matrix(1, (n-p), 1)
for(j in 1:p)
{
   Xmat = cbind(Xmat, dt[(p-j+1):(n-j)])
}

modar = lm(yt ~ -1 + Xmat)
summary(modar)

#Predictions with the AR(p) model:
k = 100 #k predictions into the future
yhat = c(dt, rep(-9999, k))
for(i in 1:k)
{
    ans = modar$coefficients[1]
    for(j in 1:p)
    {
       ans = ans + (modar$coefficients[(j+1)])*yhat[n+i-j]
    }
    yhat[(n+i)] = ans
}
predvalues = yhat[-(1:n)]

plot((1:(n+k)), yhat, type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
abline(v = n)
abline(h = mean(dt))

#Prediction Uncertainty:
k = 100
resdf = n - 2*p - 1
sighat = sqrt((sum((modar$residuals)^2))/resdf) #this is also denoted by the Residual Standard Error
Gamhat = matrix(sighat^2, 1, 1) #this is the uncertainty for the first i.e., (n+1)^th prediction
#The following vector vkp is the vector a from the lecture notes
vkp = matrix(modar$coefficients[2], 1, 1) #this is the estimate for phi1
for(i in 1:(k-1))
{
    covterm = Gamhat %*% vkp
    varterm = (sighat^2) + (t(vkp) %*% (Gamhat %*% vkp))
    Gamhat = cbind(Gamhat, covterm)
    Gamhat = rbind(Gamhat, c(t(covterm), varterm))
    if (i < p) {vkp = c(modar$coefficients[(i+2)], vkp)}
    if (i >= p) {vkp = c(0, vkp)}
}
Gamhat
predsd = sqrt(diag(Gamhat))
#Plotting predictions with uncertainty bands (+/- 2 standard deviation bounds):
predlower = predvalues - 2*predsd
predupper = predvalues + 2*predsd

yhatlower = c(dt, predlower)
yhatupper = c(dt, predupper)

plot(c((1:(n+k)), (1:(n+k)), (1:(n+k))), c(yhat, yhatlower, yhatupper), type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
points((n+1):(n+k), predlower, type = "l", col = "red")
points((n+1):(n+k), predupper, type = "l", col = "red")
abline(v = n+1)

#Simulated AR(p) data:
dt = 100 + arima.sim(n = 500, list(ar = 0.7))
plot(dt, type = "l")

n = length(dt)
p = 1
yt = dt[(p+1):n]

Xmat = matrix(1, (n-p), 1)
for(j in 1:p)
{
   Xmat = cbind(Xmat, dt[(p-j+1):(n-j)])
}

modar = lm(yt ~ -1 + Xmat)
summary(modar)

#Plot of the fitted values
plot(c(1:(n-p), 1:(n-p)), c(yt, modar$fitted.values), type = "n", xlab = "Time", ylab = "Data")
points(1:(n-p), yt, type = "l")
points(1:(n-p), modar$fitted.values, type = "l", col = "blue")

#Predictions with the AR(p) model:
k = 100 #k predictions into the future
yhat = c(dt, rep(-9999, k))
for(i in 1:k)
{
    ans = modar$coefficients[1]
    for(j in 1:p)
    {
       ans = ans + (modar$coefficients[(j+1)])*yhat[n+i-j]
    }
    yhat[(n+i)] = ans
}
predvalues = yhat[-(1:n)]

plot((1:(n+k)), yhat, type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")

#Prediction Uncertainty:
resdf = n - 2*p - 1
sighat = sqrt((sum((modar$residuals)^2))/resdf) #this is also denoted by the Residual Standard Error
Gamhat = matrix(sighat^2, 1, 1) #this is the uncertainty for the first i.e., (n+1)^th prediction
vkp = matrix(modar$coefficients[2], 1, 1) #this is the estimate for phi1
for(i in 1:(k-1))
{
    covterm = Gamhat %*% vkp
    varterm = (sighat^2) + (t(vkp) %*% (Gamhat %*% vkp))
    Gamhat = cbind(Gamhat, covterm)
    Gamhat = rbind(Gamhat, c(t(covterm), varterm))
    if (i < p) {vkp = c(modar$coefficients[(i+2)], vkp)}
    if (i >= p) {vkp = c(0, vkp)}
}
Gamhat
predsd = sqrt(diag(Gamhat))
#Plotting predictions with uncertainty bands (+/- 2 standard deviation bounds):
predlower = predvalues - 2*predsd
predupper = predvalues + 2*predsd

yhatlower = c(dt, predlower)
yhatupper = c(dt, predupper)

plot(c((1:(n+k)), (1:(n+k)), (1:(n+k))), c(yhat, yhatlower, yhatupper), type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
points((n+1):(n+k), predlower, type = "l", col = "red")
points((n+1):(n+k), predupper, type = "l", col = "red")
abline(v = n+1)




#Another Dataset
sales.raw = read.csv("MRTSSM4453USN.csv")
dt = ts(sales.raw[,2], start = c(1992, 1), end = c(2022, 7), frequency = 12)
plot(dt, type = "l", ylab = "Millions of Dollars", main = "Retail Sales: Beer, Wine and Liquor Stores")

n = length(dt)
p = 36
yt = dt[(p+1):n]

Xmat = matrix(1, (n-p), 1)
for(j in 1:p)
{
   Xmat = cbind(Xmat, dt[(p-j+1):(n-j)])
}

modar = lm(yt ~ -1 + Xmat)
summary(modar)

#Plot of the fitted values
plot(c(1:(n-p), 1:(n-p)), c(yt, modar$fitted.values), type = "n", xlab = "Time", ylab = "Data")
points(1:(n-p), yt, type = "l")
points(1:(n-p), modar$fitted.values, type = "l", col = "blue")

#Predictions with the AR(p) model:
k = 100 #k predictions into the future
yhat = c(dt, rep(-9999, k))
for(i in 1:k)
{
    ans = modar$coefficients[1]
    for(j in 1:p)
    {
       ans = ans + (modar$coefficients[(j+1)])*yhat[n+i-j]
    }
    yhat[(n+i)] = ans
}
predvalues = yhat[-(1:n)]

plot((1:(n+k)), yhat, type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
abline(v = n)
abline(h = mean(dt))

#Prediction Uncertainty:
resdf = n - 2*p - 1
sighat = sqrt((sum((modar$residuals)^2))/resdf) #this is also denoted by the Residual Standard Error
Gamhat = matrix(sighat^2, 1, 1) #this is the uncertainty for the first i.e., (n+1)^th prediction
vkp = matrix(modar$coefficients[2], 1, 1) #this is the estimate for phi1
for(i in 1:(k-1))
{
    covterm = Gamhat %*% vkp
    varterm = (sighat^2) + (t(vkp) %*% (Gamhat %*% vkp))
    Gamhat = cbind(Gamhat, covterm)
    Gamhat = rbind(Gamhat, c(t(covterm), varterm))
    if (i < p) {vkp = c(modar$coefficients[(i+2)], vkp)}
    if (i >= p) {vkp = c(0, vkp)}
}
Gamhat
predsd = sqrt(diag(Gamhat))
#Plotting predictions with uncertainty bands (+/- 2 standard deviation bounds):
predlower = predvalues - 2*predsd
predupper = predvalues + 2*predsd

yhatlower = c(dt, predlower)
yhatupper = c(dt, predupper)

plot(c((1:(n+k)), (1:(n+k)), (1:(n+k))), c(yhat, yhatlower, yhatupper), type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
points((n+1):(n+k), predlower, type = "l", col = "red")
points((n+1):(n+k), predupper, type = "l", col = "red")
abline(v = n+1)


#Another Data:
hsp.raw = read.csv("ASPUS.csv")
dt = ts(hsp.raw[,2], start = c(1963, 1), end = c(2022, 2), frequency = 4)
plot(dt, type = "l", ylab = "Dollars", main = "Average Sales Price of Houses Sold")

n = length(dt)
p = 10
yt = dt[(p+1):n]

Xmat = matrix(1, (n-p), 1)
for(j in 1:p)
{
   Xmat = cbind(Xmat, dt[(p-j+1):(n-j)])
}

modar = lm(yt ~ -1 + Xmat)
summary(modar)

#Plot of the fitted values
plot(c(1:(n-p), 1:(n-p)), c(yt, modar$fitted.values), type = "n", xlab = "Time", ylab = "Data")
points(1:(n-p), yt, type = "l")
points(1:(n-p), modar$fitted.values, type = "l", col = "blue")

#Predictions with the AR(p) model:
k = 100 #k predictions into the future
yhat = c(dt, rep(-9999, k))
for(i in 1:k)
{
    ans = modar$coefficients[1]
    for(j in 1:p)
    {
       ans = ans + (modar$coefficients[(j+1)])*yhat[n+i-j]
    }
    yhat[(n+i)] = ans
}
predvalues = yhat[-(1:n)]

plot((1:(n+k)), yhat, type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
abline(v = n)
abline(h = mean(dt))

#Prediction Uncertainty:
resdf = n - 2*p - 1
sighat = sqrt((sum((modar$residuals)^2))/resdf) #this is also denoted by the Residual Standard Error
Gamhat = matrix(sighat^2, 1, 1) #this is the uncertainty for the first i.e., (n+1)^th prediction
vkp = matrix(modar$coefficients[2], 1, 1) #this is the estimate for phi1
for(i in 1:(k-1))
{
    covterm = Gamhat %*% vkp
    varterm = (sighat^2) + (t(vkp) %*% (Gamhat %*% vkp))
    Gamhat = cbind(Gamhat, covterm)
    Gamhat = rbind(Gamhat, c(t(covterm), varterm))
    if (i < p) {vkp = c(modar$coefficients[(i+2)], vkp)}
    if (i >= p) {vkp = c(0, vkp)}
}
Gamhat
predsd = sqrt(diag(Gamhat))
#Plotting predictions with uncertainty bands (+/- 2 standard deviation bounds):
predlower = predvalues - 2*predsd
predupper = predvalues + 2*predsd

yhatlower = c(dt, predlower)
yhatupper = c(dt, predupper)

plot(c((1:(n+k)), (1:(n+k)), (1:(n+k))), c(yhat, yhatlower, yhatupper), type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
points((n+1):(n+k), predlower, type = "l", col = "red")
points((n+1):(n+k), predupper, type = "l", col = "red")


#AR(1) with negative coefficients:
dt = 100 + arima.sim(n = 300, list(ar = -0.6))
plot(dt, type = "l")

n = length(dt)
p = 1
yt = dt[(p+1):n]

Xmat = matrix(1, (n-p), 1)
for(j in 1:p)
{
   Xmat = cbind(Xmat, dt[(p-j+1):(n-j)])
}

modar = lm(yt ~ -1 + Xmat)
summary(modar)

#Plot of the fitted values
plot(c(1:(n-p), 1:(n-p)), c(yt, modar$fitted.values), type = "n", xlab = "Time", ylab = "Data")
points(1:(n-p), yt, type = "l")
points(1:(n-p), modar$fitted.values, type = "l", col = "blue")

#Predictions with the AR(p) model:
k = 100 #k predictions into the future
yhat = c(dt, rep(-9999, k))
for(i in 1:k)
{
    ans = modar$coefficients[1]
    for(j in 1:p)
    {
       ans = ans + (modar$coefficients[(j+1)])*yhat[n+i-j]
    }
    yhat[(n+i)] = ans
}
predvalues = yhat[-(1:n)]

plot((1:(n+k)), yhat, type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
abline(v = n)
abline(h = mean(dt))

#Prediction Uncertainty:
resdf = n - 2*p - 1
sighat = sqrt((sum((modar$residuals)^2))/resdf) #this is also denoted by the Residual Standard Error
Gamhat = matrix(sighat^2, 1, 1) #this is the uncertainty for the first i.e., (n+1)^th prediction
vkp = matrix(modar$coefficients[2], 1, 1) #this is the estimate for phi1
for(i in 1:(k-1))
{
    covterm = Gamhat %*% vkp
    varterm = (sighat^2) + (t(vkp) %*% (Gamhat %*% vkp))
    Gamhat = cbind(Gamhat, covterm)
    Gamhat = rbind(Gamhat, c(t(covterm), varterm))
    if (i < p) {vkp = c(modar$coefficients[(i+2)], vkp)}
    if (i >= p) {vkp = c(0, vkp)}
}
Gamhat
predsd = sqrt(diag(Gamhat))
#Plotting predictions with uncertainty bands (+/- 2 standard deviation bounds):
predlower = predvalues - 2*predsd
predupper = predvalues + 2*predsd

yhatlower = c(dt, predlower)
yhatupper = c(dt, predupper)

plot(c((1:(n+k)), (1:(n+k)), (1:(n+k))), c(yhat, yhatlower, yhatupper), type = "n", xlab = "Time", ylab = "Data")
points((1:(n+k)), yhat, type = "l")
points(1:n, dt, type = "l")
points((n+1):(n+k), predvalues, type = "l", col = "blue")
points((n+1):(n+k), predlower, type = "l", col = "red")
points((n+1):(n+k), predupper, type = "l", col = "red")























