#Sunspots Data:
sunspots.data = read.delim("SN_y_tot_V2.0_25Aug2022.txt", header = F, sep = "")
head(sunspots.data)
sunspots = sunspots.data[,1:2]
plot(sunspots[,1], sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")

#Let us remove the last few years of the data (keeping it aside as test data) so we can compare the prediction accuracy of different models. 
splitnumber = 250 
sunspots.train = sunspots[1:splitnumber,]
sunspots.test = sunspots[(nrow(sunspots.train)+1):322,]
#We shall fit some models to the training data and then predict (on the basis of these models) the future observations; and compare the predictions with the actual values from the test dataset. It turns out that the different prediction accuracies depend crucially on how the training/test split is done (for example, compare splitnumber = 250 with splitnumber = 282 and splitnumber = 285). 

tme.train = 1700:(1700+nrow(sunspots.train)-1)
n = length(tme.train)
dt = sunspots.train[,2]

#We previously fitted a single sinusoid model (sinusoid(t) + noise) to this dataset (say, via the Bayesian posterior):
#Periodogram:
plot((1:(n/2))/n, abs(fft(dt)[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram", xlab = "Fourier Frequency")
abline(h = 0)

#Bayesian posterior
grid.res = 0.001
f.val = seq(0.01, 0.49, grid.res)
X = matrix(1, nrow = n, ncol = 3)
expos = rep(-1, length(f.val)) #exact marginal posterior for omega
log.values = rep(-1, length(f.val))
log.det.term = rep(-1, length(f.val))
for(i in 1:length(f.val))
{
    X[,2] = cos(2*pi*f.val[i]*tme.train)
    X[,3] = sin(2*pi*f.val[i]*tme.train)
    mod = lm(dt ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.det.term[i] = (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(f.val, expos, type = "l")

#Best estimate of f:
fhat = f.val[which.max(expos)]
#Period corresponding to fhat: 
1/fhat #this gives the 11 year cycle

#Fit of this model to the data:
x1 = cos(2*pi*tme.train*fhat)
x2 = sin(2*pi*tme.train*fhat)
mod = lm(dt ~ x1 + x2)
summary(mod)
plot(tme.train, dt, xlab = "Year (Training Data)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")
points(tme.train, mod$fitted.values, type = "l", col = "blue")

#Data generated from the model will not look like the observed sunspots data. 
sig.est = sqrt((sum(mod$residuals^2))/(mod$df)) #from the above linear model
simul.data = mod$fitted.values + sig.est*rnorm(length(tme.train))
plot(tme.train, simul.data, type = "l")
points(tme.train, mod$fitted.values, col = "red", type = "l")
#The sunspots data obviously does not look like this:
par(mfrow = c(4, 4))
for(i in 1:10){plot(tme.train, mod$fitted.values + sig.est*rnorm(length(tme.train)), type = "l", xlab = "Time", ylab = "Data")}
plot(tme.train, sunspots.train[,2], type = "l", xlab = "Time", ylab = "Data")
for(i in 12:16){plot(tme.train, mod$fitted.values + sig.est*rnorm(length(tme.train)), type = "l", xlab = "Time", ylab = "Data")} 
par(mfrow = c(1, 1))

#Prediction performance of the sinusoidal model at best fitting period: 
tme.test = (tme.train[length(tme.train)]+1):2021
preds.mod = (cos(2*pi*tme.test*fhat))*(mod$coefficients[2]) + (sin(2*pi*tme.test*fhat))*(mod$coefficients[3]) + mod$coefficients[1]
plot(tme.test, sunspots.test[,2], type = "l", xlab = "Year (1982 to 2021)", ylab = "Yearly Sunspot Numbers", main = "Prediction Evaluation")
points(tme.test, preds.mod, type = "l", col = "blue")

#Fitting the Yule Model to the training data: 
n = nrow(sunspots.train)
yt = dt[-c(1, 2)]
x1t = dt[-c(1, n)]
x2t = dt[-c(n-1, n)]
cbind(yt, x1t, x2t)
yulemod = lm(yt ~ x1t + offset((-1)*x2t))
summary(yulemod)
slpe = yulemod$coefficients[2] #this is the estimate of 2*cos(2 pi f)
yulefhat = (acos(slpe/2))/(2*pi)
1/yulefhat #this is the period of the sinusoid
#Note that this period is somewhat smaller than 11 years

#Fitted values of Yule model
plot(c(tme.train[-c(1, 2)], tme.train[-c(1,2)]), c(yt, yulemod$fitted.values), type = "n", xlab = "Year (Training Data)", ylab = "Yearly Sunspot Numbers")
points(tme.train[-c(1, 2)], yt, type = "l")
points(tme.train[-c(1, 2)], yulemod$fitted.values, type = "l", col = "blue")

#Prediction using Yule's regression: 
preds.yulemod = rep(-1, length(sunspots.test[,2]))
val1 = dt[n] 
val2 = dt[n-1]
preds.yulemod[1] = (yulemod$coefficients[1]) + ((yulemod$coefficients[2])*val1) - val2
preds.yulemod[2] = (yulemod$coefficients[1]) + ((yulemod$coefficients[2])*(preds.yulemod[1])) - (val1)
for(i in 3:length(preds.yulemod))
{
   preds.yulemod[i] = (yulemod$coefficients[1]) + ((yulemod$coefficients[2])*(preds.yulemod[(i-1)])) - (preds.yulemod[(i-2)])
}
plot(c(tme.test, tme.test, tme.test), c(sunspots.test[,2], preds.mod, preds.yulemod), type = "n", xlab = "Year (Test Data)", ylab = "Yearly Sunspot Numbers")
points(tme.test, sunspots.test[,2], type = "l")
points(tme.test, preds.mod, type = "l", col = "blue")
points(tme.test, preds.yulemod, type = "l", col = "red")    

c(sum((sunspots.test[,2] - preds.mod)^2), sum((sunspots.test[,2] - preds.yulemod)^2)) 

#Data generated from Yule Model: 
gam = 2*(cos(2*pi*yulefhat))
ysm = c(sunspots.train[1,2], sunspots.train[2,2], rep(-1, n-2))
sgm = sqrt((sum(yulemod$residuals^2))/(yulemod$df))
for(i in 3:n)
{
   ysm[i] = ((2-gam)*mean(sunspots.train[,2])) + (gam*ysm[(i-1)]) - (ysm[(i-2)]) + sgm*rnorm(1)
}
plot(tme.train, ysm, type = "l", yaxt = "n")

par(mfrow = c(4, 4))
for(i in 1:10)
{
    for(i in 3:n)
   {
      ysm[i] = ((2-gam)*mean(sunspots.train[,2])) + (gam*ysm[(i-1)]) - (ysm[(i-2)]) + sgm*rnorm(1)
   }
   plot(tme.train, ysm, type = "l", xlab = "Time", ylab = "Data", yaxt = "n")
}
plot(tme.train, sunspots.train[,2], type = "l", xlab = "Time", ylab = "Data", yaxt = "n")
for(i in 12:16)
{
    for(i in 3:n)
   {
      ysm[i] = ((2-gam)*mean(sunspots.train[,2])) + (gam*ysm[(i-1)]) - (ysm[(i-2)]) + sgm*rnorm(1)
   }
   plot(tme.train, ysm, type = "l", xlab = "Time", ylab = "Data", yaxt = "n")
}
par(mfrow = c(1, 1))

#AR(2) Model
modar2 = lm(yt ~ x1t + x2t)
summary(modar2)

preds.modar2 = rep(-1, length(sunspots.test[,2]))
val1 = dt[n] 
val2 = dt[n-1]
preds.modar2[1] = (modar2$coefficients[1]) + ((modar2$coefficients[2])*val1) + ((modar2$coefficients[3])*val2)
preds.modar2[2] = (modar2$coefficients[1]) + ((modar2$coefficients[2])*(preds.modar2[1])) + ((modar2$coefficients[3])*val1)
for(i in 3:length(preds.modar2))
{
   preds.modar2[i] = (modar2$coefficients[1]) + ((modar2$coefficients[2])*(preds.modar2[(i-1)])) + ((modar2$coefficients[3])*preds.modar2[(i-2)])
}

plot(c(tme.test, tme.test, tme.test, tme.test), c(sunspots.test[,2], preds.mod, preds.yulemod, preds.modar2), type = "n", xlab = "Year (Test Data)", ylab = "Yearly Sunspot Numbers")
points(tme.test, sunspots.test[,2], type = "l")
points(tme.test, preds.mod, type = "l", col = "blue")
points(tme.test, preds.yulemod, type = "l", col = "red")    
points(tme.test, preds.modar2, type = "l", col = "green")

c(sum((sunspots.test[,2] - preds.mod)^2), sum((sunspots.test[,2] - preds.yulemod)^2), sum((sunspots.test[,2] - preds.modar2)^2))

