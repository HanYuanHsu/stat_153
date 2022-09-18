#In the last lecture, we studied simple models for trend and seasonality in time series datasets. These models were all fit by linear regression
#In some situations however, we are forced to deal with nonlinear regression. We shall learn how to fit a nonlinear regression model by Bayesian methods. 
#Consider the following concrete example:
sunspots.data = read.delim("SN_y_tot_V2.0_25Aug2022.txt", header = F, sep = "")
#Data description here: https://www.sidc.be/silso/infosnytot
#Column 1 is the year (2020.5 refers to the year 2020 for example); Column 2 is the yearly mean total sunspot number (this is obtained by taking a simple arithmetic mean of the daily total sunspot number over all days of the each year). Column 3 is the yearly mean standard deviation of the sunspot numbers from individual stations and Column 4 is the number of observations used to compute the yearly mean total sunspot number (-1 indicates missing value)
head(sunspots.data)
sunspots = sunspots.data[,1:2]
plot(sunspots[,1], sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")
#Some background: According to wikipedia, sunspots are temporary phenomena on the Sun's photosphere that appear as spots darker than the surrounding areas. Their number varies according to the approximately 11-year solar cycle. Data description (from http://www.sidc.be/silso/infosnytot): Yearly mean total sunspot number obtained by taking a simple arithmetic mean of the daily total sunspot number over all days of each year. (NB: in early years in particular before 1749, the means are computed on only a fraction of the days in each year because on many days, no observation is available).
#ARMA models originated through the analysis of this data (by Yule in 1923 who invented the AR(2) model while studying the Sunspots Data). 

#How do we model the seasonal trend this dataset?

#Because wikipedia says there is an approximately 11-year cycle for the sunspots, we can try to fit a sinusoid with period 11 years
tme = 1700:2021
freq = (2*pi*(1/11))
x1 = cos(tme*freq)
x2 = sin(tme*freq)
md = lm(sunspots[,2] ~ x1 + x2)
summary(md)

freq.alt = (2*pi*(1/10.5))
x1 = cos(tme*freq.alt)
x2 = sin(tme*freq.alt)
md1 = lm(sunspots[,2] ~ x1 + x2)
summary(md)

plot(tme, sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")
points(tme, md$fitted.values, type = "l", col = "red")
points(tme, md1$fitted.values, type = "l", col = "blue")

#The fit seems reasonable though not great. The following questions naturally arise here. What about other values of frequency? How about 10.5 year cycle or 11.5 year cycle? Do these periods also fit the data equally well? How was the approximate 11 year period determined in the first place? What is the uncertainty in this approximation?
#To answer these questions, it makes sense to work with a model with unknown frequency omega. This is however a nonlinear regression model. 

#Bayesian posterior for the nonlinear regression model
n = length(tme)
grid.res = 0.001
omga.val = seq(0.001, pi, grid.res)
X = matrix(1, nrow = n, ncol = 3)
expos = rep(-1, length(omga.val)) #exact marginal posterior for omega
log.values = rep(-1, length(omga.val))
log.det.term = rep(-1, length(omga.val))
for(i in 1:length(omga.val))
{
    X[,2] = cos(omga.val[i]*tme)
    X[,3] = sin(omga.val[i]*tme)
    mod = lm(sunspots[,2] ~ X)
    log.value = ((ncol(X) - n)/2)*(log(sum(mod$residuals^2))) - (0.5*(log(det(t(X) %*% X))))
    log.det.term[i] = (0.5*(log(det(t(X) %*% X))))
    log.values[i] = log.value
}
log.values = log.values - max(log.values) #scaling to remove large values 
expos = exp(log.values)
expos = (expos/sum(expos))/grid.res
plot(omga.val, expos, type = "l")

#The posterior is quite peaked. We can get a point estimate by the frequency which maximizes the posterior: 
ind.max = which.max(expos)
est.freq = omga.val[ind.max] #posterior maximizer
est.freq
#This corresponds to the following estimate of the period: 
period.est = (2*pi)/(omga.val[ind.max])
period.est

#Evaluate fit with this period:
x1 = cos(tme*est.freq)
x2 = sin(tme*est.freq)
mod = lm(sunspots[,2] ~ x1 + x2)
summary(mod)

plot(tme, sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")
points(tme, mod$fitted.values, type = "l", col = "red")

#We can also construct at uncertainty intervals containing say 95% of posterior mass (these are called credible intervals):  
postmass = function(s)
{
   return((sum(expos[(ind.max - s) : (ind.max + s)]))*grid.res)
}
lapply(1:10, postmass)
s = 2
c(omga.val[ind.max - s], omga.val[ind.max + s])
#This corresponds to the following interval for the period: 
c((2*pi)/(omga.val[ind.max + s]), (2*pi)/(omga.val[ind.max - s]))

#We can also look at the posterior mean and standard deviation: 
pm = (sum(omga.val*expos))*grid.res
pm
#Posterior standard deviation:
psd = sqrt((sum(((omga.val - pm)^2)*expos))*grid.res)
psd
c(pm - 2*psd, pm + 2*psd)
#corresponds to the following interval for the period: 
c((2*pi)/(pm+2*psd), (2*pi)/(pm-2*psd))

#Note the huge peak in the posterior. This basically means that if we are slightly off from the maximizing frequency, the fit to the data will be significantly deteriorated. For example: 
try.period = 10.7
try.freq = 2*pi*(1/try.period)
x1 = cos(tme*try.freq)
x2 = sin(tme*try.freq)
mod2 = lm(sunspots[,2] ~ x1 + x2)
summary(mod2)

plot(tme, sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")
points(tme, mod$fitted.values, type = "l", col = "red")
points(tme, mod2$fitted.values, type = "l", col = "blue")
#Even though 10.7 seems quite close to 11, the fit gets noticeably worse

#Approximate Computation (checking whether the X^T X matrix is diagonal)
omag = 1
k = 45
omag = 2*pi*(k/n) #Fourier Frequency
omag = 2*pi*(45.5/n)
X = matrix(1, nrow = n, ncol = 3)
X[,2] = cos(omag*tme)
X[,3] = sin(omag*tme)
t(X)%*%X


#Calculating the periodogram: 
grid.res = 0.001
omga.val = seq(0.201, pi, grid.res)
pgram = rep(-1, length(omga.val))
yy = sunspots[,2]
n = length(yy)
for(i in 1:length(omga.val))
{
    pgram[i] = (((sum(yy * cos(omga.val[i]*tme)))^2) + ((sum(yy * sin(omga.val[i]*tme)))^2))/n
}
plot(omga.val, pgram, type = "l")
plot(omga.val, log(pgram), type = "l")
omga.val[which.max(pgram)]


