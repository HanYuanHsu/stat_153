#Simulated Data for Illustrating Model Selection
#The problem involves figuring out the right order k for fitting k sinusoids (at frequencies 1/n,...,k/n) to the observed time series
#Generating the data:
n = 5000
ktrue = 20
sigtrue = 75
betatrue = c(0, rep(5, ktrue), rep(5, ktrue))
tme = 1:n
dt = rep(-1, n)
Xtrue = matrix(1, n, 1+(2*ktrue))
for(j in 1:ktrue)
{
    Xtrue[,(2*j)] = cos(2*pi*(j/n)*tme)
    Xtrue[,((2*j)+1)] = sin(2*pi*(j/n)*tme)
}
dt = Xtrue %*% betatrue + sigtrue*rnorm(n)
plot(dt, type = "l")

#Periodogram
plot(1:(n/2), abs(fft(dt)[2:((n/2)+1)])^2/n, type = "h", ylab = "Periodogram")
abline(h = 0)
#Note that the periodogram ordinates for the initial frequencies are quite high. 

#We  now calculate evidences using the Zellner prior (and then integration over tau):
K = 50
logevid = rep(-1, K)
Xin = matrix(1, n, 1)
for(kin in 1:K)
{
    Xin = cbind(Xin, cos(2*pi*(kin/n)*tme), sin(2*pi*(kin/n)*tme))
    mod = lm(dt ~ -1 + Xin)
    p = ncol(Xin)
    logevid[kin] = (lgamma(p/2)) - ((p/2)*(log((sum(mod$fitted.values^2))))) + (lgamma((n-p-1)/2)) - (((n-p)/2)*(log((sum(mod$residuals^2))/2)))
}
logevid = logevid - max(logevid)
evid = exp(logevid)
evid = evid/(sum(evid))
plot(1:K, evid, type = "h")
which.max(evid)
evid[which.max(evid)]

#Linear Regression Example One:
library(faraway)
data(seatpos)
help(seatpos)
names(seatpos)
pairs(seatpos)
g <- lm(hipcenter ~ ., seatpos)
summary(g)

fullX = model.matrix(g)
numvar = 8
all01 = expand.grid(replicate(numvar, 0:1, simplify = FALSE))
all01 = cbind(rep(1, nrow(all01)), all01)
colnames(all01) = colnames(fullX)
logEvid = rep(-1, nrow(all01))

for(mm in 1:nrow(all01))
{
    inds = all01[mm,]
    Xmat = fullX[,(inds == 1)]
    if(mm == 1) {Xmat = as.matrix(rep(1, nrow(seatpos)), nrow(seatpos), 1)}
    p = ncol(Xmat)
    n = nrow(Xmat)
    md = lm(seatpos$hipcenter ~ -1 + Xmat)
    logEvid[mm] = (lgamma(p/2)) - ((p/2)*(log((sum(md$fitted.values^2))))) + (lgamma((n-p-1)/2)) - (((n-p)/2)*(log((sum(md$residuals^2))/2)))
}
logEvid.scaled = logEvid - max(logEvid) #scaling so I can take exponential
Evid = exp(logEvid.scaled)
Evid = Evid/(sum(Evid))
Evid
plot(Evid, type = "h")

#High Evidence Models:
high.evid.models = all01[which(Evid > 0.02),]
cbind(high.evid.models, Evid[which(Evid > 0.02)])
#Highest Evidence: 
m1 = lm(hipcenter~Ht, data = seatpos)
summary(m1)
#You may compare to full model summary: 
summary(g)
#Next highest evidence: 
m2 = lm(hipcenter~HtShoes, data = seatpos)
summary(m2)
#Third highest evidence:
m3 = lm(hipcenter~Leg, data = seatpos)
summary(m3)

#Usually this is done via the Regsubsets() function in R: 
library(leaps)
b = regsubsets(hipcenter ~ ., seatpos)
rs = summary(b)
rs$which
names(rs)
rs$cp #Best model is the one with 1 variable - Ht
rs$adjr2 #Best model is the one with 3 variables - Age, Height, Leg.
rs$bic #Best model is the one with 1 variable - Ht.

#Linear Regression Example Two
#BodyFat Dataset
body = read.delim("bodyfat_corrected.txt", header = TRUE, sep = "")
#The descriptions of the variables:
#Density determined from underwater weighing
#Percent body fat from Siri's (1956) equation
#Age (years)
#Weight (lbs)
#Height (inches)
#Neck circumference (cm)
#Chest circumference (cm)
#Abdomen 2 circumference (cm)
#Hip circumference (cm)
#Thigh circumference (cm)
#Knee circumference (cm)
#Ankle circumference (cm)
#Biceps (extended) circumference (cm)
#Forearm circumference (cm)
#Wrist circumference (cm)

dim(body)
head(body)
tail(body)
names(body)

summary(body)
body[body$HEIGHT < 30,] #The height of this person is 29.5 which is probably a recording error. 
par(mfrow = c(4, 4))
for(i in 2:15)
{
	hist(body[,i], xlab = "", main = names(body)[i], breaks = 500)
}
par(mfrow = c(1, 1))
pairs(body[,-1])
ou1 = which(body$HEIGHT < 30)
ou2 = which(body$WEIGHT > 300)
ou3 = which(body$ANKLE > 30)
ou4 = which(body$HIP > 120)
ou = c(ou1, ou2, ou3, ou4)
ou

#Let us just remove these indices for simplicity: 
body = body[-ou,]
body = body[,-1]
head(body)
pairs(body)

#Full linear model: 
g = lm(BODYFAT ~ ., body)
summary(g)

#Bayesian model selection: 
fullX = model.matrix(g)
numvar = 13
all01 = expand.grid(replicate(numvar, 0:1, simplify = FALSE))
all01 = cbind(rep(1, nrow(all01)), all01)
colnames(all01) = colnames(fullX)
logEvid = rep(-1, nrow(all01))

for(mm in 1:nrow(all01))
{
    inds = all01[mm,]
    Xmat = fullX[,(inds == 1)]
    if(mm == 1) {Xmat = as.matrix(rep(1, nrow(body)), nrow(body), 1)}
    p = ncol(Xmat)
    n = nrow(Xmat)
    md = lm(body$BODYFAT ~ -1 + Xmat)
    logEvid[mm] = -((n/2)*(2*pi)) + (lgamma(p/2)) - ((p/2)*(log((sum(md$fitted.values^2))/2))) + (lgamma((n-p-1)/2)) - (((n-p)/2)*(log((sum(md$residuals^2))/2)))
}
logEvid.scaled = logEvid - max(logEvid) #scaling so I can take exponential
Evid = exp(logEvid.scaled)
Evid = Evid/(sum(Evid))
Evid
plot(Evid, type = "h")

head(sort(Evid, decreasing = T), 10)
high.evid.models = all01[which(Evid > 0.025),]
cbind(high.evid.models, Evid[which(Evid > 0.025)])
#Highest Evidence: 
m1 = lm(BODYFAT~WEIGHT + ABDOMEN + WRIST, data = body)
summary(m1)
#You may compare to full model summary: 
summary(g)
#Next highest evidence: 
m2 = lm(BODYFAT~HEIGHT + ABDOMEN + WRIST, data = body)
summary(m2)
#Third highest evidence:
m3 = lm(BODYFAT~WEIGHT + ABDOMEN, data = body)
summary(m3)

#Compare to regsubets
library(leaps)
b = regsubsets(BODYFAT ~ ., body, nvmax = 13)
rs = summary(b)
rs$which
names(rs)
rs$cp #Best model is the one with 7 variables - AGE, HEIGHT, NECK, CHEST, ABDOMEN, FOREARM, WRIST.
rs$adjr2 #Best model is the one with 9 variables - AGE, HEIGHT, NECK, CHEST, ABDOMEN, HIP, THIGH, FOREARM, WRIST.
rs$bic #Best model is the one with 3 variables - WEIGHT, ABDOMEN, WRIST (matches with the best evidence model above).


#Back to Time Series 
#US Accidental Deaths:
dt = USAccDeaths
n = length(dt)
tme = 1:n
plot(tme, dt, type = "l")

totfreq = 18 #total number of frequencies considered (the frequencies considered are 1/n, 2/n,...,totfreq/n)

Xfull = matrix(1, n, 1)
for(k in 1:totfreq)
{
   Xfull = cbind(Xfull, cos(2*pi*(k/n)*tme), sin(2*pi*(k/n)*tme))
}
all01 = expand.grid(replicate(totfreq, 0:1, simplify = FALSE))

logEvid = rep(-1, nrow(all01))

for(mm in 1:nrow(all01))
{
    inds = all01[mm,]
    indfreqs = rep(inds, each = 2)
    colsinclude = c(TRUE, (indfreqs == 1))
    Xmat = Xfull[,colsinclude]
    if(mm == 1) {Xmat = as.matrix(rep(1, n), n, 1)}
    p = ncol(Xmat)
    n = nrow(Xmat)
    md = lm(dt ~ -1 + Xmat)
    logEvid[mm] = (lgamma(p/2)) - ((p/2)*(log((sum(md$fitted.values^2))))) + (lgamma((n-p-1)/2)) - (((n-p)/2)*(log((sum(md$residuals^2))/2)))
}
logEvid.scaled = logEvid - max(logEvid) #scaling so I can take exponential
Evid = exp(logEvid.scaled)
Evid = Evid/(sum(Evid))
plot(Evid, type = "h")

head(sort(Evid, decreasing = T))
#Top three models:
high.evid.models = all01[which(Evid > 0.01),]
cbind(high.evid.models, Evid[which(Evid > 0.01)])
#Best model has three frequencies 1/n (this is to explain the "trend"), 6/n = 6/72 = 1/12 and 12/n = 12/72 = 2/12.
#Here is the fit of the best model:
bmd = lm(dt ~ cos(2*pi*(1/n)*tme) + sin(2*pi*(1/n)*tme) + cos(2*pi*(6/n)*tme) + sin(2*pi*(6/n)*tme) + cos(2*pi*(12/n)*tme) + sin(2*pi*(12/n)*tme))
plot(tme, dt, type = "l", xlab = "Time", ylab = "Deaths")
points(tme, bmd$fitted.values, type = "l", col = "red")

#We can also check the evidence of some specific model of interest. Suppose we only care about the model with frequencies 6/n and 12/n:
md.ind = c(Var1 = 0, Var2 = 0, Var3 = 0, Var4 = 0, Var5 = 0, Var6 = 1, Var7 = 0, Var8 = 0, Var9 = 0, Var10 = 0, Var11 = 0, Var12 = 1, Var13 = 0, Var14 = 0, Var15 = 0, Var16 = 0, Var17 = 0, Var18 = 0)
which(colSums(t(all01) == md.ind) == ncol(all01))
Evid[2081] #this model gets small evidence. This is mainly because there is a small trend in the data and the lowest frequency (1/n) is needed to fit this trend. 
md = lm(dt ~ cos(2*pi*(6/n)*tme) + sin(2*pi*(6/n)*tme) + cos(2*pi*(12/n)*tme) + sin(2*pi*(12/n)*tme))
points(tme, md$fitted.values, type = "l", col = "blue") #the red fit is clearly better than the blue fit 

#Google Trends (monthly data) for the query "Amazon" 
amazon.raw = read.delim("amazontrends25Aug2022.csv")
#the useful data is between row 2 and row 225
amazon.use = amazon.raw[2:225,1]
len = length(amazon.use)
amazon = rep(0, len)
for(i in 1:len)
{
  amazon[i] = as.numeric(unlist(strsplit(as.character(amazon.use[i]), ","))[2])
}
amazon.ts = ts(amazon, start = c(2004, 1), end = c(2022, 8), frequency = 12)
plot(amazon, type = "l", ylab = "Amazon Search Popularity", xlab = "Month", main = "Google Trends Data for Amazon") #this dataset shows both trend and seasonality.

n = length(amazon)
d = 12
tme = 1:n
scaledtme = (1:n)/n
dt = amazon
D = 8 #maximum polynomial degree
Xfull = matrix(1, n, 1)
for(i in 1:D)
{
   Xfull = cbind(Xfull, scaledtme^i)
}
totfreq = 5
for(k in 1:totfreq)
{
   Xfull = cbind(Xfull, cos(2*pi*(k/d)*tme), sin(2*pi*(k/d)*tme))
}
all01 = expand.grid(replicate((D + totfreq), 0:1, simplify = FALSE))
dim(all01)

logEvid = rep(-1, nrow(all01))
for(mm in 1:nrow(all01))
{
    inds = all01[mm,]
    indfreqs = rep(inds[(D+1):(D+totfreq)], each = 2)
    colsinclude = c(TRUE,(inds[1:D] == 1), (indfreqs == 1))
    Xmat = Xfull[,colsinclude]
    if(mm == 1) {Xmat = as.matrix(rep(1, n), n, 1)}
    p = ncol(Xmat)
    n = nrow(Xmat)
    md = lm(dt ~ -1 + Xmat)
    logEvid[mm] = (lgamma(p/2)) - ((p/2)*(log((sum(md$fitted.values^2))))) + (lgamma((n-p-1)/2)) - (((n-p)/2)*(log((sum(md$residuals^2))/2)))
}
logEvid.scaled = logEvid - max(logEvid) #scaling so I can take exponential
Evid = exp(logEvid.scaled)
Evid = Evid/(sum(Evid))
plot(Evid, type = "h")

head(sort(Evid, decreasing = T), 30)
#Top three models:
high.evid.models = all01[which(Evid > 0.01),]
cbind(high.evid.models, Evid[which(Evid > 0.01)])

#Best model (with evidence of 34 percent) is 3rd degree polynomial (with linear term set to zero) and three frequencies 1/12, 2/12, 3/12:
bmd2 = lm(dt ~ scaledtme + I(scaledtme^2) + I(scaledtme^3) + I(scaledtme^3) + I(scaledtme^4) + I(scaledtme^5) + I(scaledtme^6) + I(scaledtme^7) + I(scaledtme^8) + cos(2*pi*(1/d)*tme) + sin(2*pi*(1/d)*tme) + cos(2*pi*(2/d)*tme) + sin(2*pi*(2/d)*tme) + cos(2*pi*(3/d)*tme) + sin(2*pi*(3/d)*tme))
plot(c(tme, tme), c(dt, bmd2$fitted.values), type = "n", xlab = "Time", ylab = "Trends")
points(tme, dt, type = "l", xlab = "Time", ylab = "Trends")
points(tme, bmd2$fitted.values, type = "l", col = "red")









