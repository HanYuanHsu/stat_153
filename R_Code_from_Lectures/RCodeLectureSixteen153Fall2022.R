#US Population Dataset
uspop.raw = read.csv("POPTHM.csv")
#Data downloaded from FRED. Monthly Data. Data given for each month equals the average of the estimated population on the first day of the month and the first day of the next month. The units are thousands of dollars so 200,000 actually refers to 200 million. 
plot(uspop.raw$POPTHM, type = "l")
dt = uspop.raw$POPTHM
#Let us fit the AR(1) model:
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
#The coefficient phi1 is almost 1. So we can interpret this model as just fitting i.i.d Normals to the growth rates (growth rate = difference in population sizes of two successive years)
