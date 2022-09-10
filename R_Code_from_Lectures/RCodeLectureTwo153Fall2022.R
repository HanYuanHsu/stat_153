#Let us fit a linear trend model to the US Population Data
uspop.raw = read.csv("POPTHM.csv")
#Data downloaded from FRED. Monthly Data. Data given for each month equals the average of the estimated population on the first day of the month and the first day of the next month. The units are thousands of dollars so 200,000 actually refers to 200 million. 
plot(uspop.raw$POPTHM, type = "l")
uspop.ts = ts(uspop.raw$POPTHM, start = c(1959, 1), end = c(2022, 6), frequency = 12)
plot(uspop.ts, ylab = "Population (in thousands)", xlab = "Time (in months)", main = "US Population")

#The Linear Trend model can be fit to this data in the following way:
t = 1: length(uspop.ts)
lin.model = lm(uspop.ts ~ 1 + t)
plot(t, uspop.ts, type = "l", xlab = "Time", ylab = "US Population", main = "Population of the United States")
points(t, lin.model$fitted, type = "l", col = "red")
summary(lin.model)

#Quadratic trend model:
t = 1: length(uspop.ts)
quad.model = lm(uspop.ts ~ 1 + t + I(t^2))
plot(t, uspop.ts, type = "l", xlab = "Time", ylab = "US Population", main = "Population of the United States")
points(t, lin.model$fitted, type = "l", col = "red")
points(t, quad.model$fitted, type = "l", col = "blue")
summary(quad.model)

#Which model would you prefer for this data: linear trend or quadratic trend? We shall address this question later using Bayesian model selection. 
