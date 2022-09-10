#Here are some simple time series datasets. 
#US Population Data
uspop.raw = read.csv("POPTHM.csv")
#Data downloaded from FRED. Monthly Data. Data given for each month equals the average of the estimated population on the first day of the month and the first day of the next month. The units are thousands of dollars so 200,000 actually refers to 200 million. 
plot(uspop.raw$POPTHM, type = "l")
uspop.ts = ts(uspop.raw$POPTHM, start = c(1959, 1), end = c(2022, 6), frequency = 12)
plot(uspop.ts, ylab = "Population (in thousands)", xlab = "Time (in months)", main = "US Population")
#Natural questions based on this data: explain the trend in the data, what is the rate of growth of the population, forecasting etc.                                         

#Here is another old population dataset that is commonly used in time series textbooks
data(lynx)
help(lynx)
plot(lynx)

#Economic time series

#Unemployment Rate (from FRED)
unrate.raw = read.csv("UNRATE.csv")
unrate.ts = ts(unrate.raw$UNRATE, start = c(1948, 1), end = c(2022, 7), frequency = 12)
plot(unrate.ts)
#This data is seasonally adjusted. Unemployment rate represents the number of unemployed as a percentage of the labor force (labor force = essentially people 16 or older). 
#There seem to be certain cyclical trends in the unemployment rate (note that this data is seasonally adjusted). Understanding the nature of these "business" cycles will be useful. 

#Many other economic datasets are available at the FRED website

#Physical Time Series
#This dataset is downloaded from http://www.sidc.be/silso/datafiles#total (silso stands for sunspot index and long term solar observations)
sunspots.data = read.delim("SN_y_tot_V2.0_25Aug2022.txt", header = F, sep = "")
#Data description here: https://www.sidc.be/silso/infosnytot
#Column 1 is the year (2020.5 refers to the year 2020 for example); Column 2 is the yearly mean total sunspot number (this is obtained by taking a simple arithmetic mean of the daily total sunspot number over all days of the each year). Column 3 is the yearly mean standard deviation of the sunspot numbers from individual stations and Column 4 is the number of observations used to compute the yearly mean total sunspot number (-1 indicates missing value)
head(sunspots.data)
sunspots = sunspots.data[,1:2]
plot(sunspots[,1], sunspots[,2], xlab = "Year (1700 to 2021)", ylab = "Yearly Sunspot Numbers", type = "l", main = "Sunspot Data")
#Some background: According to wikipedia, sunspots are temporary phenomena on the Sun's photosphere that appear as spots darker than the surrounding areas. Their number varies according to the approximately 11-year solar cycle. Data description (from http://www.sidc.be/silso/infosnytot): Yearly mean total sunspot number obtained by taking a simple arithmetic mean of the daily total sunspot number over all days of each year. (NB: in early years in particular before 1749, the means are computed on only a fraction of the days in each year because on many days, no observation is available).
#ARMA models originated through the analysis of this data (by Yule in 1923 who invented the AR(2) model while studying the Sunspots Data). 

#Temperature Anomaly Data: 
#This dataset shows the change in global surface temperature compared to the long-term average from 1951 to 1980. I have downloaded this data from https://climate.nasa.gov/vital-signs/global-temperature/
tempdata = read.delim("TemperatureAnomaly25Aug2022.txt")
tempdata = tempdata[-c(1, 2, 3),]
len = length(tempdata)
temp.ts = rep(0, len)
temp.year = rep(0, len)
temp.loess = rep(0, len)

for(i in 1:len)
{
    temp.ts[i] =  as.numeric(unlist(strsplit(tempdata[i], "\\s+")))[2]
    temp.year[i] =  as.numeric(unlist(strsplit(tempdata[i], "\\s+")))[1]
    temp.loess[i] =  as.numeric(unlist(strsplit(tempdata[i], "\\s+")))[3]
    
}
plot(temp.year,temp.ts, type = "l", xlab = "Temperature Anomaly", ylab = "Year")
cbind(temp.year, temp.ts) #According to the NASA website, the year 2020 tied with 2016 for the hottest year on record since recordkeeping began in 1880. 
#A natural question with this dataset involves obtaining a smooth estimate of the underlying trend. Such an estimate is already provided by the NASA website using the method of LOESS. 
points(temp.year, temp.loess, type = "l", col = "red")
#We shall learn how to obtain such smooth estimates. 

#Southern Oscillation Index Data
#El Nino (see wiki entry) is a climate pattern that describes the unusual warming of surface waters in the eastern tropical Pacific Ocean. It is believed to occur irregularly at two to seven year intervals. 
#Closely related to the fluctuations in oceanic temperatures are large-scale changes in atmospheric pressure. El Nino events are associated with sustained negative Southern Oscillation Index (SOI) values. The SOI is computed from fluctuations in the surface air pressure difference between Tahiti and Darwin. 
soi = read.table("soi_25Aug2022.txt", header = T)
rownames(soi) = soi$YEAR
soi = soi[,-1]
soi.ts = ts(as.vector(t(soi)), start = c(1951, 1), end = c(2022, 7), frequency = 12)
plot(soi.ts, type = "h", xlab = "Year", ylab = "Southern Oscillation Index")
abline(h = 0)
abline(h = -8)
#Prolonged periods of negative SOI values correspond to El Nino. 
#Natural question of interest: currently the SOI values are positive; when will they dip back into negative (i.e., when will the next El Nino occur?

#A Financial Dataset
#Install and load the R package quantmod. 
library(quantmod)
#We can obtain the daily stock prices of Apple stock from Yahoo Finance using the following command: 
getSymbols("AAPL") #AAPL is the ticker or stock symbol for Apple Inc.
dim(AAPL)
head(AAPL)
tail(AAPL)
#Let us plot the daily closing prices in the variable AAPL.Adjusted (adjusted price takes into account stock splits etc.; this adjustment is done retroactively on Yahoo Finance). 
plot(AAPL$AAPL.Adjusted, type = "l")

#Modeling (and especially predicting) stock prices is very challenging. 
#Financial analysts also study the volatility of stock prices in terms of returns (as opposed to the stock prices themselves). How are the returns calculated? There are many ways of doing this. The most popular one is to take: R(t) = log(P(t)/P(t-1)), t = 2, 3, ... where P(t) denotes the adjusted stock price after day t. Because P(t) and P(t-1) are usually close to each other, we would have that log(P(t)/P(t-1)) equals approximately R'(t) = (P(t) - P(t-1))/P(t-1). R(t) is preferred to R'(t) as the return measure of choice. 
#Daily returns of apple stock. 
quartz()
AAPL.rtn = diff(log(AAPL$AAPL.Adjusted))
plot(AAPL.rtn, type = "l")
#Two observations from this time series: (a) There exist some large outlying observations (this means that the returns have heavy tails), (b) the returns were volatile in certain periods but stable in others (this is called volatility clustering). 
#A problem in the analysis of these types of financial data is to forecast the volatility of future returns. Models such as ARCH, GARCH and stochastic volatility models have been developed to handle these problems. Time permitting, we shall study stochastic volatility models in the later part of the class. 

#Another financial data
#This dataset gives monthly returns on the four assets (stocks of MOBIL, IBM, WEYER, CITICRP, overall market returns (average returns of all stocks listed at the New York and American exchanges) and returns of a risk free asset
capm.ts = ts(read.table("capm.txt", header = T), start = c(1978, 1), frequency = 12) #monthly returns (NOT sure if they are logged; probably not)
capm.ts
colnames(capm.ts)
plot(capm.ts)
WEYER = capm.ts[,"WEYER"] - capm.ts[,"RKFREE"]
x = capm.ts[,"MARKET"] - capm.ts[,"RKFREE"]
plot(as.vector(x), as.vector(WEYER), xlab = "Market Return", ylab = "WEYER Return", main = "Performance of WEYER relative to MARKET") #a regression will not give any information on the time varying nature of the performance.

#Google trends is another source of time series data. Go to google.com/trends and type in a query, say "amazon". We get data on the popularity of the query "amazon" in google search. Currently monthly data are being provided. For each month from January, 2004 until the most recent month, google takes the number of searches (search volume) for "amazon" in that month and divides it by the total number of searches for all possible queries in that month. The resulting data (for all months) is then scaled to lie between 0 and 100. These datasets can be downloaded. 
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

#bivariate google trends
nndl.raw = read.delim("NNDLtrends25Aug2022.csv")
#the useful data is between row 2 and row 225
nndl.use = nndl.raw[2:225,1]
len = length(nndl.use)
nndl = matrix(0, len, 2)
for(i in 1:len)
{
  nndl[i, 1] = as.numeric(unlist(strsplit(as.character(nndl.use[i]), ","))[2])
  nndl[i, 2] = as.numeric(unlist(strsplit(as.character(nndl.use[i]), ","))[3])
}
dl.ts = ts(nndl[,2], start = c(2004, 1), end = c(2022, 8), frequency = 12)
nn.ts = ts(nndl[,1], start = c(2004, 1), end = c(2022, 8), frequency = 12)

plot(dl.ts, ylab = "Neural Network vs Deep Learning Popularity", xlab = "Month", main = "Google trends data for neural network/deep learning", col = "red", ylim = c(0, 100))
points(nn.ts, type = "l")

#Correlogram Illustration:
#Consider the unemployment rate dataset:
plot(unrate.ts)
acf(unrate.ts, lag.max = 40, type = "correlation", plot = T)

#Correlogram for data generated according to the white noise model:
gwn = rnorm(200, 0, 1)
plot(gwn, type = "l", xlab = "Time", main = "Gaussian White Noise")
acf(gwn, lag.max = 40, type = "correlation", plot = T)

#Is white noise a good model for each of the following datasets?
samp.dat = arima.sim(n = 200, list(ma = 0.7), sd = 1)
plot(samp.dat, type = "l", xlab = "Time")
acf(samp.dat, lag.max = 20, type = "correlation", plot = T)

samp.dat = arima.sim(n = 200, list(ma = -0.7), sd = 1)
plot(samp.dat, type = "l", xlab = "Time")
acf(samp.dat, lag.max = 20, type = "correlation", plot = T)

samp.dat = arima.sim(n = 200, list(ar = -0.7), sd = 1)
plot(samp.dat, type = "l", xlab = "Time")
acf(samp.dat, lag.max = 20, type = "correlation", plot = T)



