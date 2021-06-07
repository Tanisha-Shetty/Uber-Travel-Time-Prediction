library(zoo)
library(fma)
library(fpp)
library(forecast)

# Read the data in the CSV file into a zoo object
# The date is taken as the index with the provided format %Y-%m-%d
uberLondon <- read.zoo(file = 'uberLondon.csv',sep=" ", header=TRUE, format = "%Y-%m-%d")

# Q1)
# Split the data into train and test sets
# Training set starts on 02-01-2016 and ends on 30-11-2017
uberTrain <- window(uberLondon,end = as.Date("2017-11-30"))
# Test set starts on 01-12-2017 and ends on 31-12-2017
uberTest <- window(uberLondon,start = as.Date("2017-12-01"))
#There are 31 observations in the test data set
nrow(uberTest)

# Q2)
# Check the summary statistics of the training set.
summary(uberTrain)

### Figure 1 in the Report ###
# Plot of the training set
plot(uberTrain, main="Training Set details", xlab="Time(Years)")

# ACF of all variables
par(mfrow=c(2,3))
acf_tmin <- acf(uberTrain$tmin, lag.max=1000, na.action = na.pass, plot=FALSE)
plot(acf_tmin, main="ACF-Minimum Temperature")
acf_tmax <- acf(uberTrain$tmax, lag.max=1000, na.action = na.pass, plot=FALSE)
plot(acf_tmax, main="ACF-Maximum Temperature")
acf_tavg <- acf(uberTrain$tavg, lag.max=1000, na.action = na.pass, plot=FALSE)
plot(acf_tavg, main="ACF-Average Temperature")
acf_prcp <- acf(uberTrain$prcp, lag.max=1000, na.action = na.pass, plot=FALSE)
plot(acf_prcp, main="ACF-Precipitation")
acf_traveltime <- acf(uberTrain$MeanTravelTimeSeconds, na.action = na.pass, lag.max=50, plot=FALSE)
plot(acf_traveltime, main="ACF-Mean Travel Time")
acf_windspeed <- acf(uberTrain$windspeed, na.action = na.pass, lag.max=50, plot=FALSE)
plot(acf_windspeed, main="ACF-Windspeed")

### Figure 2 in the Report ###
# Correlation matrix of all variables
cor(uberTrain)

# Check the outliers of all variables
outliersTravelTime <- boxplot.stats(as.vector(uberTrain$MeanTravelTimeSeconds))$out
outliersTmax <- boxplot.stats(as.vector(uberTrain$tmax))$out
outliersTmin <- boxplot.stats(as.vector(uberTrain$tmin))$out
outliersTavg <- boxplot.stats(as.vector(uberTrain$tavg))$out
outliersPrcp <- boxplot.stats(as.vector(uberTrain$prcp))$out
outliersWindspeed <- boxplot.stats(as.vector(uberTrain$windspeed))$out
length(outliersTravelTime)
length(outliersTmax)
length(outliersTmin)
length(outliersTavg)
length(outliersPrcp)
length(outliersWindspeed)

# Investigate the NA values in the dataset
colSums(is.na(uberLondon))
uberLondon [ !complete.cases(uberLondon), ]

# Q3)
# We first fix the NA values using linear interpolation
# in the original dataset and split to train and test datasets
# We can interpolate the missing values
uberLondon$tmin <- na.interp(uberLondon$tmin)
uberLondon$tmax <- na.interp(uberLondon$tmax)
uberLondon$prcp <- na.interp(uberLondon$prcp)

uberTrain <- window(uberLondon,end = as.Date("2017-11-30"))
uberTest <- window(uberLondon,start = as.Date("2017-12-01"))

# Check for stationarity
Box.test(uberTrain$MeanTravelTimeSeconds, type="Ljung-Box")
kpss.test(uberTrain$MeanTravelTimeSeconds)
adf.test(uberTrain$MeanTravelTimeSeconds)
# Data appears to be non-stationary

# Analyse the ACF of the mean travel times
acf(coredata(uberTrain$MeanTravelTimeSeconds), lag.max=100)

# There are seasonal variations in the ACF every 7 lags
# We can set this frequency as 7 for the time series of mean travel times
travel_times_train <- ts(uberTrain$MeanTravelTimeSeconds, frequency=7)
travel_times_test <- ts(uberTest$MeanTravelTimeSeconds, start = tail(index(travel_times_train), n=1), frequency=7)

# Decompose the time series of mean travel times
decomposed_train <- decompose(travel_times_train)
plot(decomposed_train)

# There appears to be seasonality and a complex trend in the decomposed series
travel_sdiff <- diff(travel_times_train,7)

# On seasonal differencing the series as above(lag 7), we can see white noise in the acf
par(mfrow=c(1,1))
acf(coredata(travel_sdiff), lag.max=20)
pacf(coredata(travel_sdiff), lag.max = 20)

# The ACF appears to cut off at lag 5. This indicates a non-seasonal MA(5) component
# The significant spike at lag 7 indicates a seasonal MA(1) component
# The PACF appears to cut off at lag 3
# There are significant spikes at every 7th and 8th lag on the PACF
# These spikes show exponential decay in successive seasonal lags
# The significant spikes at lag 7 and 8 indicate that we need seasonal AR(1) component
# We can try an ARIMA(3,0,5)(1,1,1)[7]

# Arima(3,0,5)(1,1,1)[7]
# AIC=7817.91   AICc=7818.29   BIC=7867.84
fit<-Arima(travel_times_train,c(3,0,5),c(1,1,1))

# Auto.arima for comparision
# ARIMA(0,0,4)(0,1,1)[7] 
# AIC=7902.6   AICc=7902.72   BIC=7929.83
fit_auto<-auto.arima(travel_times_train, seasonal=TRUE, stepwise = FALSE, approximation=FALSE)

# Similar ARIMA models with various combinations of regression errors:

# Arima(3,0,5)(1,1,1)[7] with tavg as predictors
# AIC=7819.78   AICc=7820.24   BIC=7874.25
fit_tavg <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=uberTrain$tavg)

# Arima(3,0,5)(1,1,1)[7] with tmin as predictors
# AIC=7819.59   AICc=7820.05   BIC=7874.07
fit_tmax <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=uberTrain$tmax)

# Arima(3,0,5)(1,1,1)[7] with tmax as predictors
fit_tmin <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=uberTrain$tmin)
# AIC=7819.89   AICc=7820.35   BIC=7874.36

# Arima(3,0,5)(1,1,1)[7] with prcp as predictors
# AIC=7817.18   AICc=7817.64   BIC=7871.66
fit_prcp <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=uberTrain$prcp)

# Arima(3,0,5)(1,1,1)[7] with windspeed as predictors
# AIC=7819.31   AICc=7819.77   BIC=7873.79
fit_windspeed <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=uberTrain$windspeed)

# tmin and tavg showed higher correlation with mean travel times in Q2
# We try a combination of the two to see if the prediction is better
covariates_temp_train <- cbind(tmin=uberTrain$tmin,
                               tavg=uberTrain$tavg)

# Arima(3,0,5)(1,1,1)[7] with tmin,tavg as predictors
# AIC=7821.4   AICc=7821.94   BIC=7880.41
fit_multixreg_temp <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=covariates_temp_train)


# We also check using all the other variables as predictors
covariates_all_train <- cbind(tmax=uberTrain$tmax,
                              tmin=uberTrain$tmin,
                              windspeed=uberTrain$windspeed,
                              tavg = uberTrain$tavg,
                              prcp = uberTrain$prcp)

# Arima(3,0,5)(1,1,1)[7] with all other weather data as predictors
# AIC=7824.13   AICc=7824.94   BIC=7896.77
fit_multixreg_all <- Arima(travel_times_train,c(3,0,5),c(1,1,1), xreg=covariates_all_train)

# The AIC and BIC are better for our model compared to the ARIMA(0,0,4)(0,1,1)[7] model chosen by auto.arima 
# The regressors don't have a significant effect in improving the SARIMA model fit

# In the decomposed series(the variable decomposed_train), we saw a complex trend
# This may also indicate non-constant variance
# We try the damped Holt-Winters model on the time series
forecast_holt <- hw(travel_times_train, damped = TRUE, seasonal="multiplicative", h=31)

# Neural network models based on the predictors chosen for the SARIMA models

# Neural Network without predictors
set.seed(12345)
nn_fit = nnetar(travel_times_train)

# Neural Network with the tavg and tmin as predictors
set.seed(12345)
nn_tmin_tmax  = nnetar(travel_times_train, xreg = covariates_temp_train)

# Neural Network with the all variables as predictors
set.seed(12345)
nn_fit_all  = nnetar(travel_times_train, xreg = covariates_all_train)

# We compute the RMSE scores of the SARIMA with the best AIC and BIC, the Holt Winters model
# and all neural network models.

covariates_temp_test <- cbind(tmin=uberTest$tmin,
                             tavg = uberTest$tavg)

covariates_all_test <- cbind(tmax=uberTest$tmax,
                            tmin=uberTest$tmin,
                            windspeed=uberTest$windspeed,
                            tavg = uberTest$tavg,
                            prcp = uberTest$prcp)

# Arima(3,0,5)(1,1,1)[7]
forecast_sarima <- forecast(fit, h=31)
# Arima(3,0,5)(1,1,1)[7] with tavg and tmin as predictors
forecast_sarima_temp <- forecast(fit_multixreg_temp, h=31, xreg=covariates_temp_test)
# Arima(3,0,5)(1,1,1)[7] with all weather data as predictors
forecast_sarima_all <- forecast(fit_multixreg_all, h=31, xreg=covariates_all_test)
# NNAR(22,1,12)[7]
forecast_nn <- forecast(nn_fit,PI=TRUE, h=31)
# NNAR(22,1,12)[7] with tavg and tmin as predictors
forecast_nn_temp <- forecast(nn_tmin_tmax,PI=TRUE, h=31, xreg=covariates_temp_test)
# NNAR(22,1,14)[7] with all weather data as predictors 
forecast_nn_all <- forecast(nn_fit_all,PI=TRUE, h=31, xreg=covariates_all_test)

rmse_sarima<-rmse(as.ts(travel_times_test), forecast_sarima$mean)
rmse_sarima_temp<-rmse(as.ts(travel_times_test), forecast_sarima_temp$mean)
rmse_sarima_all<-rmse(as.ts(travel_times_test), forecast_sarima_all$mean)
rmse_holt <- rmse(as.ts(travel_times_test),forecast_holt$mean)
rmse_nn <- rmse(as.ts(travel_times_test),forecast_nn$mean)
rmse_nn_temp <- rmse(as.ts(travel_times_test),forecast_nn_temp$mean)
rmse_nn_all <- rmse(as.ts(travel_times_test),forecast_nn_all$mean)
# nn_fit_all appears to be the best fit from these models since it 
# has the lowest RMSE score of  165.2068
# This is a NNAR(22,1,14)[7] model with P=22, p=1 and k=(P+p+1)/2=14

### Figure 3 in the Report ###
# The chosen model is NNAR(22,1,14)[7] with weather data as predictors. We now check the residuals of the fit
checkresiduals(nn_fit_all)
# Check for stationarity
Box.test(nn_fit_all$residuals, lag = 10, type="Ljung-Box")

### Figure 4 in the Report ###
# Plot of forecast of our chosen model - NNAR(22,1,14)[7]
plot(forecast_nn_all,main="Forecast of NNAR(22,1,14)[7]",xlab="Time", ylab="Trip Time(Seconds)")
lines(forecast_nn_all$fitted, col='red')
lines(travel_times_test,col="black")
legend("topleft", legend=c("Forecast","Actual", "Fitted"),
       col=c("blue", "black", "red"), lty=1)
