# Script for Assignment ATSA, Alessandro Salvatori r0926525

library(CADFtest)
library(forecast)
library(fGarch)
library(vars)
library(urca)

# Upload csv file
dataIT <- read.csv("GDP_IT.csv", header = TRUE, dec = ",")
summary(dataIT)
head(dataIT)

### UNIVARIATE TIME SERIES ANALYSIS ###

data_GDP_num <- as.numeric(dataIT$Value)
GDP_ts <- ts(data_GDP_num, frequency = 4, start = c(1996,1), end = c(2016,4)) 

# Using log time series because we are interested in learning the 
# relative differences 
logGDP_ts <- log(GDP_ts)
ts.plot(logGDP_ts, col="green")
acf(logGDP_ts)
pacf(logGDP_ts)
# The correlogram shows strong persistency.

# Testing for white noise and stationarity
max.lag <- round(sqrt(length(logGDP_ts)))
Box.test(logGDP_ts, lag = max.lag, type = "Ljung-Box")
CADFtest(logGDP_ts, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# The time series is not a white noise and the p-value obtained from the ADF test (0.136) is greater than 0.05, 
# so we do not reject H0 and conclude that the time-series is not stationary and we go to differences.

# Italian GDP in differences
dlogGDP_ts <- diff(logGDP_ts)
ts.plot(dlogGDP_ts, col="brown")
acf(dlogGDP_ts)
pacf(dlogGDP_ts)
# From the autocorrelations we can see that there might be stationarity as the 
# autocorrelations decline after lag 1, which suggests an MA(2). While, the 
# partial correlations drop after lag 1, so it might suggest an AR(1) model.

# Testing for white noise and stationarity
Box.test(dlogGDP_ts, lag = max.lag, type = "Ljung-Box")
CADFtest(dlogGDP_ts, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# The time series is not a white noise and the p-value obtained (0.001) is less than 0.05, so we reject H0 
# and conclude that the time-series is stationary.

# Checking for seasonality
monthplot(dlogGDP_ts)
# No seasonality

# Fitting first model: AR(1)
fit_ar <- arima(logGDP_ts, order = c(1,1,0), seasonal = list(order = c(0,0,0)))
summary(fit_ar)
par(mfrow=c(2,1))
ts.plot(fit_ar$residuals)
acf(fit_ar$residuals)
Box.test(fit_ar$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value (0.306) indicates that we do not reject H0 and we conclude that the
# residuals are white-noise, so the model is validated.

# Fitting second model: MA(2)
fit_ma2 <- arima(logGDP_ts, order = c(0,1,2), seasonal = list(order = c(0,0,0)))
summary(fit_ma2)
par(mfrow=c(2,1))
ts.plot(fit_ma2$residuals)
acf(fit_ma2$residuals)
Box.test(fit_ma2$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value (0.788) indicates that we do not reject H0 and we conclude that the
# residuals are white-noise, so the model is validated.

# Comparing the models (AR(1), MA(2)) using AIC and SIC
AIC(fit_ar); AIC(fit_ma2)
AIC(fit_ar,k=log(length(logGDP_ts))); AIC(fit_ma2,k=log(length(logGDP_ts)))
# The AIC measure indicates that the MA(2) model is slightly better, but the SIC measure indicates that AR(1) model is 
# slightly lower, hence better, so we prefer that. Also the AR(1) model is more parsimonious.
# Since results are still really close we still try to forecast using both models and then compare results.
par(mfrow=c(1,1))

### FORECASTING ###

# Forecast AR(1)
myforecastAR<-predict(fit_ar,n.ahead=8)
expected<-myforecastAR$pred
# The confidence bounds of the 95% prediction interval:
lower<-myforecastAR$pred-qnorm(0.975)*myforecastAR$se
upper<-myforecastAR$pred+qnorm(0.975)*myforecastAR$se
plot.ts(logGDP_ts,xlim=c(2005,2019),ylim=c(4.4,4.8))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")

# Forecast MA(2)
myforecastMA<-predict(fit_ma2,n.ahead=8)
expected<-myforecastMA$pred 
# The confidence bounds of the 95% prediction interval:
lower<-myforecastMA$pred-qnorm(0.975)*myforecastMA$se
upper<-myforecastMA$pred+qnorm(0.975)*myforecastMA$se
plot.ts(logGDP_ts,xlim=c(2005,2019),ylim=c(4.4,4.8))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")
#Forecast of the two models is pretty much the same.

# Now it is used an expanding-window approach to forecast values of logGDP_ts for 1-period.
# The loops are necessary to compare the forecasts for the ARIMA(1,1,0) model and the ARIMA(0,1,2) model.
y<-logGDP_ts
S=round(0.75*length(y))
h=1
error1.h <- c()
for (i in S:(length(y)-h))
  {
  mymodel.sub<-arima(y[1:i], order = c(1,1,0),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h<-c(error1.h,y[i+h]-predict.h)
}
error2.h <- c()
for (i in S:(length(y)-h))
  {
  mymodel.sub<-arima(y[1:i], order = c(0,1,2),seasonal=c(0,0,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}

# To evaluate the forecasting performance of the two models the Mean Absolute Error (MAE) is calculated:
MAE1 <- mean(abs(error1.h))
MAE2 <- mean(abs(error2.h))
MAE1; MAE2
# The AR(1) model presents a better (lower) Mean Absolute Error (0.00301)

# The same is now done using the squared value loss.
MSE1 <- mean(abs(error1.h^2))
MSE2 <- mean(abs(error2.h^2))
MSE1; MSE2
# The AR(1) obtains the lowest MSE (1.410299e-05).

# Now a Diebold-Mariano test is performed, to see if the forecast performance of the two models is significantly different
dm.test(error1.h,error2.h,h=h,power=1)
dm.test(error1.h,error2.h,h=h,power=2)
# Diebold-Mariano test:
# P-value = 0.065 > 5%, thus H0 is not rejected and it can be concluded  that the forecast performance
# of the two models, using the squared value loss, is not significantly different.

# Checking for heteroscedasticity
acf(fit_ar$residuals^2)
Box.test(fit_ar$residuals^2,lag=max.lag,type="Ljung-Box")
# There are some significant autocorrelations in the squared residuals. We might be in presence
# of heteroscedasticity so we decide to fit a garch model

# Garch Model
fit_garch<-garchFit(~arma(1,0)+garch(1,1),data=logGDP_ts)
summary(fit_garch)
plot(fit_garch)
# The the Jarque-Bera and the Shapiro-Wilk tests have p-value = 0.264 < 5% and p-value = 0.0955 < 5%,
# respectively. In both tests, we do not reject H0 and conclude that the standardized residuals are normal.
# From the acf of the standardized residuals and of the squared standardized residuals
# we observe that there are 2 slightly significant autocorrelation in the first one and 0 in the second one
# and conclude that the model could be valid.
# The Q-tests on the standardized residuals and on the squared standardized residuals have
# p-values > 5% for Q(10), Q(15) and Q(20) (only for standardized residuals Q(10) we have a 0.046). 
# Thus, we conclude that there is no structure left in the (squared) standardized residuals. Hence, the model is valid.
# From the conditional standard deviation we observe that the standard deviation is not constant over time 
# and that there are clusters of high volatility.


### MULTIVARIATE TIME SERIES ANALAYSIS ###

# Upload csv file
dataEU <- read.csv("GDP_EU.csv", header = TRUE, dec = ",")
summary(dataEU)
head(dataEU)

data_EU_GDP_num <- as.numeric(dataEU$Value)
GDP_EU_ts <- ts(data_EU_GDP_num, frequency = 4, start = c(1996,1), end = c(2016,4)) 

# Using log time series because we are interested in learning the 
# relative differences 
logGDP_EU_ts <- log(GDP_EU_ts)
ts.plot(logGDP_EU_ts)
ts.plot(logGDP_ts, logGDP_EU_ts, col=(c("red","blue")))
ts.plot(logGDP_ts-logGDP_EU_ts)
acf(logGDP_EU_ts)
pacf(logGDP_EU_ts)
# The correlogram shows strong persistency.
# The partial correlations drop after lag 1.

# Testing for white noise and stationarity
library(CADFtest)
max.lag <- round(sqrt(length(logGDP_EU_ts)))
Box.test(logGDP_EU_ts, lag = max.lag, type = "Ljung-Box")
CADFtest(logGDP_EU_ts, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# The time series is not white noise and the p-value obtained from the ADF test (0.495) is greater than 0.05, so we do not reject H0 and conclude 
# that the time-series is not stationary and we go to differences.

#EU GDP in differences
dlogGDP_EU_ts <- diff(logGDP_EU_ts)
ts.plot(dlogGDP_EU_ts, col="brown")
acf(dlogGDP_EU_ts)
pacf(dlogGDP_EU_ts)

# Testing for autocorrelation and stationarity (white noise, box test)
Box.test(dlogGDP_EU_ts, lag = max.lag, type = "Ljung-Box")
CADFtest(dlogGDP_EU_ts, type = "drift", criterion = "BIC", max.lag.y = max.lag)
# The time series is not white noise and the p-value obtained (0.003) is less than 0.05, so we reject H0 and conclude 
#that the time-series is stationary.


# Multivariate TSA: ADLM(1)

# Estimation of an Autoregressive Dynamic model of order 1 for dlogGDP_ts (ADLM(1))
lag <- 1
n <- length(dlogGDP_ts)
dlogGDP.0 <- dlogGDP_ts[(lag+1):n]
dlogGDP.1 <- dlogGDP_ts[lag:(n-1)]
dlogGDP_EU.1 <- dlogGDP_EU_ts[lag:(n-1)]
fit_adlm <- lm(dlogGDP.0 ~ dlogGDP.1+dlogGDP_EU.1)
acf(fit_adlm$residuals)
sqrt(length(dlogGDP_ts))
Box.test(fit_adlm$residuals, lag = max.lag, type = "Ljung-Box")
summary(fit_adlm)
# R-squared:  0.3197, 31.97% of the variance of dlogGDP.0 is explained by the model. 
# The overall F-statistics has p-value (2.469e-07) < 5%, it is possible to conclude that the regressors are jointly significant.  
# The only significant variable is dlogGDP_EU.1, that has a positive value.

# Testing for Granger causality comparing the ADLM(1) with the model without lagged explanatory variables:
fit_adlm_nox <- lm(dlogGDP.0 ~ dlogGDP.1)
anova(fit_adlm,fit_adlm_nox)
# The p-value = 0.0326 < 5%, thus we reject H0 of no Granger Causality. We conclude that
# dlogGDP_EU has incremental explanatory power in predicting dlogGDP.


# Multivariate TSA: VAR Model

data<-data.frame(dlogGDP_ts, dlogGDP_EU_ts)
names(data) <- c("dlogGDP","dlogGDP_EU")
VARselect(data,lag.max=10,type="const")
# The maximum number of lags is specified to be 10, a constant is included in each equation.
# The order selected by the Schwarz information criterion is 1.

# Estimating the VAR model
fit_varautom <- VAR(data,type="const",p=1)
summary(fit_varautom)
varautom_residuals<-resid(fit_varautom)
acf(varautom_residuals[,1])
acf(varautom_residuals[,2])
ccf(varautom_residuals[,1], varautom_residuals[,2])
Box.test(varautom_residuals[,1], lag = max.lag, type="Ljung-Box")
Box.test(varautom_residuals[,2], lag = max.lag, type="Ljung-Box")
# The R2 = 0.4453, thus 44.53% of the variance of dlogGDP is explained by the lagged observations
# of dlogGDP and of dlogGDP_EU at lag 1.
# The F-statistics has p-value (7.789e-11) < 5%, thus we reject H0 and conclude that the regressors
# are jointly significant.
# In the correlograms we observe no significant correlations, while a strong significant cross-correlation 
# without lag is observed in the cross-correlogram. There is no problem with contemporaneous correlations 
# in a VAR model, as it does not allow for explicit modeling of contemporaneous interdependence. 
# Thus, the residuals look multivariate white noise and the model is valid.

# Impulse response function
irf_var <- irf(fit_varautom, ortho=F, boot=T)
plot(irf_var)
# The IRFs provide an easy way to interpret the estimated coefficients of the VAR model. 
# Given a unitary impulse in dlogGDP at time t, we observe a positive response of dlogGDP_EU at t + 1.
# Given a unitary impulse in dlogGDP_EU at time t, we observe a positive response of dlogGDP at time t + 1 
# but not significant because the 0 is contained in the prediction interval.

# Testing for cointegration with Engle-Granger Test
fit_ci <- lm(logGDP_ts ~ logGDP_EU_ts)
res_fit_ci <- fit_ci$residuals
CADFtest(res_fit_ci,type="drift",criterion="BIC",max.lag.y=max.lag)
# Test statistics: -1.12, which is bigger than the Engle-Granger ADF test statistics for one explanatory variable −3.41. 
# Thus, H0 of no cointegration is not rejected and it is possible to conclude that logGDP_ts and logGDP_EU_ts are not cointegrated.


data2<-data.frame(logGDP_ts, logGDP_EU_ts)
#Testing for cointegration with Johansen Test
trace_test <- ca.jo(data2,type="trace",K=4,ecdet="const",spec="transitory")
summary(trace_test)
# For r=0, the test statistics is larger than then the critical value (29.16 > 19.96), thus there
# is at least one cointegrating relation.
# The cointegrating equation is 0.002 + log(GDPt) − 0.922 log(GDP_EUt) = δt
maxeigen_test <- ca.jo(data2,type="eigen",K=4,ecdet="const",spec="transitory")
summary(maxeigen_test)
# After repeating the cointegration test using the Johansen’s maximum eigenvalue statistics
# we conclude that logGDP and logGDP_EU are cointegrated.


# Multivariate TSA: Vector Error Correction Model
fit_vecm1 <- cajorls(trace_test,r=1)
fit_vecm1
fit_vecm2 <- cajorls(maxeigen_test,r=1)
fit_vecm2
# The cointegrating equation, which corresponds to the one while testing with the johansen procedure, is a stationary
# linear combination of logGDP and logGDP_EU.

fit_var<-vec2var(trace_test,r=1)
myforecast<-predict(fit_var,n.ahead=8)
par(mfrow=c(2,1))
logGDP_ts_forecast<-ts(myforecast$fcst$logGDP_ts[,1],frequency=4,start=c(2017,1))
logGDP_ts_lower<-ts(myforecast$fcst$logGDP_ts[,2],frequency=4,start=c(2017,1))
logGDP_ts_upper<-ts(myforecast$fcst$logGDP_ts[,3],frequency=4,start=c(2017,1))
ts.plot(logGDP_ts_forecast,logGDP_ts_lower,logGDP_ts_upper,col=c("blue","red","red"))
title(main = "8-step-ahead forecast of logGDP_EU_ts")

logGDP_EU_ts_forecast<-ts(myforecast$fcst$logGDP_EU_ts[,1],frequency=4,start=c(2017,2))
logGDP_EU_ts_lower<-ts(myforecast$fcst$logGDP_EU_ts[,2],frequency=4,start=c(2017,2))
logGDP_EU_ts_upper<-ts(myforecast$fcst$logGDP_EU_ts[,3],frequency=4,start=c(2017,2))
ts.plot(logGDP_EU_ts_forecast,logGDP_EU_ts_lower,logGDP_EU_ts_upper,col=c("blue","red","red"))
title(main = "8-step-ahead forecast of logGDP_EU_ts")

