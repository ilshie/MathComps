#New data 
TGlobal_runs <- TGlobal
plot.ts(TGlobal_runs[,3])
spec(TGlobal_runs)

#Fit model using last 500 yrs of coupled control
length(TGlobal_c)
coupled_500 <- TGlobal_c[1301:1801]
acf(coupled_500, lag.max = 35)
pacf(coupled_500, lag.max = 35) 
plot.ts(coupled_500)
arima_auto_c <- auto.arima(coupled_500, max.d=0,max.D = 0, max.Q = 0,stepwise = F, ic = "aic")
arima_auto_c

#Start with AR(3)
coupled_5_1 <- arima(coupled_500, order = c(3,0,0))
plot(rstandard(coupled_5_1),ylab ='Standardized Residuals',type='o'); abline(h=0)
qqnorm(residuals(coupled_5_1)); qqline(residuals(coupled_5_1))
acf(residuals(coupled_5_1))
tsdiag(coupled_5_1, gof = 15, omit.initial = F)

#add ma term
coupled_5_2 <- arima(coupled_500, order = c(3,0,1))
plot(rstandard(coupled_5_2),ylab ='Standardized Residuals',type='o'); abline(h=0)
qqnorm(residuals(coupled_5_2)); qqline(residuals(coupled_5_2))
acf(residuals(coupled_5_2))
tsdiag(coupled_5_2, gof = 15, omit.initial = F)


#compare spectral density with theoretical spectral density for ARMA(3,1)
sp_c_5 <- spec(coupled_500 , log = "no", xlab = "Frequency", ylab = "Sample Spectral Density", sub = "")
k = kernel('daniell', m = 5)
sp_c_5_2 <- spec(coupled_500, kernel = k, log = 'no', sub = '', xlab='Frequency',ylab='Smoothed Sample Spectral Density')
arma.spec.obj_c_5 <- arma.spec(ar = c(1.2881, -0.7897, 0.3592), ma = c(-0.6173), var.noise = 0.006226)
lines(arma.spec.obj_c_5$freq, arma.spec.obj_c_5$spec)

#Keep ARMA(3,1) for now
#Simulate 50 year time series from model, fit linear trend, store slope coeff
#0.00182
n <- 1000
coeffs <- c()
for (i in (1:n)){
  series <- arima.sim(model = list(ar = c(1.2881, -0.7897, 0.3592), ma = c(-0.6173)), sd = sqrt(0.006226), n = 50)
  index <- c(1:50)
  linearcoef <- lm(series ~ index)$coef[2]
  coeffs[i] <- linearcoef
}
summary(coeffs)
sd(coeffs)


#----------------------------------------------------------------------------
#compare spectra of AR(1) and ARMA(3,1)
arma.spec.obj_c_5_2 <- arma.spec(ar =c(0.5064), var.noise = 0.007198)
arma.spec.obj_c_5 <- arma.spec(ar = c(1.2881, -0.7897, 0.3592), ma = c(-0.6173), var.noise = 0.006226)
lines(arma.spec.obj_c_5_2$freq, arma.spec.obj_c_5_2$spec)

#------------------------------------------------------------------
#Fit a linear trend to each of the 40 ensemble members (2011-2061)
#0.00156
#Find 2011-2016
Runs_50 <- TGlobal_runs[92:141, ]
length(Runs_50[ ,1])

coeffs_ens <- c()

for (i in (1:40)){
  index <- c(1:50)
  linearcoef <- lm(Runs_50[,i] ~ index)$coef[2]
  coeffs_ens[i] <- linearcoef
}
summary(coeffs_ens)
sd(coeffs_ens)


#__________________________________________________________________
#Compare to AR(1) model - not great fit for data but estimate close to sd of ensemble members
#0.001588159
ar1_coupled <- arima(coupled_500, order = c(1,0,0))
ar1_coupled
plot(rstandard(ar1_coupled),ylab ='Standardized Residuals',type='o'); abline(h=0)
qqnorm(residuals(ar1_coupled)); qqline(residuals(ar1_coupled))
acf(residuals(ar1_coupled))
tsdiag(ar1_coupled, gof = 15, omit.initial = F)

n <- 1000
coeffs <- c()
for (i in (1:n)){
  series <- arima.sim(model = list(ar = c(0.5064)),sd = sqrt(0.007198), n = 50)  
  index <- c(1:50)
  linearcoef <- lm(series ~ index)$coef[2]
  coeffs[i] <- linearcoef
}
summary(coeffs)
sd(coeffs)

#--------------------------------------------------------------------------------
#Compared to ARMA(3,1) model - seems to overestimate sd slightly 
#0.001681623
ar31_coupled <- arima(coupled_500, order = c(3,0,1))
ar31_coupled

plot(rstandard(ar31_coupled),ylab ='Standardized Residuals',type='o'); abline(h=0)
qqnorm(residuals(ar31_coupled)); qqline(residuals(ar31_coupled))
acf(residuals(ar31_coupled))
tsdiag(ar31_coupled, gof = 15, omit.initial = F)
  
n <- 1000
coeffs31 <- c()
for (i in (1:n)){
  series <- arima.sim(model = list(ar = c(1.2881,-0.7897, 0.3592), ma = c(-0.6173)),sd = sqrt(0.006226), n = 50)  
  index <- c(1:50)
  linearcoef <- lm(series ~ index)$coef[2]
  coeffs31[i] <- linearcoef
}
summary(coeffs31)
sd(coeffs31)

  
#----------------------------------------------------------------------------------------------
#REPEAT FOR 10 YEAR CHUNKS
#Why would AR(1) overestimate the sd? Shouldn't it underestimate it if our explanation was correct?

#Linear trends for ensemble members: sd =  0.01241956
Runs_10 <- TGlobal_runs[132:141, ]
length(Runs_10[ ,1])

coeffs_ens_10 <- c()

for (i in (1:40)){
  index <- c(1:10)
  linearcoef <- lm(Runs_10[,i] ~ index)$coef[2]
  coeffs_ens_10[i] <- linearcoef
}
summary(coeffs_ens_10)
sd(coeffs_ens_10)


#AR(1) 10 YEARS - 0.01464884
n <- 1000
coeffs_10 <- c()
for (i in (1:n)){
  series <- arima.sim(model = list(ar = c(0.5064)),sd = sqrt(0.007198), n = 10)  
  index <- c(1:10)
  linearcoef <- lm(series ~ index)$coef[2]
  coeffs_10[i] <- linearcoef
}
summary(coeffs_10)
sd(coeffs_10)

#AR(3,1) 10 YEARS - 0.01264501
n <- 1000
coeffs31_10 <- c()
for (i in (1:n)){
  series <- arima.sim(model = list(ar = c(1.2881,-0.7897, 0.3592), ma = c(-0.6173)),sd = sqrt(0.006226), n = 10)  
  index <- c(1:10)
  linearcoef <- lm(series ~ index)$coef[2]
  coeffs31_10[i] <- linearcoef
}
summary(coeffs31_10)
sd(coeffs31_10)


#-------------------------------------------------------------------------------------------------------------
#Write loop to calculate percentage of sd that AR(1) and ARMA(3,1) simulations find for different time lengths

sd_31 <- c()
sd_1 <- c()
sd_ensemble <- c()
for (j in (10:50)){
    coeffs31_j <- c()
    coeffs1_j <- c()
  for (i in (1:n)){
    #AR(1) model
    series1 <- arima.sim(model = list(ar = c(1.2881,-0.7897, 0.3592), ma = c(-0.6173)),sd = sqrt(0.006226), n = j)  
    index <- c(1:j)
    linearcoef <- lm(series1 ~ index)$coef[2]
    coeffs31_j[i] <- linearcoef
    
    #ARM(3,1) model 
    series2 <- arima.sim(model = list(ar = c(0.5064)),sd = sqrt(0.007198), n = j)  
    index <- c(1:j)
    linearcoef <- lm(series2 ~ index)$coef[2]
    coeffs1_j[i] <- linearcoef
  }
  sd_31[j-10] <- sd(coeffs31_j)
  sd_1[j-10] <- sd(coeffs1_j)
  
  #Ensemble members
  Runs_j <- TGlobal_runs[(132 - (j-9) + 1):141, ]
  coeffs_ens_j <- c()
  
  for (k in (1:40)){
    index <- c(1:j)
    linearcoef <- lm(Runs_j[,k] ~ index)$coef[2]
    coeffs_ens_j[k] <- linearcoef
  }
  
  sd_ensemble[j-10] <- sd(coeffs_ens_j)
  
}

plot(sd_31)
lines(sd_1)
plot(sd_ensemble)
summary(sd_ensemble)
summary(sd_31)

#plot percentages
#Why is there are spike at 30???******
percent_1 <- (sd_1/sd_ensemble)*100
percent_31 <- (sd_31/sd_ensemble)*100

plot(percent_1, type = "l", col = 'blue', ylim = c(80, 150))
lines(percent_31, type = "l", col = "red")
abline(100,0)
lines(percent_1)
#---------------------------------------------------------------------------------------------


#Just Ensemble members
sd_ensemble <- c()
for (j in (10:50)){
    Runs_j <- TGlobal_runs[(132 - (j-9) + 1):141, ]
    coeffs_ens_j <- c()

  for (k in (1:40)){
    index <- c(1:j)
    linearcoef <- lm(Runs_j[,k] ~ index)$coef[2]
    coeffs_ens_j[k] <- linearcoef
    print(coeffs_ens_j)
}

sd_ensemble[j-10] <- sd(coeffs_ens_j)

}


#-------------------------------------------------------------------------------------
#Pulling out one of ensembles members at time and comparing to average over 40


