model_31 <- arima(TGlobal_coupled[1302:1801],c(3,0,1))
model_1 <- arima(TGlobal_coupled[1302:1801],c(1,0,0))
model_iid <- arima(TGlobal_coupled[1302:1801],c(0,0,0))

sd_ratios <- data.frame(year=NA,model=NA,ratio=NA,lowerBound=NA,upperBound=NA)

#original years used in paper. 10-50 times span
for(j in 10:50){
  ensemble_slopes <- rep(NA,40)
  for(i in 1:40){
    ensemble_slopes[i]<- lm(TGlobal_runs[(142-j):141,i]~seq(1:j))$coefficients[2]
  }

  sim_slopes <-rep(NA,1000)
  for(i in 1:1000){
    sim_slopes[i] <- lm(arima.sim(list(ar=model_31$coef[1:3], ma=model_31$coef[4]), sd=sqrt(model_31$sigma2), n=j)~seq(1:j))$coefficients[2]
  }

  arma31Boot <- rep(NA,1000)
  for(i in 1:1000){
    arma31Boot[i] <-sd(sim_slopes)/sd(sample(ensemble_slopes,40, replace=T))
  }

  sd_ratios<- rbind(sd_ratios, c(j,"arma31",sd(sim_slopes)/sd(ensemble_slopes),quantile(arma31Boot,c(.025)),quantile(arma31Boot,c(.975))))
  

  sim_ar1_slopes <-rep(NA,1000)
  for(i in 1:1000){
    sim_ar1_slopes[i] <- lm(arima.sim(list(ar=model_1$coef[1]), sd=sqrt(model_1$sigma2), n=j)~seq(1:j))$coefficients[2]
  }


  arma1Boot <- rep(NA,1000)
  for(i in 1:1000){
    arma1Boot[i] <-sd(sim_ar1_slopes)/sd(sample(ensemble_slopes,40, replace=T))
  }

  sd_ratios<-rbind(sd_ratios, c(j,"ar1",sd(sim_ar1_slopes)/sd(ensemble_slopes),quantile(arma1Boot,c(.025)),quantile(arma1Boot,c(.975))))
  
  #iid stuff, uncomment if we care about it
  sim_slopes_iid <-rep(NA,1000)
  # for(i in 1:1000){
  #   sim_slopes_iid[i] <- lm(rnorm(n=j,sd=sqrt(model_iid$sigma2))~seq(1:j))$coefficients[2]
  # }
  # sd_ratios[j-9,4] <- sd(sim_slopes_iid)/sd(ensemble_slopes)
}

library(ggplot2)
library(tidyr)
#sd_ratios_long <- gather(sd_ratios, "model", "ratio", ar1:arma31, "lowerbound")
sd_ratios <- sd_ratios[complete.cases(sd_ratios),]
sd_ratios$ratio<- as.numeric(sd_ratios$ratio)
sd_ratios$year<- as.numeric(sd_ratios$year)
limits <- aes(ymax = as.numeric(upperBound), ymin=as.numeric(lowerBound))

p <- ggplot(sd_ratios,aes(x=year,col=model,y=ratio))
p+geom_point()+geom_line()+geom_errorbar(limits)


##exact same code but can be tweeked to use different end years
##runs faster since only every 3rd year

sd_ratios <- data.frame(year=NA,model=NA,ratio=NA,lowerBound=NA,upperBound=NA)

for(j in 3:17){
  j <- 3*j
  ensemble_slopes <- rep(NA,40)
  for(i in 1:40){
    ensemble_slopes[i]<- lm(TGlobal_runs[(162-j):161,i]~seq(1:j))$coefficients[2]
  }
  
  sim_slopes <-rep(NA,1000)
  for(i in 1:1000){
    sim_slopes[i] <- lm(arima.sim(list(ar=model_31$coef[1:3], ma=model_31$coef[4]), sd=sqrt(model_31$sigma2), n=j)~seq(1:j))$coefficients[2]
  }
  
  arma31Boot <- rep(NA,1000)
  for(i in 1:1000){
    arma31Boot[i] <-sd(sim_slopes)/sd(sample(ensemble_slopes,40, replace=T))
  }
  
  sd_ratios<- rbind(sd_ratios, c(j,"arma31",sd(sim_slopes)/sd(ensemble_slopes),quantile(arma31Boot,c(.025)),quantile(arma31Boot,c(.975))))
  
  
  sim_ar1_slopes <-rep(NA,1000)
  for(i in 1:1000){
    sim_ar1_slopes[i] <- lm(arima.sim(list(ar=model_1$coef[1]), sd=sqrt(model_1$sigma2), n=j)~seq(1:j))$coefficients[2]
  }
  
  
  arma1Boot <- rep(NA,1000)
  for(i in 1:1000){
    arma1Boot[i] <-sd(sim_ar1_slopes)/sd(sample(ensemble_slopes,40, replace=T))
  }
  
  sd_ratios<-rbind(sd_ratios, c(j,"ar1",sd(sim_ar1_slopes)/sd(ensemble_slopes),quantile(arma1Boot,c(.025)),quantile(arma1Boot,c(.975))))
  
  sim_slopes_iid <-rep(NA,1000)
  # for(i in 1:1000){
  #   sim_slopes_iid[i] <- lm(rnorm(n=j,sd=sqrt(model_iid$sigma2))~seq(1:j))$coefficients[2]
  # }
  # sd_ratios[j-9,4] <- sd(sim_slopes_iid)/sd(ensemble_slopes)
}

library(ggplot2)
library(tidyr)
#sd_ratios_long <- gather(sd_ratios, "model", "ratio", ar1:arma31, "lowerbound")
sd_ratios <- sd_ratios[complete.cases(sd_ratios),]
sd_ratios$ratio<- as.numeric(sd_ratios$ratio)
sd_ratios$year<- as.numeric(sd_ratios$year)
limits <- aes(ymax = as.numeric(upperBound), ymin=as.numeric(lowerBound))

p <- ggplot(sd_ratios,aes(x=year,col=model,y=ratio))
p+geom_point()+geom_line()+geom_errorbar(limits)

