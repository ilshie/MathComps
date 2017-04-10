
#get model parameters from 500 years of coupled data
model_31 <- arima(TGlobal_coupled[1302:1801],c(3,0,1))
model_1 <- arima(TGlobal_coupled[1302:1801],c(1,0,0))

#set up dataframe to store ratios
sd_ratios <- data.frame(span=NA,model=NA,ratio=NA,lowerBound=NA,upperBound=NA)

#consider spans of 10 to 40 years
for(span in 10:50){
  print(span)
  
  #save ensemble variance at each end year
  varSlopeEndYears <- rep(NA,40)
  
  #get simulated arma3,1 sd
  sim_slopes_31 <-rep(NA,1000)
  for(i in 1:1000){
    sim_slopes_31[i] <- lm(arima.sim(list(ar=model_31$coef[1:3], ma=model_31$coef[4]), sd=sqrt(model_31$sigma2), n=span)~seq(1:span))$coefficients[2]
  }
  sd31 <- sd(sim_slopes_31)

  #get simulated ar1 sd
  sim_slopes_1 <-rep(NA,1000)
  for(i in 1:1000){
    sim_slopes_1[i] <- lm(arima.sim(list(ar=model_1$coef[1]), sd=sqrt(model_1$sigma2), n=span)~seq(1:span))$coefficients[2]
  }
  sd1 <- sd(sim_slopes_1)
  
  #store ensemble slopes 
  ensembleSlopes <-array(dim=c(35,40))
  
  #go through each end year
  for(endYear in 1:40){ #go through each endyear (i.e. ending in '61 to 2100)
    for(ensemble in 1:35){ #go through each ensemble
      ensembleSlopes[ensemble,endYear]<- lm(TGlobal_runs[(141-span+endYear):(140+endYear),ensemble]~seq(1:span))$coefficients[2] #2100 - 2060
    }
    varSlopeEndYears[endYear] <- var(ensembleSlopes[,endYear]) #store variance across ensembles for each given end year
  }
  
  #calculate ratio statistic
  ratios_31<-sd31/sqrt(mean(varSlopeEndYears))
  ratios_1<-sd1/sqrt(mean(varSlopeEndYears))
  
  
  #get confidence intervals
  boot_1 <- rep(NA,1000)
  boot_31<- rep(NA,1000)
  
  for(i in 1:1000){
    ensembles<-(sample.int(35,35,replace=TRUE))
    variance <- rep(NA,40)
    for(endYear in 1:40){
      variance[endYear] <- var(ensembleSlopes[ensembles,endYear]) #bootstrapped variances, 2100 - 2060
    }
    boot_1[i]<-sd1/sqrt(mean(variance))
    boot_31[i]<-sd31/sqrt(mean(variance))
  }
  
  #make confidence intervals from boot strap
  upConf_1 <-quantile(boot_1,c(.975))
  lowConf_1 <-quantile(boot_1,c(.025))
  upConf_31 <-quantile(boot_31,c(.975))
  lowConf_31 <-quantile(boot_31,c(.025))
  
  #add to dataframe
  sd_ratios <- rbind(sd_ratios,c(span, "ar1",ratios_1,lowConf_1,upConf_1)) 
  sd_ratios <- rbind(sd_ratios,c(span, "arma31",ratios_31,lowConf_31,upConf_31)) 
  
}

#plot
library(ggplot2)
sd_ratios <- sd_ratios[complete.cases(sd_ratios),]
sd_ratios$ratio<- as.numeric(sd_ratios$ratio)
sd_ratios$span<- as.numeric(sd_ratios$span)
limits <- aes(ymax = as.numeric(upperBound), ymin=as.numeric(lowerBound))

p <- ggplot(sd_ratios,aes(x=span,col=model,y=ratio))
p+geom_point()+geom_line()+geom_errorbar(limits)
