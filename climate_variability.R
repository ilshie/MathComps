library(dplyr)
library(tidyr)
library(ggplot2)

########################## Parameters ############################################
# number of repetitions for parametric bootstrap
n <- 10000
# number of repetitions for non-parametric bootstrap
m <- 10000
# minimum and maximum trend lengths to consider
min_trend_length <- 10
max_trend_length <- 50
##################################################################################

# Load coupled control temperature data
if (!exists("TGlobal")) load("~/Classes/MATH400/CoupledControl.RData")
coupled_TGlobal <- TGlobal
remove(lat, lon, TRegional, TLocal, TGlobal)

# Subset last 500 years of global temperature data
temperatures <- tail(coupled_TGlobal, 500)

# Fit models to coupled control global temperatures
(ar1 <- arima(temperatures, order = c(1,0,0)))
(arma31 <- arima(temperatures, order = c(3,0,1)))

# Fit white noise model to coupled control global temperatures
(white_noise <- arima(temperatures, order = c(0,0,0)))

# Function that carries out the parametric bootstrap for a given list of ARMA 
# coefficients, standard error, and trend length; returns the standard error 
# across all computed slopes
parametric_bootstrap <- function(ARMA_coeff, ARMA_SE, trend_length, num_rep) {
  slopes <- sapply(1:num_rep, function(i) {
      # simulate trend_length time steps from fitted model
      sim <- arima.sim(model=ARMA_coeff, sd=ARMA_SE, n=trend_length)
      # fit linear trend to simulated time series
      sim.lm <- lm(sim ~ c(1:trend_length))
      # return estimated slope
      return(sim.lm$coefficients[2])
  })
  # return standard error across all computed slopes
  return(sd(slopes))
}

# obtain vector of standard errors across different trend lengths for both models
ar1_SE <- sapply(min_trend_length:max_trend_length, function(i) {
  parametric_bootstrap(list(ar=0.5045), sqrt(ar1$sigma2), i, n)
}) 
arma31_SE <- sapply(min_trend_length:max_trend_length, function(i) {
  parametric_bootstrap(list(ar=c(1.2958, -0.7934, 0.3583), ma=-0.6278), 
                       sqrt(arma31$sigma2), i, n)
})

# obtain vector of standard errors across different trend lengths for white noise
white_noise_SE <- sapply(min_trend_length:max_trend_length, function(i) {
  parametric_bootstrap(list(ar=0), sqrt(white_noise$sigma2), i, n)
})

############################################################################################

# Function that fits a linear trend of given length trend_length to overlapping 
# intervals of the last 500 years of coupled control temperature data; returns
# all 450 computed slopes
overlapping_trends <- function(trend_length) {
    slopes <- sapply(1:450, function(i) {
        # fit linear trend of length trend_length
        trend.lm <- lm(temperatures[i:(i + trend_length - 1)] ~ c(1:trend_length))
        # return estimated slope
        return(trend.lm$coefficients[2])
    })
    # return standard error across all computed slopes
    return(as.vector(slopes)) 
}

# obtain matrix of 450 slopes per trend length
slope_matrix <- data.frame(sapply(min_trend_length:max_trend_length, function(i) {
  overlapping_trends(i)
})) 

# for each trend length, bootstrap m times on the 450 slopes to compute the mean
# and bootstrap standard error for the slope variances
df <- slope_matrix %>%
  gather(X1:X41, key=trend_length, value=slope) %>%
  mutate(trend_length = extract_numeric(trend_length)+9) %>%
  group_by(trend_length) %>%
  do({
    SE_vector <- sapply(1:m, function(i) {
                    bootstrap_sample <- sample_n(., size=450, replace=TRUE)
                    return(sd(bootstrap_sample$slope))
                 })
    data.frame(trend_length=.[1,1], true_mean=mean(SE_vector), 
               true_lower=quantile(SE_vector, 0.025), true_upper=quantile(SE_vector, 0.975))
  }) %>%
  ungroup()

############################################################################################

# add SE derived from ARMA models to data frame
df$AR1 <- ar1_SE
df$ARMA31 <- arma31_SE
df$white <- white_noise_SE

# plot computed standard errors against trend length
df %>% gather(true_mean:white, key=type, value=SE) %>%
  ggplot(aes(x=trend_length, y=SE)) +
    geom_line(aes(color=type, linetype=type)) +
    scale_color_manual(values=c(2,3,1,1,1,4)) +
    scale_linetype_manual(values=c(1,1,2,1,2,1)) +
    labs(x="\nTrend length (year)", y="Estimated standard error\n") +
    scale_y_continuous(breaks=seq(0, 0.015, 0.003)) +
    theme(axis.title=element_text(size=16), axis.text=element_text(size=14),
          legend.title=element_text(size=16), legend.text=element_text(size=14))

# plot ratio of ARMA standard error to true standard error against trend length
df2 <- df %>% transmute(trend_length, 
    AR1_mean=AR1/true_mean, AR1_upper=AR1/true_upper, AR1_lower=AR1/true_lower, 
    ARMA31_mean=ARMA31/true_mean, ARMA31_upper=ARMA31/true_upper, ARMA31_lower=ARMA31/true_lower,
    white_mean=white/true_mean, white_upper=white/true_upper, white_lower=white/true_lower)

p <- df2 %>% gather(AR1_mean:white_lower, key=model, value=SE) %>%
  ggplot(aes(x=trend_length, y=SE)) +
    geom_line(aes(color=model, size=model)) +
    geom_line(y=1, linetype=2) +
    scale_color_manual(breaks=c("AR1_mean", "ARMA31_mean", "white_mean"), 
                       values=rep(c(2,3,4), each=3), labels=c("AR1", "ARMA31", "white noise")) +
    scale_size_manual(values=rep(c(0.2,1.2,0.2),3), guide=FALSE) +
    labs(x="\nTrend length (year)", y="Ratio of estimated to true standard error\n") +
    theme(axis.title=element_text(size=16), axis.text=element_text(size=14),
          legend.title=element_text(size=16), legend.text=element_text(size=14))

df3 <- df2 %>% select(-c(AR1_mean, ARMA31_mean))
  
p + geom_ribbon(data=df3, inherit.aes=FALSE,
                aes(x=trend_length, ymin=AR1_lower, ymax=AR1_upper), 
                fill=2, alpha=0.1) +
  geom_ribbon(data=df3, inherit.aes=FALSE,
              aes(x=trend_length, ymin=ARMA31_lower, ymax=ARMA31_upper),
              fill=3, alpha=0.1) +
  geom_ribbon(data=df3, inherit.aes=FALSE,
              aes(x=trend_length, ymin=white_lower, ymax=white_upper),
              fill=4, alpha=0.1)