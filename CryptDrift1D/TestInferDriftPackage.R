library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DriftR)
theme_set(theme_bw(24))

transf.data <- function(louis.data)
{
  fraction.split <- max(louis.data, na.rm = T)
  all.measure <- matrix(0, fraction.split, ncol(louis.data))
  for(i in 1:ncol(louis.data))
  {
    dd <- table(louis.data[,i])
    for(j in 1:length(dd)) all.measure[as.numeric(names(dd)[j]), i] <- dd[j]
  }
  all.measure
}
data_raw       = read.table("~/PostDoc/ANNA/Arid1A_24-Feb/clone size data overview_Prox.csv", sep = ",", stringsAsFactors = F, na.strings = "", skip = 2)

# Fraction of crypt 
fractions_data = data_raw[,c(1,3,5,7,9)]/data_raw[,c(1,3,5,7,9)+1]
fractions_data = ceiling(fractions_data*8)
x              = transf.data(fractions_data)
time.interval  = c(4, 7, 10, 14, 21)
condition      = "Arid1A"

WT_params = list(lambda = 0.0949, N = 5, tau = 0.8134)
fit_out   = fitBiasedDrift(x, time.interval, WT_params)
  
plotsConvergence(fit_out)
plotPosterior_Pr(fit_out)
Pr_est = getBiasParams(fit_out)
pp = plotsBiasDrift_Fit(x, time.interval, WT_params, Pr_est, max_x = 100, condition)
plot(pp)



