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
#data_raw       = read.table("~/PostDoc/ANNA/Arid1A_24-Feb/clone size data overview_Prox.csv", sep = ",", stringsAsFactors = F, na.strings = "", skip = 2)

#This is for an unbiased case
rawDataUnbiased<-read.csv('rawdata.csv')
time.interval  = c(4, 7, 10, 14, 21)
condition      ="LGR5_YFP"

#WT_params = list(lambda = 0.0949, N = 5, tau = 0.8134)\
fit_out<-fitNeutralDrift(rawData, time.interval)

#Wild type paramaters
WT_est = getNeutralDirftParams(fit_out)

fit_out   = fitBiasedDrift(x, time.interval, WT_est)

plotsConvergence(fit_out)
plotPosterior_Pr(fit_out)
Pr_est = getBiasParams(fit_out)
pp = plotsBiasDrift_Fit(x, time.interval, WT_params, Pr_est, max_x = 100, condition)
plot(pp)



