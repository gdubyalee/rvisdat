library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DriftR)
library(InferCryptDrift)
theme_set(theme_bw(24))

rawData<-read.csv('rawdata.csv')
time.interval = c(4, 7, 10, 14, 21)

# Fraction of crypt 
condition      = "LGR5_YFP"

fit_out = fitNeutralDrift(rawData, time.interval)

plotsConvergence_Neutral(fit_out)
pdf('neutralconvergence1.pdf',14,9)

plotPosterior_Neutral(fit_out)

WT_est = getNeutralDirftParams(fit_out)
pp     = plotsNeutralDrift_Fit(rawData, time.interval, WT_est, max_x = 100, condition)
plot(pp + theme(plot.margin = unit(c(4,1,4,1), "cm")))

dev.off()


