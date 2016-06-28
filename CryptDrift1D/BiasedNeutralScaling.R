library(reshape2)
library(ggplot2)
# source('~/Dropbox/BiasedDrift/BiasedDrift.R')
# source("MasterEq.functions.R")
Ns     <- 32
lambda <- 4
Pr.1   <- 0.4
Pr.2   <- 0.6
# alpha  <- 3.5
# beta   <- 0.5
time_intervals <- seq(3, 9, length.out=3)
neutral        <- NeutralDrift(lambda, Ns, time_intervals)
pers.neutral   <- t(apply(neutral[2:(Ns+1),], 1, function(x, y = neutral[1,]) x/(1-y)))
pers.biased1   <- biased_drift(alpha = 2*(1-Pr.1)*lambda, beta= 2*Pr.1*lambda, Ns, 
                                time_intervals, persisting = T, splitNum = 0)
pers.biased2   <- biased_drift(alpha = 2*(1-Pr.2)*lambda, beta= 2*Pr.2*lambda, Ns, 
                                time_intervals, persisting = T, splitNum = 0)
day_names      <- paste("day", round(time_intervals, 0))
colnames(pers.neutral) <- colnames(pers.biased1) <- colnames(pers.biased2) <- day_names

neutral_melt       <- cbind(melt(pers.neutral), Drift_type = "Neutral drift")
biased_melt1       <- cbind(melt(pers.biased1), Drift_type = paste("Biased Pr =", Pr.1))
biased_melt2       <- cbind(melt(pers.biased2), Drift_type = paste("Biased Pr =", Pr.2))
all_data           <- rbind(neutral_melt, biased_melt1, biased_melt2)
colnames(all_data) <- c("CloneSize", "TimePoint", "Frequency", "DriftType")
all_data$TimePoint <- factor(as.character(all_data$TimePoint), levels = day_names)
data_proc          <- ddply(all_data, .(TimePoint, DriftType), mutate, avge_n = rep(Frequency%*%1:Ns, Ns), 
                            avge_n_overn = CloneSize/avge_n)
x <- seq(0, 5, length.out=50)
scaling_function <- 0.5*pi*x*exp(-0.25*pi*x^2)
scaling_function <- data.frame(x = x, y = scaling_function)
pdf("Scaling_biasedDrift.pdf")
pp <- ggplot(data_proc, aes(y = avge_n*Frequency, x = avge_n_overn)) + geom_point(size = 4, alpha = 0.5) +
  geom_line(data=scaling_function, aes(x = x, y = y), col = "red", lty = 2) +
  facet_grid(DriftType~TimePoint) + xlab("<n>/n") + ylab("<n> x Freq ")
plot(pp)
dev.off()








