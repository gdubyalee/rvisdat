library(Rcpp)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(plyr)

sourceCpp("gillespie_Linear.cpp")
calc.N_noTAC <- function(partials, mu, crypts_counted) {
  part_mean <- mean(colSums(partials))/crypts_counted
  (1+sqrt( 1+8*part_mean/mu))/2
}
calc.lambda  <- function(y, Age, mu, crypts_counted) coef(lm(y~Age))[2]/(mu*crypts_counted)

sim_contLabelling <- function(mu, lambda, Pr, Ns, time_points, numSim)
{
  beta  <- 2*lambda*Pr
  alpha <- 2*lambda*(1-Pr)
  #   Number of chain reactions
  M      <- Ns
  # Initial state vector
  x0             <- c(numSim, rep(0,M)) 
  names(x0)      <- paste("x",seq(M+1),sep="") 
  nu             <- matrix(rep(0,(M*(M+1))),ncol=M)
  diag(nu)       <- -1
  diag(nu[2:M,]) <- +1
  nu[M+1,M]      <- +1
  nu             <- cbind(nu, -1*nu)
  nu             <- nu[,-1*c(ncol(nu))]
  a_mat          <- cbind(c(lambda*Ns*mu, rep(beta,M-1), rep(alpha,M-1)), c(seq(M),2:M))
  
  LinearGillespie(numSim=1, nu, a_mat, x0, time_points) 
}

Pr     <- 0.7
lambda <- 0.5
Ns     <- 10
mu     <- 1.1*10^(-4) 
cryptCounted <- 1e7
time_points <- seq(1, 200, length.out=200) #0.5*(1:200)

pdf("ContinousLabelling_clonesize.pdf", 9 ,7)
theme_set(theme_bw(14))

check_timepoints <- c(10, 30, 200)
neutral_cont <- sim_contLabelling(mu, lambda, 0.5, Ns, time_points, numSim = cryptCounted)[-1,]
pp1 <- ggplot(data.frame(time = time_points, partial = colSums(neutral_cont[1:(Ns-1),]), full = neutral_cont[Ns,]), 
       aes(x = time, y = full)) + geom_line() + geom_line(aes(y = partial, x = time_points)) + 
       geom_vline(xintercept = time_points[check_timepoints], col = "red", lty = 2) + ylab("Clones")+ 
  ggtitle("Neutral Drift")


neutral_cont <- apply(neutral_cont, 2, function(x) x[-length(x)]/sum(x[-length(x)]))
colnames(neutral_cont) <- paste("t=",round(time_points,2))
size_neutral <- data.frame(neutral_cont[, check_timepoints], frac_crypt = 1:(Ns-1)/Ns)
size_neutral <- melt(size_neutral, id = "frac_crypt")
pp2 <- ggplot(size_neutral, aes(x = frac_crypt, y = value, col = variable)) + geom_point() + geom_line()+
  xlab("Fraction of the crypt") + ylab("Percentage of the clone")+ 
  ggtitle("Neutral Drift")


# predicted N and lambda 
check_timepoints <- c(10, 30, 200)
biased_cont <- sim_contLabelling(mu, lambda, Pr, Ns, time_points, numSim = cryptCounted)[-1,]
pp3 <- ggplot(data.frame(time = time_points, partial = colSums(biased_cont[1:(Ns-1),]), full = biased_cont[Ns,]), 
              aes(x = time, y = full)) + geom_line() + geom_line(aes(y = partial, x = time_points)) + 
  geom_vline(xintercept = time_points[check_timepoints], col = "red", lty = 2) + ylab("Clones") + 
  ggtitle(paste("Bias Pr=", Pr))

biased_cont <- apply(biased_cont, 2, function(x) x[-length(x)]/sum(x[-length(x)]))
colnames(biased_cont) <- paste("t=",round(time_points,2))
size_neutral <- data.frame(biased_cont[, check_timepoints], frac_crypt = 1:(Ns-1)/Ns)
size_neutral <- melt(size_neutral, id = "frac_crypt")
pp4 <- ggplot(size_neutral, aes(x = frac_crypt, y = value, col = variable)) + geom_point() + geom_line()+
  xlab("Fraction of the crypt") + ylab("Percentage of the clone")+
  ggtitle(paste("Bias Pr=", Pr))
grid.arrange(pp1,pp2,pp3,pp4, ncol = 2)
  

x_sizes <- seq(0, 1, length.out=500)
seq_pos_fracs <- as.matrix( seq(1/Ns, (Ns-1)/Ns, by = 1/Ns))
qq <- apply(seq_pos_fracs, 1, function(x, x_sizes) dnorm(x=x_sizes, mean= x, sd = 0.01), x_sizes)
colnames(qq) <- 1:(Ns-1)
qq2 <- melt(data.frame(qq, x_pos = x_sizes), id.vars="x_pos")

pp1 <- ggplot(qq2, aes(x = x_pos, y = value, col=variable)) + geom_line()
size_noise <- ddply(size_neutral, .(variable), summarise, size_dens = qq%*%value, x_sizes=x_sizes)
pp2 <- ggplot(size_noise, aes(x = x_sizes, y = size_dens, col=variable)) + geom_line()

qq <- apply(seq_pos_fracs, 1, function(x, x_sizes) dnorm(x=x_sizes, mean= x, sd = 0.05), x_sizes)
colnames(qq) <- 1:(Ns-1)
qq2 <- melt(data.frame(qq, x_pos = x_sizes), id.vars="x_pos")
pp3 <- ggplot(qq2, aes(x = x_pos, y = value, col=variable)) + geom_line()

size_noise <- ddply(size_neutral, .(variable), summarise, size_dens = qq%*%value, x_sizes=x_sizes)
pp4 <- ggplot(size_noise, aes(x = x_sizes, y = size_dens, col=variable)) + geom_line()
grid.arrange(pp1,pp2,pp3,pp4, ncol = 2) 
  
dev.off()




# pdf("Biased_nonBiased.pdf")
# plot(y= sim_data[6,], x = time_points, type = "l", xlab = "Mouse age", ylab = "Clones", 
#      main = "Same slope and partials\nBlack lines: lambda = 0.1, N=5, Pr = 0.32
#      Red lines: lambda = 0.013, N=4, Pr = 0.5",cex.main=0.8)
# lines(y= colSums(sim_data[2:5,]), x = time_points)
# lines(y= sim_data_infrd[N_infrd+1,], x = time_points, type = "l", col = "red")
# lines(y= colSums(sim_data_infrd[2:N_infrd,]), x = time_points, col = "red")
# dev.off()
# 
# # plot(sim_data[6,]/colSums(sim_data[2:5,]), type = "l")
# # lines(sim_data_infrd[N_infrd+1,]/colSums(sim_data_infrd[2:N_infrd,]), col ="red")
# uu          <- NeutralDrift(0.3, 7, 1:100)
  