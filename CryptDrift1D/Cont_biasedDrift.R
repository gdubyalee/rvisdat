library(Rcpp)
sourceCpp("gillespie_Linear.cpp")
# source("MasterEq.functions.R")
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

Pr     <- 0.32
lambda <- 0.1
Ns     <- 5
mu     <- 1.1*10^(-4) 
cryptCounted <- 1e7
time_points <- seq(1, 1000, length.out=100) #0.5*(1:200)

sim_data <- sim_contLabelling(mu, lambda, Pr, Ns, time_points, numSim = cryptCounted)
plot(y= sim_data[6,], x = time_points, type = "l")
lines(y= colSums(sim_data[2:5,]), x = time_points)
abline(a = 0, b = mu*lambda*cryptCounted, col = "red", lty = 2)
abline(a = cryptCounted*0.5*mu*Ns*(Ns-1), b = 0, col = "red", lty = 2)

# predicted N and lambda 
N_infrd      <- calc.N_noTAC(sim_data[2:5,20:ncol(sim_data)], mu, cryptCounted)
lambda_infrd <- calc.lambda(sim_data[6,20:ncol(sim_data)], time_points[20:ncol(sim_data)], mu, cryptCounted) 
print(paste(round(lambda_infrd,4), round(N_infrd,4)))
sim_data_infrd <- sim_contLabelling(mu, lambda_infrd, 0.5, round(N_infrd,0), time_points, numSim = cryptCounted)

pdf("Biased_nonBiased.pdf")
plot(y= sim_data[6,], x = time_points, type = "l", xlab = "Mouse age", ylab = "Clones", 
     main = "Same slope and partials\nBlack lines: lambda = 0.1, N=5, Pr = 0.32
     Red lines: lambda = 0.013, N=4, Pr = 0.5",cex.main=0.8)
lines(y= colSums(sim_data[2:5,]), x = time_points)
lines(y= sim_data_infrd[N_infrd+1,], x = time_points, type = "l", col = "red")
lines(y= colSums(sim_data_infrd[2:N_infrd,]), x = time_points, col = "red")
dev.off()

# plot(sim_data[6,]/colSums(sim_data[2:5,]), type = "l")
# lines(sim_data_infrd[N_infrd+1,]/colSums(sim_data_infrd[2:N_infrd,]), col ="red")
uu          <- NeutralDrift(0.3, 7, 1:100)
  