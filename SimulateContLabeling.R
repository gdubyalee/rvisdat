library(Rcpp)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DriftR)
theme_set(theme_bw(14))

sim_contLabelling_revert <- function(mu, lambda, Ns, time_points, numSim, Pr = 0.5, raw_sims = F, revert_mu = 0)
{
  beta  <- 2*lambda*Pr
  alpha <- 2*lambda*(1-Pr)
  #   Number of chain reactions
  M      <- Ns
  # Initial state vector
  x0             = c(numSim, rep(0,M)) 
  names(x0)      = paste("x",seq(M+1),sep="") 
  nu             = matrix(rep(0,(M*(M+1))),ncol=M)
  diag(nu)       = -1
  diag(nu[2:M,]) = +1
  nu[M+1,M]      = +1
  nu             = cbind(nu, -1*nu)
  # nu             = nu[,-1*c(ncol(nu))]

  a_mat          = cbind(c(lambda*Ns*mu, rep(beta,M-1), rep(alpha,M-1), lambda*Ns*revert_mu), 
                         c(seq(M), 2:(M+1)))
  
  sims <- LinearGillespie(numSim=1, nu, a_mat, x0, time_points) 
  if(!raw_sims){
    return_val = rbind(colSums(sims[2:Ns,])/numSim, sims[Ns+1,]/numSim)
  }else{return_val = sims}
  return_val
}

runSims = function(mu = 1.1*10^(-4), revert_mu = 0, simName="",
                   Pr = 0.5, lambda = 0.3, Ns = 7, cryptCounted = 1e7,time_points  = seq(1, 200, length.out=200))
{
  sims_out = sim_contLabelling_revert(mu, lambda, Ns, time_points, numSim = cryptCounted, Pr, raw_sims = F, revert_mu = revert_mu)
  if(revert_mu>0) revertant = "Revertant"
  else revertant = "NoRevertants"
  data.frame(mono = sims_out[2,], partial = sims_out[1,], time = time_points, type = simName, revertant = revertant, stringsAsFactors = F)
}

convert2counts = function(sims_out, time_points, type = "")
{

  data.frame(counts_monoclonal = qbinom(p = 0.5, sims_out),
             ymin_mono         = mon
             )
}

# ## Colon parameters
# mu_CA30  = 1.1*10^(-4)
# all_sim_freq = rbind(runSims(mu =      mu_CA30, revert_mu = 0, simName="CA30"),
#       runSims(mu =   10*mu_CA30, revert_mu = 0, simName="CA30x10"),
#       runSims(mu =  100*mu_CA30, revert_mu = 0, simName="CA30x100"),
#       runSims(mu =  50*mu_CA30, revert_mu = 0, simName="CA30x50"),
#       runSims(mu = 1000*mu_CA30, revert_mu = 0, simName="CA30x1000"),
#       runSims(mu =      mu_CA30, revert_mu = mu_CA30, simName="CA30"),
#       runSims(mu =   10*mu_CA30, revert_mu = 10*mu_CA30, simName="CA30x10"),
#       runSims(mu =  100*mu_CA30, revert_mu = 100*mu_CA30, simName="CA30x100"),
#       runSims(mu =  50*mu_CA30, revert_mu = 50*mu_CA30, simName="CA30x50"),
#       runSims(mu = 1000*mu_CA30, revert_mu = 1000*mu_CA30, simName="CA30x1000"))   
# 
# pp1 = ggplot(all_sim_freq, aes(x = time, y = 100*mono)) + geom_line(col = "blue", size = 1.5) + 
#   geom_line(aes(y = 100*partial, x = time), col = "blue", alpha = 0.6, size = 1.5) + 
#                    ylab("Fraction of crypts (%)") + xlab("Days") + #ggtitle("Neutral Drift") + 
#   facet_grid(type~revertant,scales = "free")
# pdf("ContinousLabel_alpha.pdf")
# plot(pp1)
# dev.off()


## Colon parameters
mu_CA30   = 1.1*10^(-4)
mu        = 300*mu_CA30
revert_mu = 0.5*mu
sims_out = sim_contLabelling_revert(mu, 0.3, 7, seq(1, 200, length.out=200), numSim = 1e5, 0.5, raw_sims = T, revert_mu = revert_mu)
sims_out = t(sims_out[-1,]/sum(sims_out[,1]))
colnames(sims_out) = paste("N", 1:7, sep ="")
print(100*sims_out[102,7])

sim_list     = data.frame(mono = sims_out[,7], partial = rowSums(sims_out[,-7]), time = time_points, 
                      type = "Simulation", stringsAsFactors = F)

sim_list_sep = data.frame(sims_out, time_points=time_points) %>% gather(key = MutSC, value = Fraction, - time_points)
  

pp1 = ggplot(sim_list, aes(x = time, y = 100*mono)) + geom_line(col = "blue", size = 1.5) +
  geom_line(aes(y = 100*partial, x = time), col = "blue", alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(28, 70, 102, 121), lty = 2, alpha = 0.4, col = "red") + 
  ylab("Fraction of crypts (%)") + xlab("Days") #ggtitle("Neutral Drift") 
# plot(pp1)

pp2 = ggplot(sim_list_sep, aes(x = time_points, y = 100*Fraction, col = MutSC)) + geom_line(size = 1.5) +
    geom_vline(xintercept = c(28, 70, 102, 121), lty = 2, alpha = 0.4, col = "red") + 
  ylab("Fraction of crypts (%)") + xlab("Days") #ggtitle("Neutral Drift") 
# plot(pp1)

# pdf("ContinousLabel_alpha.pdf")
grid.arrange(pp1, pp2)
# dev.off()
# 
