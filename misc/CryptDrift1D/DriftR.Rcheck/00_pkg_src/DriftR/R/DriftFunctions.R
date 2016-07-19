sim_contLabelling <- function(mu, lambda, Ns, time_points, numSim, Pr = 0.5, raw_sims = F)
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
  
  sims <- LinearGillespie(numSim=1, nu, a_mat, x0, time_points) 
  if(!raw_sims){
    return_val = rbind(colSums(sims[2:Ns,])/numSim, sims[Ns+1,]/numSim)
  }else{return_val = sims}
  return_val
}

th_PulseChase <- function(lambda, Ns, tau, time_points, persisting = T, Pr = 0.5, splitNum = 0, alpha_beta = NULL)
{
  if (time_points[1] < tau) stop("The first time point must be larger than tau!!!!")
  if(!is.null(alpha_beta))
  {
    alpha = alpha_beta[1]
    beta  = alpha_beta[2]
  }else{
    beta   = 2*lambda*Pr
    alpha  = 2*lambda*(1-Pr)    
  }
  Crypt_drift_c(alpha, beta, Ns, time_points - tau, persisting, splitNum)  
}


CA30Function_MonoClonal <- function(alpha, lambda, Ns, time.intervals)
{
  theory.timeseries <- rep(0, length(time.intervals))
  m       <- 1:(Ns-1)
  intercept <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2))
  aplha_lamb <- alpha*lambda
  for(i in 1:length(time.intervals))
  {
    t <- time.intervals[i]
    theory.timeseries[i] <- 0.5*alpha*sum((-1)^(m+1)*(tan(0.5*pi*m/Ns))^(-2)*(1-exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2)))
  }  
  aplha_lamb*time.intervals - theory.timeseries
}

CA30Function_PartialClones <- function(alpha, lambda, Ns, time.intervals)
{
  theory.timeseries   <- rep(0, length(time.intervals))
  theory.timeseries.n <- rep(0, length(time.intervals))
  m       <- 1:(Ns-1)
  n.all   <- 1:(Ns-1)
  for(n in n.all)
  {
    for(i in 1:length(time.intervals))
    {
      t                      <- time.intervals[i]
      theory.timeseries.n[i] <- 0.5*alpha*sum(sin(pi*m/Ns)*sin(pi*m*n/Ns)/(sin(0.5*pi*m/Ns))^2*exp(-4*lambda*t*(sin(0.5*pi*m/Ns))^2))
    }  
    theory.timeseries <- theory.timeseries+theory.timeseries.n   
  }
  # Return value
  alpha*Ns*(Ns-1)/2 - theory.timeseries
}


th_contLabelling <- function(alpha, lambda, Ns, time.intervals)
{
  rbind(CA30Function_PartialClones(alpha, lambda, Ns, time.intervals), 
        CA30Function_MonoClonal(alpha, lambda, Ns, time.intervals))  
}
  
