unscaled.log.posterior_bias <- function(x, time.intervals, lambda, Ns, p, tau)
{
  if(lambda<0 | Ns<2 | tau > time.intervals[1]| p > 1 | p < 0)
  {
    unsc.log.post <- -999999999999
  }else{
    size_measure  <- nrow(x)  
    log_all_preds <- log(th_PulseChase(lambda, Ns, tau, time.intervals, persisting = T, Pr = p, splitNum = size_measure))
    log_all_preds <- sum(log_all_preds*x)
    if(!is.finite(log_all_preds))
    {
      all_preds     <- th_PulseChase(lambda, Ns, tau, time.intervals, persisting = T, Pr = p, splitNum = size_measure)
      all_preds[all_preds==0] <- 1 ## Correct for zero values, by discarding that data value
      log_all_preds <- sum(log(all_preds)*x)	
    }
    shape.val     <- 0.0001 # 1 
    rate.val      <- 0.0001 # 2
    # To correct case where delay tau is switched off tau = 0
    if(tau == 0){aux_tau <- 0}else{aux_tau <- dgamma(tau, shape = shape.val, rate = rate.val, log = T)}
    log.prior    <- dbeta(p, shape1 = 0.5, shape2 = 0.5, log = T) + aux_tau
    # N has a constant prior
    unsc.log.post <- log_all_preds + log.prior
  }
  unsc.log.post
}

propse.p <- function(p, sd = 0.1)
{
  sd <- sample(c(0.01, 0.25), 1, prob = c(0.8, 0.2))
  rnorm(n = 1, mean = p, sd = sd)
}


MH_pulse_chase_Biased <- function(x, time.intervals, max.iter, Ns.WT, lambda.WT, tau.fix)
{
  ## Make sure x is a matrix! 
  x <- as.matrix(x)
  ten.percent       <- max.iter/10
  p.mcmc            <- rep(0, max.iter)
  unsc.post.mcmc    <- rep(0, max.iter)
  p.mcmc[1]         <- rbeta(1,0.5,0.5)
  log.unsc.post     <- unscaled.log.posterior_bias(x, time.intervals, lambda.WT, Ns.WT, p.mcmc[1], tau.fix)
  
  p.switch          <- 0
  
  unsc.post.mcmc[1] <- log.unsc.post
  for(i in 1:(max.iter-1))
  {
    ## Update p =================================
    p.new             <- propse.p(p.mcmc[i])
    log.unsc.post.new <- unscaled.log.posterior_bias(x, time.intervals, lambda.WT, Ns.WT, p.new, tau.fix)
    h.ratio           <- min(0,log.unsc.post.new - log.unsc.post)
    rnum.i            <- log(runif(1))
    if(rnum.i < h.ratio)
    {
      p.mcmc[i+1]   <- p.new
      log.unsc.post <- log.unsc.post.new
      p.switch      <- p.switch + 1 
    }
    else{
      p.mcmc[i+1]  <- p.mcmc[i]
    }
    
    unsc.post.mcmc[i+1] <- log.unsc.post
    if(i%%ten.percent==0) message(paste("Done ", 100*i/max.iter, "%"))
  }  
  list(p = p.mcmc, unsc.post = unsc.post.mcmc, 
       switching.MH = c(p = p.switch/max.iter))
}
MH_pulse_chase_Biased <- cmpfun(MH_pulse_chase_Biased)


