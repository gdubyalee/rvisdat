unscaled.log.posterior <- function(x, time.intervals, lambda, Ns, tau)
{
  if(lambda<0 | Ns<2 | tau > time.intervals[1])
  {
    unsc.log.post <- -999999999999
  }else{
    size_measure  <- nrow(x)
    
    log_all_preds <- log(th_PulseChase(lambda, Ns, tau, time.intervals, persisting = T, Pr = 0.5, splitNum = size_measure))
    log_all_preds <- sum(log_all_preds*x)
    if(!is.finite(log_all_preds))
    {
      all_preds     <- th_PulseChase(lambda, Ns, tau, time.intervals, persisting = T, Pr = 0.5, splitNum = size_measure)
      all_preds[all_preds==0] <- 1 ## Correct for zero values, by discarding that data value
      log_all_preds <- sum(log(all_preds)*x)	
    }
    shape.val     <- 0.0001 # 1 
    rate.val      <- 0.0001 # 2
    # To correct case where delay tau is switched off tau = 0
    #       if(tau == 0){aux_tau <- 0}else{aux_tau <- dgamma(tau, shape = shape.val, rate = rate.val, log = T)}
    log.prior    <- dgamma(lambda, shape = shape.val, rate = rate.val, log = T) #+ aux_tau
    # N has a constant prior
    unsc.log.post <- log_all_preds + log.prior
  }
  unsc.log.post
}

propse.N <- function(Ns, range.jump = 1)
{
  jump.interval <- 1:range.jump
  N.choose      <- c(Ns-jump.interval, Ns+jump.interval)
  jump.porb     <- dnorm(N.choose, mean = Ns, sd = range.jump/3)
  jump.porb     <- jump.porb/sum(jump.porb)
  sample(N.choose, size = 1, prob = jump.porb)
}

# propse.N.2 <- function(Ns, prob.jump = 0.2)
# {
#   N.choose      <- c(Ns-1, 0, Ns+1)
#   sample(N.choose, size = 1, prob = c(prob.jump, 1-2*prob.jump, prob.jump))
# }


propse.lambda <- function(lambda, sd = 0.1)
{
  sd <- sample(c(0.01, 0.25), 1, prob = c(0.8, 0.2))
  rnorm(n = 1, mean = lambda, sd = sd)
}

propse.tau <- function(tau, sd = 0.1)
{
  sd <- sample(c(1, 3), 1, prob = c(0.8, 0.2))
  rnorm(n = 1, mean = tau, sd = sd)
}


jointLambdaN_MH <- function(curr.N, curr.lambda, sd.MH = 0.001)
{
  new.N              <- propse.N(curr.N)
  curr.s             <- curr.lambda/curr.N^2
  new.lambda         <- rnorm(n = 1,  mean = curr.s*new.N^2, sd = sd.MH*new.N^2)
  new.s              <- new.lambda/new.N^2
  q.old.g.new        <- dnorm(curr.lambda, mean = new.s*curr.N^2, sd = sd.MH*curr.N^2)
  q.new.g.old        <- dnorm(new.lambda,  mean = curr.s*new.N^2, sd = sd.MH*new.N^2)
  log.ratio.proposal <- log(q.old.g.new/q.new.g.old)
  #   new.lambda         <- curr.lambda/curr.N^2*new.N^2
  #   log.ratio.proposal <- 0
  c(new.N = new.N, new.lambda = new.lambda, log.ratio.proposal = log.ratio.proposal)
}


MH_pulse_chase <- function(x, time.intervals, max.iter, N.interval, tau.interval, lambda.interval)
{
  tau.fix <- 0
  ## Make sure x is a matrix! 
  x <- as.matrix(x)
  ten.percent       <- max.iter/10
  Ns.mcmc           <- rep(0, max.iter)
  lambda.mcmc       <- rep(0, max.iter)
  tau.mcmc          <- rep(0, max.iter)
  unsc.post.mcmc    <- rep(0, max.iter)
  Ns.mcmc[1]        <- sample(N.interval[1]:N.interval[2], 1) #5, 40
  lambda.mcmc[1]    <- 1/16^2*Ns.mcmc[1]^2 ## Start in reasonable area #runif(1, lambda.interval[1], lambda.interval[2])
  tau.mcmc[1]       <- tau.fix #runif(1, tau.interval[1], tau.interval[2])
  log.unsc.post     <- unscaled.log.posterior(x, time.intervals, lambda.mcmc[1], Ns.mcmc[1], tau.mcmc[1])
  N.switch          <- 0
  lambda.switch     <- 0
  tau.switch        <- 0
  
  unsc.post.mcmc[1] <- log.unsc.post
  for(i in 1:(max.iter-1))
  {
    ## Update N =====================================
    new.vals <- jointLambdaN_MH(Ns.mcmc[i], lambda.mcmc[i], sd.MH = 8*10^(-5))
    new.N    <- new.vals[1]; lambda.new <- new.vals[2]; log.ratio.proposal <- new.vals[3]
    # Check if in interval
    if(new.N > N.interval[1] & new.N < N.interval[2] & lambda.new > lambda.interval[1] & lambda.new < lambda.interval[2])
    {
      log.unsc.post.new <- unscaled.log.posterior(x, time.intervals, lambda.new, new.N, tau.mcmc[i]) 
      h.ratio           <- min(0, log.unsc.post.new - log.unsc.post + log.ratio.proposal)
      rnum.i            <- log(runif(1))
      if(rnum.i < h.ratio)
      {
        Ns.mcmc[i+1]      <- new.N
        lambda.mcmc[i+1]  <- lambda.new
        log.unsc.post     <- log.unsc.post.new
        N.switch          <- N.switch + 1 
      }else{
        Ns.mcmc[i+1]      <- Ns.mcmc[i]
        lambda.mcmc[i+1]  <- lambda.mcmc[i]
      }
    }else{
      Ns.mcmc[i+1]      <- Ns.mcmc[i]      
      lambda.mcmc[i+1]  <- lambda.mcmc[i]
    }
    ## Update lambda =================================
    lambda.new        <- propse.lambda(lambda.mcmc[i+1], 0.1)
    if(lambda.new > lambda.interval[1] & lambda.new < lambda.interval[2])
    {
      log.unsc.post.new <- unscaled.log.posterior(x, time.intervals, lambda.new, Ns.mcmc[i+1], tau.mcmc[i])
      h.ratio           <- min(0,log.unsc.post.new - log.unsc.post)
      rnum.i            <- log(runif(1))
      if(rnum.i < h.ratio)
      {
        lambda.mcmc[i+1]  <- lambda.new
        log.unsc.post     <- log.unsc.post.new
        lambda.switch     <- lambda.switch + 1 
      }
    }
    ## Update tau =================================
    tau.new <- propse.tau(tau.mcmc[i])
    if(tau.new > tau.interval[1] & tau.new < tau.interval[2])
    {
      log.unsc.post.new <- unscaled.log.posterior(x, time.intervals, lambda.mcmc[i+1], Ns.mcmc[i+1], tau.new)
      h.ratio           <- min(0,log.unsc.post.new - log.unsc.post)
      rnum.i            <- log(runif(1))
      if(rnum.i < h.ratio)
      {
        tau.mcmc[i+1]  <- tau.new
        log.unsc.post  <- log.unsc.post.new
        tau.switch     <- tau.switch + 1 
      }else{
        tau.mcmc[i+1]  <- tau.mcmc[i]
      }
    }else{
      tau.mcmc[i+1]  <- tau.mcmc[i]
    }     
    
    unsc.post.mcmc[i+1] <- log.unsc.post
    if(i%%ten.percent==0) message(paste("Done ", 100*i/max.iter, "%"))
  }  
  list(lambda = lambda.mcmc, N = Ns.mcmc, tau = tau.mcmc, unsc.post = unsc.post.mcmc, 
       switching.MH = c(tau = tau.switch/max.iter, lambda = lambda.switch/max.iter, Ns = N.switch/max.iter))
}
MH_pulse_chase <- cmpfun(MH_pulse_chase)

