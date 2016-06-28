fitBiasedDrift = function(x, time.interval, WT_params, max.iter = 80000, n.cpu = detectCores()-1, 
                          burn.in = 5000, thin = 50)
{
  Ns.WT     = WT_params$N
  lambda.WT = WT_params$lambda
  tau.WT    = WT_params$tau
  
  registerDoMC(cores=n.cpu)  
  all.runs = foreach(mcmc.run = 1:n.cpu, .combine= rbind)%dopar%{
    mcmc_out = MH_pulse_chase_Biased(x, time.interval, max.iter, Ns.WT, lambda.WT, tau.WT)
#     message("Finished!")
#     message(paste(round(mcmc_out$switching.MH, 4), collapse = " "))
    cbind(mcmc_out$p, mcmc_out$unsc.post)  
  }  
  ## All chains have been concatenated so to do burn in and thin need to do this 
  remove.me = NULL
  for(i in 1:n.cpu) remove.me <- c(remove.me, (i-1)*max.iter + (1:burn.in) )
  all.runs.subs = all.runs[-1*remove.me, ]
  all.runs.subs = all.runs.subs[seq(1,nrow(all.runs.subs), thin),]
  all.runs.subs
}

simulateBiasedDriftData = function(lambda, Ns, tau, p, time.samples, size_measure, num.crypts)
{
  Pn    <- th_PulseChase(lambda, Ns, tau, time.samples, persisting = T, Pr = p, splitNum = size_measure)
  xn    <- matrix(0, size_measure, length(time.samples))
  for(i in 1:length(time.samples))
  {
    x.i          <- sample(1:size_measure, num.crypts, replace = T, prob = Pn[,i])
    t.i          <- table(x.i)
    indx.i       <- as.numeric(names(t.i))
    xn[indx.i,i] <- t.i    
  }
  xn
}


plotsConvergence = function(all.runs.subs)
{
  sub.p         = all.runs.subs[,1]
  sub.unsc.post = all.runs.subs[,2]  
  pp1 = qplot(y = sub.p, x = 1:length(sub.p), geom = "line") + xlab("Iteration (thinned)") + ylab("Pr")
  pp2 = qplot(y = sub.unsc.post, x = 1:length(sub.unsc.post), geom = "line") + xlab("Iteration (thinned)") + 
    ylab("Unscaled Posterior")
  ac.p = acf(sub.p, plot=F)
  acd  = data.frame(lag=ac.p$lag, acf=ac.p$acf)
  pp3  = ggplot(acd, aes(x=lag, y=acf)) + geom_bar(stat="identity") +
    geom_hline(yintercept=c(0.05, -0.05), linetype="dashed") 

  len.1.chain = round(length(sub.p)/2)
  aux_df = as.data.frame(qqplot(sub.p[1:len.1.chain], sub.p[(len.1.chain+1):length(sub.p)], 
                                plot.it=FALSE))
  pp4 = ggplot(aux_df) + geom_point(aes(x=x, y=y)) + geom_abline(slope=1,colour="red") + xlab("First half of run") +
    ylab("Second half of run")
  
  grid.arrange(pp1, pp2,
               pp3, pp4, nrow=2)
  
}

plotPosterior_Pr = function(all.runs.subs)
{
  all.mcmc   = data.frame(prob = all.runs.subs[,1])  
  vals.cols2 = c("#999999", "darkgreen", "steelblue" ,"#66CC00")
  pp =  ggplot(all.mcmc, aes(x = prob)) +
    geom_histogram(binwidth = 0.002, col = vals.cols2[3], fill = vals.cols2[3]) + 
    labs(x = "Mutant cell fitness", y = "Probability density") +
    geom_vline(xintercept = 0.5, lty =2, size = 1.5) + #ylim(0, 600)+
    ggtitle("Distribution for bias parameter") + xlim(0,1)
  #     ggtitle(bquote("Arid1a"^"-/-"))
  pp
}


plotsBiasDrift_Fit = function(x, time.interval, WT_params, Pr_est, max_x = 100, condition)
{
  # Transform format for plotting
  data_list = format_tbl2list(x, time.interval, condition = condition)
  uu2       = format_exp_data(data_list, clone_fractions = nrow(x) )
  col_p = "steelblue"
  pp1 = ggplot(uu2, aes(y = value, x = Age)) + geom_point(size = 3, col = col_p) + 
    geom_errorbar(aes(ymin = low_lim, ymax = hi_lim), col = col_p) + 
    ggtitle("Biased drift fit") + #theme(legend.position="bottomright") +
    labs(x = "Time post labelling (days)", y = "Clone fraction") +
    xlim(0, max_x) + ylim(0,1) + theme_bw(24) +
    facet_grid(Condition ~ Nths)
  # plot(pp1)
  
  time_sim        = WT_params$tau:max_x
  sim_WT          = th_PulseChase(lambda = WT_params$lambda, Ns = WT_params$N, tau = WT_params$tau,
                                  time_points = time_sim, splitNum = nrow(x))
  sim_mut   = th_PulseChase(lambda = WT_params$lambda, Ns = WT_params$N, tau = WT_params$tau,
                            Pr = Pr_est[2], time_points = time_sim, splitNum = nrow(x))
  sim_WT    = format_th_data(th_x = sim_WT, time_points = time_sim, condition)
  sim_mut   = format_th_data(th_x = sim_mut, time_points = time_sim, condition)
  cols = c("WT"="orange","Mutant"=col_p, "Dummy = #111111")
  pp <- pp1 + scale_colour_manual(name="Type", values=cols)  +
    geom_line(data = sim_WT, aes(col = "WT"), lty = 2)   +
    geom_line(data = sim_mut, aes(col = "Mutant"))
  pp  
}

getBiasParams = function(all.runs.subs)
{
  sub.p        = all.runs.subs[,1]
  quantile(sub.p, probs = c(0.025, 0.5, 0.975))
}

