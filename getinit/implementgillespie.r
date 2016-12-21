
library(foreach)
library(doMC)
registerDoMC(8)
library(Rcpp)
sourceCpp('gillespie_replaceAndPush.cpp')

sim1DCrypt<-function(startConf,lambda,numSim,times){
  Ns = length(startConf)
  # Use 0.5 as 1/2 the time it goes right or left (as two different updates)
  lambda_cc <- lambda/4 # There are 4 possible ways to do this move (due to border displacement)
  lambda_bb <- 0 
  lambda_cb <- 0
  lambda_bc <- 0
  
  #   Number of reactions
  numb_states = 2*Ns # 2 rings
  # Initial state vector
  x0          = rep(0,2*Ns)
  # if(length(Init)!=2) stop("Init must be a vector of type Init = c(centre=0,border=1)")
  # if(Init[1]>0) x0[1:Init[1]]      <- 1
  # if(Init[2]>0) x0[Ns+(1:Init[2])] <- 1
  x0[1:Ns]    = startConf
  ##  Reaction matrix ======================================================
  # Reactions cc           -----------------------------------------
  # cc right 
  nu_cc_r            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_cc_r)      <- 1
  diag(nu_cc_r[-1,]) <- -1
  nu_cc_r[1,Ns]      <- -1
  # cc left
  nu_cc_l            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_cc_l)      <- 1
  diag(nu_cc_l[,-1]) <- -1
  nu_cc_l[Ns,1]      <- -1
  # Replacements at border -----
  # Starting with right centre move
  nu_border_aux_l            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_border_aux_l)      <- -2
  nu_border_aux_r            <- matrix(0, ncol = Ns, nrow = Ns)  
  diag(nu_border_aux_r[-1,]) <- -2
  nu_border_aux_r[1,Ns]      <- -2
  nu_cc_r_aux1<- rbind(cbind(nu_cc_r, nu_cc_r), cbind(nu_border_aux_l, nu_border_aux_r))
  
  # Starting with left centre move
  nu_border_aux_l            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_border_aux_l[-(1:(Ns-2)),])      <- -2
  diag(nu_border_aux_l[,-(1:2)]) <- -2
  nu_border_aux_r            <- matrix(0, ncol = Ns, nrow = Ns)  
  diag(nu_border_aux_r[,-1]) <- -2
  nu_border_aux_r[Ns,1]      <- -2
  nu_cc_l_aux1               <- rbind(cbind(nu_cc_l, nu_cc_l), cbind(nu_border_aux_l, nu_border_aux_r))
  nu_cc <- cbind(nu_cc_l_aux1, nu_cc_r_aux1)
  
  # Reactions bb           -----------------------------------------
  ring_movement <- cbind(nu_cc_l,nu_cc_r)
  aux_mat <- matrix(0, ncol = 2*Ns, nrow = Ns)
  nu_bb   <- rbind(aux_mat, ring_movement)
  
  # Reactions bc           -----------------------------------------
  # Right move
  nu_bc_r_aux1            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_bc_r_aux1)      <- 1
  diag(nu_bc_r_aux1[-1,]) <- -2
  nu_bc_r_aux1[1,Ns]      <- -2
  
  nu_bc_r_aux2            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_bc_r_aux2[-1,]) <- -1
  nu_bc_r_aux2[1,Ns]      <- -1
  nu_bc_r                <- rbind(nu_bc_r_aux2, nu_bc_r_aux1)
  
  # Left move
  nu_bc_l_aux1            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_bc_l_aux1)      <- 1
  diag(nu_bc_l_aux1[,-1]) <- -2
  nu_bc_l_aux1[Ns,1]      <- -2
  
  nu_bc_l_aux2            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_bc_l_aux2)      <- -1
  nu_bc_l                 <- rbind(nu_bc_l_aux2, nu_bc_l_aux1)
  
  nu_bc                   <- cbind(nu_bc_l, nu_bc_r)
  
  # Reactions cb           -----------------------------------------
  # Right move
  nu_cb_aux1            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_cb_aux1)      <- 1
  
  nu_cb_r_aux2            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_cb_r_aux2)      <- -1
  nu_cb_r                <- rbind(nu_cb_aux1, nu_cb_r_aux2)
  
  # Left move
  nu_cb_l_aux2            <- matrix(0, ncol = Ns, nrow = Ns)
  diag(nu_cb_l_aux2[,-1]) <- -1
  nu_cb_l_aux2[Ns,1]      <- -1
  nu_cb_l                 <- rbind(nu_cb_aux1, nu_cb_l_aux2)
  
  nu_cb                   <- cbind(nu_cb_l, nu_cb_r)
  
  # Put all nu matrices together
  nu      = cbind(nu_cc, nu_bb, nu_bc, nu_cb)
  
  # Reaction rate constants
  a_mat   = cbind( c(rep(lambda_cc, 4*Ns), rep(lambda_bb, 2*Ns), rep(lambda_bc, 2*Ns), rep(lambda_cb, 2*Ns)),    0)  
  
  # Run simulations
  sim_out = replaceAndPushGillespie(numSim, nu, a_mat, x0, time_points)
  sim_out[1:(Ns+1),]
}


#startConf  = c(1, 0, 0, 1, 0, 0, 0)
nStartingConfig<-1000
numSim<-10000
nStem<-7
lambda<-.2
finTime<-30
time_points<-1:finTime
qq=c()
fracs<-seq(0,.5,.05)
for(fracExpressing in fracs){
  qq <-append(
    qq,
    list(foreach(i = 1:nStartingConfig, .combine = "+")%dopar%{
      print(paste0('Run ',i))
      #We should sample _with_ replacement...
      #Better still just give a proportion of cells expressing
      startConf  <- runif(nStem)<fracExpressing
      sim1DCrypt(startConf, lambda, numSim, time_points)
    }/(numSim*nStartingConfig))
  )
}

fin<-lapply(1:length(fracs),function(i){qq[i][[1]][,finTime]})
propPar<-unlist(lapply(fin,function(entry){sum(entry[2:nStem])}))
propClone<-unlist(lapply(fin,function(entry){sum(entry[nStem+1])}))
print(head(propPar))
print(head(propClone))

library(ggplot2)
dat<-data.frame(propPar=propPar,propClone=propClone,fracs=fracs)

plt<-ggplot(dat,aes(y=propPar,x=fracs))+geom_line()
print(plt)

