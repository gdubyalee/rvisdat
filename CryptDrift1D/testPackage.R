## Cont labeling sim/theory
alpha       = 1e-4
lambda      = 0.1
Ns          = 5 
time_points = seq(1, 500, length.out = 100)
uu1 = th_contLabelling(alpha, lambda, Ns, time_points )
uu2 = sim_contLabelling(alpha, lambda, Ns, time_points, 100000000, Pr = 0.5, raw_sims = F)

plot(time_points, uu1[1,])
lines(time_points, uu2[1,])


time_points = seq(1, 100, length.out = 100)
source("../MasterEq.functions.R")
uu1 = NeutralDrift(lambda, Ns, time_points)
uu2 = th_PulseChase(lambda, Ns, time_points, persisting = F, splitNum = 0)

