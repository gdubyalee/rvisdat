source("R/format_data.R")
time_sim        = 1:100
sim_1           = th_PulseChase(lambda = 0.1, Ns = 5,tau = 1,time_points = time_sim, splitNum = 8)
sim_2           = th_PulseChase(lambda = 1, Ns = 16,tau = 1,time_points = time_sim, splitNum = 8)
model.list      = rbind(format_th_data(th_x = sim_1, time_points = time_sim, condition = "Mine"),
                        format_th_data(th_x = sim_2, time_points = time_sim, condition = "Bens"))

# fractions.split = nrow(sim)
# col_names       = c("Age", paste(1:fractions.split, ".", fractions.split, sep = ""), "Condition")
# 
# ## Sim process  
# model.preds           = data.frame(Age = time_sim, t(sim), Condition = "aa")
# colnames(model.preds) = col_names
# model.list            = model.preds %>% gather(Nths, value, -Age, -Condition)

pp <- ggplot(data = model.list, aes(x=Age, y=value, col = Condition)) + geom_line() + 
  labs(x = "Time post labelling (days)", y = "Clone fraction") +
  theme_bw() + facet_wrap(~Nths, nrow = 1) + ggtitle("") # facet_grid(Location~variable) 
plot(pp)

