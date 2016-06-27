library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(dplyr)
library(stringr)
library(DriftR)
library(InferCryptDrift)

theme_set(theme_bw(18))

raw_data2 = data.frame(t(read.csv("LGK974_day4_7_10_21_complete_raw_data 2.csv", stringsAsFactors = F, header = F)))
WT_betcat = data.frame(read.csv("proc_WT_pooledReps.csv", stringsAsFactors = F, header = T)) 
WT_betcat = format_tbl2list(WT_betcat, time_points = c(4,10,14,21), condition = "WT") %>% select(Condition, Day, CloneSize)
betcat = data.frame(read.csv("proc_Bcat_pooledReps.csv", stringsAsFactors = F, header = T)) 
betcat = format_tbl2list(betcat, time_points = c(4,10,14,21), condition = "BetaMut") %>% select(Condition, Day, CloneSize)

raw_data = raw_data2 %>% separate(X1, c("Condition", "Day"), sep ="_day") %>% gather(dummy, CloneSize, -Condition, -Day) %>%
  select(-dummy) %>% filter(CloneSize!="", CloneSize!=0) %>% 
  mutate(Condition = str_replace(Condition, "LGK", "LGK974")) %>%
  mutate(Condition = str_replace(Condition, "WT", "VEH")) %>% 
  mutate(Condition = str_replace(Condition, "VEH", "Veh")) %>% 
  mutate(CloneSize = as.numeric(CloneSize))

all_data = rbind(raw_data, WT_betcat) #, betcat)



## Neutral drift clone size dist
summary_means = all_data %>% group_by(Day, Condition) %>% 
  summarise(meanCloneSize = mean(CloneSize), se = sqrt(var(CloneSize)/length(CloneSize)))

# Append zero
summary_means = rbind(summary_means, 
                      data.frame(Day = 0, Condition = "LGK974", meanCloneSize = 0, se = 0),
                      data.frame(Day = 0, Condition = "Veh", meanCloneSize = 0, se = 0))


pp_means = ggplot(summary_means, aes(x = as.numeric(Day), y = meanCloneSize, col = Condition)) + geom_point(size = 4) + 
  geom_line(size = 1.2) + 
  geom_errorbar(aes(ymin = meanCloneSize - 1.96*se, ymax = meanCloneSize + 1.96*se), size = 1.2, width = 1) +
  ylab("Average clone size") + xlab("Day") + ylim(0, 8) + xlim(0, 22)+ 
  ggtitle("Change of average clone size with time")
# plot(pp_means)



fitAndPlot_NeutralDrift = function(data_x, time.interval, condition, plot4ths = F, max_x = 100)
{
  if(plot4ths)
  {
    data_x_4ths = rbind(colSums(data_x[1:2,]),
                       colSums(data_x[3:4,]),
                       colSums(data_x[5:6,]),
                       colSums(data_x[7:8,]))
  }  
  # condition = "WT"
  fit_out   = fitNeutralDrift(data_x, time.interval)
  plotsConvergence_Neutral(fit_out)
  plotPosterior_Neutral(fit_out)
  WT_est  = getNeutralDirftParams(fit_out)
  pp_neut = plotsNeutralDrift_Fit(data_x, time.interval, WT_est, max_x = max_x, condition)
  plot(pp_neut+ theme(plot.margin = unit(c(6,1,6,1), "cm")))
  if(plot4ths){
    pp_neut4 = plotsNeutralDrift_Fit(data_x_4ths, time.interval, WT_est, max_x = max_x, condition)
    plot(pp_neut4 + theme(plot.margin = unit(c(6,1,6,1), "cm")))
  }
}

## Fit neutral drift
x_VEH = format_list2tbl(all_data, "Veh")[,-1]
x_WT  = format_list2tbl(all_data, "WT")[,-1]
x_LGK = format_list2tbl(all_data, "LGK974")[,-1]
## Collapse to 4ths
x_LGK_4ths = rbind(colSums(x_LGK[1:2,]),
                   colSums(x_LGK[3:4,]),
                   colSums(x_LGK[5:6,]),
                   colSums(x_LGK[7:8,]))
## Collapse to 4ths
x_WT_4ths =  rbind(colSums(x_WT[1:2,]),
                    colSums(x_WT[3:4,]),
                    colSums(x_WT[5:6,]),
                    colSums(x_WT[7:8,]))

## Collapse to 4ths
x_VEH_4ths =  rbind(colSums(x_VEH[1:2,]),
                    colSums(x_VEH[3:4,]),
                    colSums(x_VEH[5:6,]),
                    colSums(x_VEH[7:8,]))

time.interval_LGK = as.numeric(colnames(x_VEH))
time.interval_WT  = as.numeric(colnames(x_WT))


### Plots
pdf("David_All_NeutDrift.pdf", 14, 9)
plot(pp_means)
# WT
grid.arrange(textGrob("Fits for WT"))
fitAndPlot_NeutralDrift(     x_WT, time.interval_WT, condition = "WT", plot4ths = T, max_x = 100)
fitAndPlot_NeutralDrift(x_WT_4ths, time.interval_WT, condition = "WT", plot4ths = F, max_x = 100)
grid.arrange(textGrob("Fits for VEH"))
fitAndPlot_NeutralDrift(     x_VEH, time.interval_LGK, condition = "VEH", plot4ths = T, max_x = 100)
fitAndPlot_NeutralDrift(x_VEH_4ths, time.interval_LGK, condition = "VEH", plot4ths = F, max_x = 100)
grid.arrange(textGrob("Fits for LGK"))
fitAndPlot_NeutralDrift(     x_LGK, time.interval_LGK, condition = "LGK", plot4ths = T, max_x = 100)
fitAndPlot_NeutralDrift(x_LGK_4ths, time.interval_LGK, condition = "LGK", plot4ths = F, max_x = 100)
dev.off()







