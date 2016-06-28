library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
source("R/format_data.R")

# Read both types of data
raw_data = read.csv("~/PostDoc/Sansom/Data_clones_David.csv", stringsAsFactors = F)
# raw_data$Day = factor(raw_data$Day)
raw_data$Condition = factor(raw_data$Condition)
raw_data           = filter(raw_data, CloneSize!=0)

# Format for Ploting
uu = format_exp_data(raw_data, time_points = NULL, clone_fractions = NULL)
pp1 = ggplot(uu, aes(col = Condition, y = value, x = Day)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin = low_lim, ymax = hi_lim)) + 
  ggtitle("David data") + theme(legend.position="bottomright") +
  labs(x = "Time post labelling (days)", y = "Clone fraction") +
  xlim(0, 100) + ylim(0,1) + 
  facet_grid(Condition ~ Nths)
plot(pp1)

Louis_raw         = data.frame(read.table("~/PostDoc/LouisFiles/LGR5_Louis/Louis_LGR5.WT_SI_8frac_mat.data", header = T))
time_points_Louis = data.frame(read.table("~/PostDoc/LouisFiles/LGR5_Louis/Louis_LGR5.WT_SI_8frac_mat.data", header = F))[1,]
Louis_list        = format_tbl2list(Louis_raw, time_points_Louis, condition = "WT_SI_Lgr5")
uu2               = format_exp_data(Louis_list, time_points = NULL, clone_fractions = NULL)

pp1 = ggplot(uu2, aes(col = Condition, y = value, x = Age)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin = low_lim, ymax = hi_lim)) + 
  ggtitle("David data") + theme(legend.position="bottomright") +
  labs(x = "Time post labelling (days)", y = "Clone fraction") +
  xlim(0, 100) + ylim(0,1) + 
  facet_grid(Condition ~ Nths)
plot(pp1)


## Make tabel from data Louis_list and raw_data (and back in the latter case?)
tblFormat_list      = format_list2tbl(raw_data, "Veh")
tblFormat_Louislist = format_list2tbl(Louis_list, "WT_SI_Lgr5")






