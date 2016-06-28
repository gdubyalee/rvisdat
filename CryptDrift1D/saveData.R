library(stringr)
file_name         = "data/Louis_LGR5.APCmin_8frac_mat.data"

data_x            = read.table(file_name, header = T)
colnames(data_x)  = str_replace(colnames(data_x), "X", "") #c(4, 7, 14, 28, 56, 126, 210)

Louis_LGR5.APCmin = data_x
save(Louis_LGR5.APCmin, file = "data/Louis_LGR5.APCmin.Rdata")

