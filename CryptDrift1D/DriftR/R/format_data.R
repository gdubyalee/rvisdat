#' Transform tbl format to list format. Useful for ggplot plotting
format_tbl2list = function(tbl_x, time_points, condition = "Condition1")
{
  list_data = data.frame(CloneSize = NULL, Day = NULL, Condition = NULL)
  for(i in 1:nrow(tbl_x))
  {
    for(j in 1:ncol(tbl_x)) 
    {
      if(tbl_x[i,j]!=0) list_data = rbind(list_data, 
                                          data.frame(CloneSize = rep(i, tbl_x[i,j]),
                                                     Day       = as.character(time_points[j]),
                                                     Condition = condition,                                                                    
                                                     stringsAsFactors = F))
    }  
  }
  list_data
}

#' Transform list format to table format with clone size along rows and time points
#' along cols. Will return time points as colnames. This is the format used for inference
#' If a clone size does not appear in any data point it will not be included (!!!)
format_list2tbl = function(list_x, condition_filter = NULL, clone_fractions = NULL)
{
  # Check what type of data
  if(is.null(clone_fractions))
  {
    clone_fractions = max(as.numeric(list_x$CloneSize))
    message(paste("Assuming clones have been counted in 1/", clone_fractions,"s",sep=""))    
  }
  
#   # Filter and remove 
#   list_filt = filter(list_x, Condition == condition_filter) %>% select(-Condition)
#   # Make sure order is correct and all fractions present
#   list_filt$CloneSize = factor(list_filt$CloneSize, levels = 1:clone_fractions) 
#   list_filt$Day       = as.character(list_filt$Day)
#   ordered_levels      = sort(as.numeric(unique(list_filt$Day)))
#   list_filt$Day       = factor(list_filt$Day, levels = ordered_levels)
#   
#   ret_tbl = list_filt                            %>% 
#             group_by(Day, CloneSize)             %>%
#             summarise(Counts = n())              %>%
#             spread(Day, Counts, drop = F)
# Filter and remove Condition col
  if(!is.null(condition_filter)) list_x = filter(list_x, Condition == condition_filter)

  list_filt = list_x %>% select(CloneSize, Day)%>%
  mutate(CloneSize = factor(CloneSize, levels = 1:clone_fractions)) %>% 
  mutate(num_Age = factor(as.numeric(str_extract(Day, pattern = "[0-9]+"))))# Get numerical age

  ret_tbl = list_filt                            %>% 
            group_by(num_Age, CloneSize)         %>%
            summarise(Counts = n())              %>%
            spread(num_Age, Counts, drop = F)

  ret_tbl[is.na(ret_tbl)] = 0 # Remove NAs
  ret_tbl
}

#' Add clones that haven't been counted 
#' If clones haven't been counted they don't appear rather than 
#' count as zero
addMissingVals = function(data_Nths, frac_size)
{
  clone_sizes_all = paste(1:frac_size, frac_size, sep = ".")
  # Data frame with missing vals in third col
  missing_vals = data_Nths %>% group_by(Condition, Age) %>% 
    summarise(missing = paste(clone_sizes_all[!(clone_sizes_all%in%Nths)], collapse = ","))
  # Loop and add missing vals
  for(i in 1:nrow(missing_vals))
  {
    if(missing_vals$missing[i]!= "") 
    {
      missing_i = str_split(missing_vals$missing[i], pattern = ",")
      for(jj in missing_i) data_Nths = rbind(data_Nths, data.frame(Day       =        missing_vals$Age[i],
                                                                   Nths      =                         jj, 
                                                                   Age       =        missing_vals$Age[i], 
                                                                   Condition =  missing_vals$Condition[i], 
                                                                   value     =                          0, 
                                                                   low_lim   =                          0, 
                                                                   hi_lim    =                          0))
    }
  }
  data_Nths
}


#' Calculate proportions with 95\% credible interval (confidence interval).
#' If there are replicates they will be pooled(!). It will work for both count table format and list format.
#' For it to work for table format give time vector
format_exp_data = function(data_x, time_points = NULL, clone_fractions = NULL, condition = "Condition1")
{
  # If time points given assume it's in table format
  if(!is.null(time_points)) data_x = format_tbl2list(data_x, time_points, condition)
  
  # Check right headers are there, ie CloneSize, Condition, Day
  necessary_names = c("CloneSize", "Condition", "Day")
  if(!( all(necessary_names %in% colnames(data_x)))) stop(paste("Can't find", paste(necessary_names, collapse=","), "in colnames"))
  
  # Check what type of data
  if(is.null(clone_fractions))
  {
    clone_fractions = max(as.numeric(data_x$CloneSize))
    message(paste("Assuming clones have been counted in 1/", clone_fractions,"s",sep=""))    
  }
  
  Frac_summarise = data_x %>% mutate(Nths = paste(CloneSize, clone_fractions, sep = "."))
  
  data_Nths = Frac_summarise %>% 
    group_by(Condition, Nths, Day) %>% summarise(Counts = n()) %>%
    group_by(Condition, Day)       %>%    mutate(TotalCounts = sum(Counts), 
                                                 value       = qbeta(  0.5, Counts, TotalCounts-Counts),
                                                 low_lim     = qbeta(0.025, Counts, TotalCounts-Counts), 
                                                 hi_lim      = qbeta(0.975, Counts, TotalCounts-Counts)) %>%          
    mutate(variable = Nths, Age = as.numeric(str_replace(Day, "_days",""))) %>%
    select(Nths, Age, Condition, value, low_lim, hi_lim)
  
  # Add clones that haven't been counted, ie counts of zero don't appear, so add manually if missing
  data_Nths = addMissingVals(data_Nths, clone_fractions)  
  data_Nths  
}

#' Adjust theory drift data for plotting
format_th_data = function(th_x, time_points, condition = "Condition1")
{
  fractions.split = nrow(th_x)
  col_names       = c("Age", paste(1:fractions.split, ".", fractions.split, sep = ""), "Condition")
  
  ## Sim process  
  model.preds           = data.frame(Age = time_points, t(th_x), Condition = condition)
  colnames(model.preds) = col_names
  model.list            = model.preds %>% gather(Nths, value, c(-Age, -Condition))
  model.list
}


