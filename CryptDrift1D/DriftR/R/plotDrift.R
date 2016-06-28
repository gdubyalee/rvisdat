# library(ggplot2)
# library(reshape2)
# library(plyr)

normalise.counts   <- function(x.mat, prior.param.dir = 0.05)
{
  apply(x.mat,2,function(x){(x+prior.param.dir)/(sum(x) + prior.param.dir*length(x))})
}

se.counts.multinom <- function(x.mat, prior.param.dir = 0.05)
{
  x.norm <- normalise.counts(x.mat)
  sd.aux <- x.norm*(1-x.norm)
  A <- prior.param.dir*nrow(x.mat)
  n <- colSums(x.mat)
  for(i in 1:ncol(x.mat)) sd.aux[,i] <- sd.aux[,i]/(A+n[i]+1)
  post.sd <- sqrt(sd.aux)
  2*post.sd
}

expandData <- function(x)
{
  x.exp <- NULL
  for(i in 1:length(x))
  {
    x.exp <- c(x.exp, rep(i, x[i]))
  }  
  x.exp
}

calcMean <- function(x.mat) apply(x.mat, 2, function(x) x%*%( 1:length(x) ) )

makeFractions <- function(x, fractions.split = NULL, time.interval)
{
  if(is.null(fractions.split)) fractions.split <- nrow(x)
  if((fractions.split == 4)&(nrow(x)!=4)) x  <- rbind(colSums(x[1:2,]), colSums(x[3:4,]), colSums(x[5:6,]),  colSums(x[7:8,]))
  
  #   fractions.split <- nrow(x)
  ## Plot biased drift fit
  clone.quarters  <- data.frame(Age = time.interval, t(normalise.counts(x)), Location = "Data")
  colnames(clone.quarters)[2:(1+fractions.split)] <- paste(1:fractions.split, ".", fractions.split, sep = "")
  colnames(clone.quarters)[1] <- "Age"
  quarters.melt <- melt(clone.quarters, id = c("Age", "Location"))
  # SE calculations
  se.clone.quarters <- data.frame(Age = time.interval, t(se.counts.multinom(x)), Location = "Data")  
  colnames(se.clone.quarters)[2:(1+fractions.split)] <- paste(1:fractions.split, ".", fractions.split, sep = "")
  colnames(se.clone.quarters)[1] <- "Age"
  se.melt   <- melt(se.clone.quarters, id = c("Age", "Location"))
  quarters.melt <- cbind(quarters.melt, se = se.melt[,4])
  quarters.melt
}


plotPulseChase <- function(sim, time_sim, x = NULL, time_x=NULL, tau=0)
{
  time_sim        = time_sim + tau
  fractions.split = nrow(sim)
  col_names       = c("Age", paste(1:fractions.split, ".", fractions.split, sep = ""))
  
  ## Sim process  
  model.preds           = data.frame(Age = time_sim, t(sim))
  colnames(model.preds) = col_names
  model.melt            = melt(model.preds, id = "Age")

  if(!is.null(x)){
    ## Data process
    clone.quarters           = data.frame(Age = time_x, t(normalise.counts(x)))
    colnames(clone.quarters) = col_names
    quarters.melt            = melt(clone.quarters, id = "Age")
    # SE calculations
    se.clone.quarters           = data.frame(Age = time_x, t(se.counts.multinom(x)))  
    colnames(se.clone.quarters) = col_names
    se.melt                     = melt(se.clone.quarters, id = "Age")
    quarters.melt               = cbind(quarters.melt, se = se.melt[,4])
    plot_data = geom_point(data = quarters.melt, aes(x=Age, y=value), size = 2.2) + 
    geom_errorbar(data = quarters.melt, aes(ymax = value + se, ymin=value - se), width=3) 
  }else{plot_data=NULL}
#   browser()  
    pp <- ggplot(data = model.melt, aes(x=Age, y=value)) + geom_line() + plot_data +
    labs(x = "Time post labelling (days)", y = "Clone fraction") +
    theme_bw() + facet_wrap(~variable, nrow = 1) + ggtitle("") # facet_grid(Location~variable) 
  pp
}


# ## Plot Anna data
# x.all <- read.table("All_Anna_Data_withWT.csv", header = T, sep = "\t")
# time.interval <- c(4,7,10,14,21)
# 
# x_hom_prox  <- as.matrix(subset(x.all, Het_hom == "Homozygous" & Location == "Proximal", select=1:5))
# x_hom_dist  <- as.matrix(subset(x.all, Het_hom == "Homozygous" & Location == "Distal", select=1:5))
# x_hom_colon  <- as.matrix(subset(x.all, Het_hom == "Homozygous" & Location == "Colon", select=1:5))
# x_het_prox  <- as.matrix(subset(x.all, Het_hom == "Heterozygous" & Location == "Proximal", select=1:5))
# x_het_dist  <- as.matrix(subset(x.all, Het_hom == "Heterozygous" & Location == "Distal", select=1:5))
# x_het_colon  <- as.matrix(subset(x.all, Het_hom == "Heterozygous" & Location == "Colon", select=1:5))
# 
# factor.levels <- c( "Arid1A Het Prox", "Arid1A Hom Prox", "Arid1A Het Dist", "Arid1A Hom Dist", "Arid1A Het Colon", "Arid1A Hom Colon")
# 
# 
# 
# 
# pdf("Anna_allData.pdf", 10, 8)
# clone.fracs <- 4
# max.time    <- 220
# num.tp      <- length(time.interval)
# num_elem    <- num.tp*clone.fracs
# all.data.nd <- rbind(makeFractions(x_hom_prox,  clone.fracs, time.interval),
#                      makeFractions(x_hom_dist,  clone.fracs, time.interval),
#                      makeFractions(x_hom_colon, clone.fracs, time.interval),
#                      makeFractions(x_het_prox,  clone.fracs, time.interval),
#                      makeFractions(x_het_dist,  clone.fracs, time.interval),
#                      makeFractions(x_het_colon, clone.fracs, time.interval))
# all.data.nd <- data.frame(all.data.nd, dataSet = factor(c(rep("Arid1A Hom Prox", num_elem), 
#                                                           rep("Arid1A Hom Dist", num_elem), 
#                                                           rep("Arid1A Hom Colon", num_elem),
#                                                           rep("Arid1A Het Prox", num_elem), 
#                                                           rep("Arid1A Het Dist", num_elem), 
#                                                           rep("Arid1A Het Colon", num_elem)),
#                                                         levels = factor.levels))
# 
# 
# pp <- ggplot(data = all.data.nd, aes(x=Age, y=value, colour = dataSet)) +
#   geom_point(size = 2.5) + geom_errorbar(aes(ymax = value + se, ymin=value - se), width=1) + #xlim(c(0,13)) +
#   labs(x = "Time post labelling (days)", y = "Clone fraction") +
#   theme_bw() + theme(legend.position = "none") + xlim(0,30) + 
#   facet_grid(dataSet~variable)
# print(pp)

#  plotContLabel


