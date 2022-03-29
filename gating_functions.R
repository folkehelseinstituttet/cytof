#' find_gate_lower_noise_plus_high_low_selected_cells find position for gating for each subdataset
#' @param data fcs_data sets
#' @param channel which channel used for gating
#' @return value of gating for the given channel in each subdataset. 

find_gate_lower_noise_plus_high_low_selected_cells <- function(data, channel, include, mark, positiv = TRUE){
  column <- which(colnames(data[[1]]) == channel)
  results <- rep(NA, length(data))
  for(i in 1:length(data)){
    if(positiv){
      res <-  find_gate_lower_noise_per_file_plus_high_low(xx = data[[i]][include[[i]][,mark], column])
    } else {
      res <-  find_gate_lower_noise_per_file_plus_high_low(xx = data[[i]][!include[[i]][,mark], column])
    }
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates))  
}




#' find_gate_lower_noise_plus_high_low find position for gating for each subdataset
#' @param data fcs_data sets
#' @param channel which channel used for gating
#' @return value of gating for the given channel in each subdataset. 

find_gate_lower_noise_plus_high_low <- function(data, channel){
  column <- which(colnames(data[[1]]) == channel)
  results <- rep(NA, length(data))
  for(i in 1:length(data)){
    results[i] <- find_gate_lower_noise_per_file_plus_high_low(data[[i]][, column])
  }
  return(results)
}


#' find_gate_lower_noise_per_file_plus_high_low function
#' @param xx vector of values
#' @return the value that correspond to the percentage in the density plot.

find_gate_lower_noise_per_file_plus_high_low <- function(xx){
  dens <- density(xx)
  #  ts_y<-ts(smooth(dens$y))
  #  tp <- pastecs::turnpoints(ts_y)
  #  top1 <- dens$x[tp$peaks][1]
  top1 <- dens$x[dens$y == max(dens$y)]
  firstGate <- top1 + (top1 - dens$x[1])
  
  xxx <- xx[xx > firstGate]
  dens <- density(xxx)
  browser()
  
  return(value)
}




#' find_gate_lower_noise find position for gating for each subdataset
#' @param data fcs_data sets
#' @param channel which channel used for gating
#' @return value of gating for the given channel in each subdataset. 

find_gate_lower_noise <- function(data, channel){
  column <- which(colnames(data[[1]]) == channel)
  results <- rep(NA, length(data))
  for(i in 1:length(data)){
    results[i] <- find_gate_lower_noise_per_file_based_on_x_top(data[[i]][, column])
  }
  return(results)
}


#' find_gate_lower_noise_per_file function
#' @param xx vector of values
#' @return the value that correspond to the percentage in the density plot.

find_gate_lower_noise_per_file_based_on_x_top <- function(xx){
  dens <- density(xx)
#  ts_y<-ts(smooth(dens$y))
#  tp <- pastecs::turnpoints(ts_y)
#  top1 <- dens$x[tp$peaks][1]
  top1 <- dens$x[dens$y == max(dens$y)]
 value <- top1 + (top1 - dens$x[1])
  return(value)
}




#' find_gate_perc_upper_noise find position for gating for each subdataset
#' @param data fcs_data sets
#' @param channel which channel used for gating
#' @param upper_perc how many percentages that are assumed to be noise
#' @return value of gating for the given channel in each subdataset. 

find_gate_perc_height_upper_noise <- function(data, channel, upper_perc_height = 0.001){
  column <- which(colnames(data[[1]]) == channel)
  results <- rep(NA, length(data))
  for(i in 1:length(data)){
    results[i] <- find_gate_upper_noise_per_file(data[[i]][, column], upper_perc_height = upper_perc_height)
  }
  return(results)
}


#' find_gate_upper_noise_per_file function
#' @param xx vector of values
#' @param upper_perc how many percentage to trow away.
#' @return the value that correspond to the percentage in the density plot.

find_gate_upper_noise_per_file <- function(xx, upper_perc_height){
  dens <- density(xx)
  posible <- dens$y > max(dens$y) * upper_perc_height
  xposible <- which(posible[2:length(posible)] - posible[1:(length(posible)-1)] == -1)[1]
  value <- dens$x[xposible]
  return(value)
}


#' find_gate_gaussian_first_top_per_file function
#' @param xx vector of values
#' @param upper_perc how many percentage to trow away.
#' @return the value that correspond to the percentage in the density plot.
find_gate_gaussian_first_top_per_file <- function(xx, perc_included = 0.9995){
  fit <- density(xx)
  x.new <- rnorm(10000, sample(xx, size = 10000, replace = TRUE), fit$bw)
  #plot(fit)
  #lines(density(x.new), col = "red")
  fit <- VGAM::vglm(x.new ~ 1, 
                    VGAM::mix2normal(eq.sd = F), iphi=0.5, imu= 0, isd1=2, imu2=5, 
                    isd2=1)
  #  fit2 <- vglm(x.new ~ 1, uninormal(), lmean = 0, lsd = 1)
  
  # Calculated parameters
  pars <- as.vector(coef(fit))
  w <- VGAM::logitlink(pars[1], inverse=TRUE, )
  m1 <- pars[2]
  sd1 <- exp(pars[3])
  m2 <- pars[4]
  sd2 <- exp(pars[5])
  gate <- m1 + qnorm(perc_included) * sd1
  return(gate)
}


#' find_gate_gaussian_first_top find position for gating for each subdataset
#' @param data fcs_data sets
#' @param channel which channel used for gating
#' @param perc_included how many percentages that are assumed to be noise
#' @return value of gating for the given channel in each subdataset. 

find_gate_gaussian_first_top <- function(data, channel, perc_included = 0.9995){
  column <- which(colnames(data[[1]]) == channel)
  results <- rep(NA, length(data))
  for(i in 1:length(data)){
    results[i] <- find_gate_gaussian_first_top_per_file(data[[i]][, column], perc_included = perc_included)
  }
  return(results)
}




#' events_to_keep, find which events to keep in each subdataset based on lower and/or upper gate.
#' @param data, transformed data 
#' @param channel, which channel to plot
#' @param lower_gate, vector with values for lower gate, a number (same lower gate for all subset) or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate, a number (same upper gate for all subset)  or NA (no upper gating)
#' @return list of vectors of true/false for each events (true means keep), one vector for each subdataset. 


events_to_keep <- function(data, channel, lower_gate = NA, upper_gate = NA){
  number_of_files <- length(data)
  column <- which(colnames(data[[1]]) == channel)
  if(length(lower_gate) == 1){
    lower_gate <- rep(lower_gate, number_of_files)
  }
  if(length(upper_gate) == 1){
    upper_gate <- rep(upper_gate, number_of_files)
  }
  kept_events <- NULL
  if(is.na(lower_gate[1])){
    for (i in 1:number_of_files){
      kept_events[[i]] <- data[[i]][,column] < upper_gate[i]
    }
  } else { 
    if(is.na(upper_gate[1])){
      for (i in 1:number_of_files){
        kept_events[[i]] <- data[[i]][,column] > lower_gate[i]
      }
    } else {     
      for (i in 1:number_of_files){
        kept_events[[i]] <- data[[i]][,column] < upper_gate[i] & data[[i]][,column] > lower_gate[i]
      }
    }
  }  
  return(kept_events)
}


#' percent_to_keep_this_gating
#' @param data, data 
#' @param kept_events, list of vectors of true/false for alle events in each subdataset. Those events that are true will be kept
#' @param file_names, default NA will give the name 1,2,3 etc.
#' @return percentage kept this gating
#' 
percent_to_keep_this_gating <- function(kept_events, file_names = NA){
  number_of_files <- length(kept_events)
  if(is.na(file_names[1])){
    file_names <- 1:number_of_files
  }
  file_names <- as.character(file_names)
  percent <- rep(NA, number_of_files)
  for (i in 1:number_of_files){
    percent[i] <- table(kept_events[[i]])["TRUE"]/length(kept_events[[i]])
  }
  names(percent) <- file_names
  return(percent)
}


#' update_data_based_on_events_to_keep
#' @param data, data 
#' @param kept_events, list of vectors of true/false for alle events in each subdataset. Those events that are true will be kept
#' @return new dataset with only those events that we want to keep. 

update_data_based_on_events_to_keep <- function(data, kept_events){
  number_of_files <- length(data)
  for (i in 1:number_of_files){
    data[[i]] <- data[[i]][kept_events[[i]],]
  }
  return(data)
}




#' find_gate_second_top, find the gaussian gates of the second top for a vector xx
#' @param xx, vector of numbers 
#' @param lower_gate_prop, propotions for lower gate
#' @param upper_gate_prop, propotions for upper gate
#' @return list of lower and upper gates 

# find_gate_second_top <- function(xx, lower_gate_prop, upper_gate_prop, perc_included){
#   dens <- density(xx)
#   ts_y<-ts(smooth(dens$y))
#   tp <- pastecs::turnpoints(ts_y)
#   bunn1 <- dens$x[tp$pits][1]
#   xx[xx < bunn1] <- NA
#   dens <- density(xx[!is.na(xx)])  
#   lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop])
#   upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop])
#   return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
# }




# 
# 
# find_gate_second_top <- function(xx, lower_gate_prop, upper_gate_prop, perc_included, main_top_to_left = F){
# #  browser()
#   dens <- density(xx)
#   ts_y<-ts(smooth(dens$y))
#   tp <- pastecs::turnpoints(ts_y)
#   bunn1 <- dens$x[tp$pits][1]
#   xx[xx < bunn1] <- NA
#   dens <- density(xx[!is.na(xx)])
#   if(!is.na(lower_gate_prop[1])){
#     lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop])
#     upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop])
#   } else{
# 
#       xx <- xx[!is.na(xx)]
#       fit <- density(xx)
#       x.new <- rnorm(10000, sample(xx, size = 10000, replace = TRUE), fit$bw)
#        #plot(fit)
#          #lines(density(x.new), col = "red")
#          fit <- VGAM::vglm(x.new ~ 1,
#                           VGAM::mix2normal(eq.sd = F), iphi=0.5, imu= 0, isd1=2, imu2=5,
#                           isd2=1)
#        #  fit2 <- vglm(x.new ~ 1, uninormal(), lmean = 0, lsd = 1)
# 
#            # Calculated parameters
#            pars <- as.vector(coef(fit))
#            w <- VGAM::logitlink(pars[1], inverse=TRUE)
#            m1 <- pars[2]
#            sd1 <- exp(pars[3])
#            m2 <- pars[4]
#            sd2 <- exp(pars[5])
#            if(main_top_to_left == TRUE){
#              lower_gate <- m1 - qnorm(perc_included) * sd1
#              upper_gate <- m1 + qnorm(perc_included) * sd1
#            } else {
#              lower_gate <- m2 - qnorm(perc_included) * sd2
#              upper_gate <- m2 + qnorm(perc_included) * sd2
# 
#            }
# 
#   }
#   return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
# }
# 
# 
# 
# 
# 
# 













vgam_gates <- function(xx, lower_gate_prop, upper_gate_prop, perc_included, main_top_to_left){
  xx <- xx[!is.na(xx)]
  dens <- density(xx)
  x_max_dens <- dens$x[dens$y == max(dens$y)]
  x.new <- rnorm(10000, sample(xx, size = 10000, replace = TRUE), dens$bw)
  fit <- try(VGAM::vglm(x.new ~ 1,
                        VGAM::mix2normal(eq.sd = F), iphi=0.5, imu= 0, isd1=2, imu2=5,
                        isd2=1, epsilon = 1e-5))
  # Calculated parameters
  pars <- as.vector(coef(fit))
  w <- VGAM::logitlink(pars[1], inverse=TRUE)
  m1 <- pars[2]
  sd1 <- exp(pars[3])
  m2 <- pars[4]
  sd2 <- exp(pars[5])
  if(main_top_to_left == TRUE){
    lower_gate <- m1 - qnorm(perc_included) * sd1
    upper_gate <- m1 + qnorm(perc_included) * sd1
    if(lower_gate > x_max_dens | upper_gate < x_max_dens){
      lower_gate <- m2 - qnorm(perc_included) * sd2
      upper_gate <- m2 + qnorm(perc_included) * sd2
    }
  } else {
    lower_gate <- m2 - qnorm(perc_included) * sd2
    upper_gate <- m2 + qnorm(perc_included) * sd2
    if(lower_gate > x_max_dens | upper_gate < x_max_dens){
      lower_gate <- m1 - qnorm(perc_included) * sd1
      upper_gate <- m1 + qnorm(perc_included) * sd1
    }
  }
  if(lower_gate < 0 & upper_gate > 10 * x_max_dens){ #if only one normal dist. 
    temp <- find_gate_highest_top(xx, lower_gate_prop = 0.06, upper_gate_prop = 0.06)
    lower_gate <- temp$lower_gate
    upper_gate <- temp$upper_gate
  }
  
  return(list(lower_gate = lower_gate, upper_gate = upper_gate))
}


gaus_gates <- function(dens,  lower_gate_prop, upper_gate_prop){
  lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop])
  upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop])
  return(list(lower_gate = lower_gate, upper_gate = upper_gate))
  
}




find_gate_second_top <- function(xx, lower_gate_prop, upper_gate_prop, perc_included, main_top_to_left = F, minimum = NA){
  #browser()
  dens <- density(xx)
  ts_y<-ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  #bunn1 <- dens$x[tp$pits][1]
  # xx[xx < bunn1] <- NA
  # dens <- density(xx[!is.na(xx)])
  if(!is.na(minimum)){
    mini <- which(dens$x > minimum)[1]
  } else {
    mini <- 0
  }
  bunn_n <- min(max(which(tp$pits)[1], mini),  length(dens$y) - 1, na.rm = T)
 # print(bunn_n)
 ###slett# top_n <-  bunn_n + min(which(tp$peaks[bunn_n:length(dens$y)])[1], length(bunn_n:length(dens$y)), na.rm = T) - 1
  top_n <-  bunn_n + min(which(dens$y[bunn_n:length(dens$y)] == max(dens$y[bunn_n:length(dens$y)]))[1],  length(bunn_n:length(dens$y)), na.rm = T) - 1 
  #tp$peaks[top]
  h <- dens$y[top_n] - dens$y[bunn_n]
  if(!is.na(lower_gate_prop[1])){
    cutoff_h_lower <- dens$y[bunn_n] + lower_gate_prop*h
    cutoff_n_lower <- bunn_n + which(dens$y[bunn_n:length(dens$y)] > cutoff_h_lower)[1] - 1
    lower_gate <- min(dens$x[cutoff_n_lower])
    upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop])
  } else{
    # temp <- try(vgam_gates(xx, perc_included = perc_included, main_top_to_left = main_top_to_left))
    # if(inherits(temp, "try-error")){
    #   temp <- gaus_gates(xx = xx, lower_gate_prop = lower_gate_prop, upper_gate_prop =  upper_gate_prop)
    # }
    
    temp <- vgam_gates(xx, lower_gate_prop = lower_gate_prop, upper_gate_prop =  upper_gate_prop, perc_included = perc_included, main_top_to_left = main_top_to_left)
    lower_gate <- temp$lower_gate
    upper_gate <- temp$upper_gate
    
  }
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}






#' find_gaussian_gates_second_top_selected_cells, find the split between the first and secound top for all subsets
#' @param data, data 
#' @param data, data 
#' @param channel, which channel to plot
#' @param include, list with matrices over which cells to include'
#' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' @param mark, which column from include to use
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset) 
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)  
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_second_top_top_selected_cells <- function(data, channel, lower_gate_percent = NA, upper_gate_percent = NA, perc_included, main_top_to_left ,  minimum = 0, include, mark, positiv = TRUE){
# browser()
   column <- which(colnames(data[[1]]) == channel)
  if(!is.na(lower_gate_percent[[1]])){
    if(lower_gate_percent >= 1){
      lower_gate_prop <- lower_gate_percent/100
    } else {
      lower_gate_prop <- lower_gate_percent
    }
    if(upper_gate_percent >= 1){
      upper_gate_prop <- upper_gate_percent/100
    } else {
      upper_gate_prop <- upper_gate_percent
    }
  } else {
    lower_gate_prop <- NA
    upper_gate_prop <- NA
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
 for(i in 1:length(data)){
    if(positiv){
     res <-  find_gate_second_top(xx = data[[i]][include[[i]][,mark], column], lower_gate_prop = lower_gate_prop, upper_gate_prop = upper_gate_prop, perc_included = perc_included, main_top_to_left = main_top_to_left, minimum = minimum)
    } else {
      res <-  find_gate_second_top(xx = data[[i]][!include[[i]][,mark], column], lower_gate_prop = lower_gate_prop, upper_gate_prop = upper_gate_prop, perc_included = perc_included, main_top_to_left = main_top_to_left, minimum = minimum)
    }
     lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gate = lower_gates, upper_gate = upper_gates))
}



#' find_gaussian_gates_second_top, find the gaussian gates of the second top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset) 
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)  
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_second_top <- function(data, channel, lower_gate_percent = NA, upper_gate_percent = NA, perc_included, main_top_to_left , minimum = NA){
  
  column <- which(colnames(data[[1]]) == channel)
  if(!is.na(lower_gate_percent[[1]])){
    if(lower_gate_percent >= 1){
      lower_gate_prop <- lower_gate_percent/100
    } else {
      lower_gate_prop <- lower_gate_percent
    }
    if(upper_gate_percent >= 1){
      upper_gate_prop <- upper_gate_percent/100
    } else {
      upper_gate_prop <- upper_gate_percent
    }
  } else {
    lower_gate_prop <- NA
    upper_gate_prop <- NA
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
  for(i in 1:length(data)){
    res <-  find_gate_second_top(xx = data[[i]][, column], lower_gate_prop = lower_gate_prop, upper_gate_prop = upper_gate_prop, perc_included = perc_included, main_top_to_left = main_top_to_left, minimum = minimum)
    
    
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates))
}




#' find_split, find the first bottom
#' @param xx, vector of numbers 
#' @return vector of splits

find_split <- function(xx, minimum){
  dens <- density(xx)
  ts_y <- ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  bunn1 <- min(dens$x[tp$pits][dens$x[tp$pits] > minimum])
  return(bunn1)
}

#' find_split_first_second_top, find the split between the first and secound top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @return vector of splits

find_split_first_second_top <- function(data, channel, minimum = 0){
  column <- which(colnames(data[[1]]) == channel)
  splits <- rep(NA, length(data))
  for(i in 1:length(data)){
    splits[[i]] <-  find_split(data[[i]][, column], minimum)
  }
  return(splits)
}






#' find_split_first_second_top_selected_cells, find the split between the first and secound top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @param include, list with matrices over which cells to include'
#' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' @param mark, which column from include to use
#' @return vector of splits

find_split_first_second_top_selected_cells <- function(data, channel, minimum = 0, include, mark, positiv = TRUE){
  column <- which(colnames(data[[1]]) == channel)
  splits <- rep(NA, length(data))
  for(i in 1:length(data)){
    if(positiv){
      splits[[i]] <-  find_split(data[[i]][include[[i]][,mark], column], minimum)
    } else {
      splits[[i]] <-  find_split(data[[i]][!include[[i]][,mark], column], minimum)
    }
  }
  return(splits)
}



#' find_split_first_second_top, find the split between the first and secound top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @return vector of splits

find_split_neg_low_high <- function(data, channel, neg = 0.05, minLowHigh = 0.1){
    column <- which(colnames(data[[1]]) == channel)
    splits <- rep(NA, length(data))
    neg <- rep(NA, length(data))
    lower_gates <- splits
    upper_gates <- splits
    for(i in 1:length(data)){
      xx <- data[[i]][,column]
      neg[[i]] <- find_split(xx, minimum = 0)
      #     neg[[i]] <- find_gate_first_top(xx, lower_gate_prop = negProp, upper_gate_prop = negProp)$upper_gate
      splits[[i]] <-  find_split(xx[xx > neg[[i]]], minimum = minLowHigh)
      
    }
    
    return(list(neg_splits = neg, low_high_splits = splits))
    
    
    
  }
  
  


#' find_split_first_second_top, find the split between the first and secound top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @param include, list with matrices over which cells to include'
#' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' @param mark, which column from include to use

#' @return vector of splits

find_split_neg_low_high_selected_cells <- function(data, channel, include, mark, positiv = TRUE, neg = 0.05, minLowHigh = 0.1){
  column <- which(colnames(data[[1]]) == channel)
  splits <- rep(NA, length(data))
  neg <- rep(NA, length(data))
  lower_gates <- splits
  upper_gates <- splits
  for(i in 1:length(data)){
    if(positiv){
      xx <- data[[i]][include[[i]][,mark],column]
    } else{
      xx <- data[[i]][!include[[i]][,mark],column]
    }
    neg[[i]] <- find_split(xx, minimum = 0)
    #    neg[[i]] <- find_gate_first_top(xx, lower_gate_prop = negProp, upper_gate_prop = negProp)$upper_gate
        splits[[i]] <-  find_split(xx[xx > neg[[i]]], minimum = minLowHigh)
  }
  return(list(neg_splits = neg, low_high_splits = splits))
}



#' find_gate_first_top, find the gaussian gates of the first top for a vector xx
#' @param xx, vector of numbers 
#' @param lower_gate_prop, propotions for lower gate
#' @param upper_gate_prop, propotions for upper gate
#' @return list of lower and upper gates 

find_gate_first_top <- function(xx, lower_gate_prop, upper_gate_prop, min_upper_gate = NA){
  dens <- density(xx)
  ts_y <- ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  bunn1 <- dens$x[tp$pits][1]
  xx[xx > bunn1] <- NA
  dens <- density(xx[!is.na(xx)])  
  lower_gate <- max(min(dens$x[dens$y > max(dens$y) * lower_gate_prop]))
  upper_gate <- max(c(dens$x[dens$y > max(dens$y) * upper_gate_prop], min_upper_gate), na.rm = T)
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}




#' find_gaussian_gates_first_top, find the gaussian gates of the first top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset) 
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)  
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_first_top <- function(data, channel, lower_gate_percent, upper_gate_percent, min_upper_gate = NA){
  column <- which(colnames(data[[1]]) == channel)
  if(lower_gate_percent > 1){
    lower_gate_prop <- lower_gate_percent/100
  } else {
    lower_gate_prop <- lower_gate_percent
  }
  if(upper_gate_percent > 1){
    upper_gate_prop <- upper_gate_percent/100
  } else {
    upper_gate_prop <- upper_gate_percent
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
  for(i in 1:length(data)){
    res <-  find_gate_first_top(data[[i]][, column], lower_gate_prop = lower_gate_prop, upper_gate_prop = upper_gate_prop, min_upper_gate = min_upper_gate)
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates ))
}






#' find_gate_highest_top, find the gaussian gates of the main top for a vector xx
#' @param xx, vector of numbers 
#' @param lower_gate_prop, propotions for lower gate
#' @param upper_gate_prop, propotions for upper gate
#' @return list of lower and upper gates 

find_gate_highest_top <- function(xx, lower_gate_prop, upper_gate_prop, min_upper_gate = NA){
  dens <- density(xx)
  ts_y <- ts(smooth(dens$y))
  lower_gate <- max(min(dens$x[dens$y > max(dens$y) * lower_gate_prop]))
  upper_gate <- max(c(dens$x[dens$y > max(dens$y) * upper_gate_prop], min_upper_gate), na.rm = T)
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}




#' find_gaussian_gates_highest_top, find the gaussian gates of the highest top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset) 
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)  
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_highest_top <- function(data, channel, lower_gate_percent, upper_gate_percent, min_upper_gate = NA){
  column <- which(colnames(data[[1]]) == channel)
  if(lower_gate_percent > 1){
    lower_gate_prop <- lower_gate_percent/100
  } else {
    lower_gate_prop <- lower_gate_percent
  }
  if(upper_gate_percent > 1){
    upper_gate_prop <- upper_gate_percent/100
  } else {
    upper_gate_prop <- upper_gate_percent
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
  for(i in 1:length(data)){
    res <-  find_gate_highest_top(data[[i]][, column], lower_gate_prop = lower_gate_prop, upper_gate_prop = upper_gate_prop, min_upper_gate = min_upper_gate)
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates ))
}



      


