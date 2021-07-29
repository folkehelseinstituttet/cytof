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
  fit <- vglm(x.new ~ 1, 
             mix2normal(eq.sd = F), iphi=0.5, imu= 0, isd1=2, imu2=5, 
             isd2=1)
#  fit2 <- vglm(x.new ~ 1, uninormal(), lmean = 0, lsd = 1)
  
  # Calculated parameters
  pars <- as.vector(coef(fit))
  w <- logitlink(pars[1], inverse=TRUE)
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




#' cells_to_keep, find which cells to keep in each subdataset based on lower and/or upper gate.
#' @param data, transformed data 
#' @param channel, which channel to plot
#' @param lower_gate, vector with values for lower gate, a number (same lower gate for all subset) or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate, a number (same upper gate for all subset)  or NA (no upper gating)
#' @return list of vectors of true/false for each cells (true means keep), one vector for each subdataset. 


cells_to_keep <- function(data, channel, lower_gate = NA, upper_gate = NA){
  number_of_files <- length(data)
  column <- which(colnames(data[[1]]) == channel)
  if(length(lower_gate) == 1){
    lower_gate <- rep(lower_gate, number_of_files)
  }
  if(length(upper_gate) == 1){
    upper_gate <- rep(upper_gate, number_of_files)
  }
  kept_cells <- NULL
  if(is.na(lower_gate[1])){
    for (i in 1:number_of_files){
      kept_cells[[i]] <- data[[i]][,column] < upper_gate[i]
    }
  } else { 
    if(is.na(upper_gate[1])){
      for (i in 1:number_of_files){
        kept_cells[[i]] <- data[[i]][,column] > lower_gate[i]
      }
    } else {     
      for (i in 1:number_of_files){
        kept_cells[[i]] <- data[[i]][,column] < upper_gate[i] & data[[i]][,column] > lower_gate[i]
      }
    }
  }  
  return(kept_cells)
}


#' update_data_based_on_cells_to_keep
#' @param data, data 
#' @param kept_cells, list of vectors of true/false for alle cells in each subdataset. Those cells that are true will be kept
#' @return new dataset with only those cells that we want to keep. 

update_data_based_on_cells_to_keep <- function(data, kept_cells){
  number_of_files <- length(data)
  for (i in 1:number_of_files){
    data[[i]] <- data[[i]][kept_cells[[i]],]
  }
  return(data)
}




#' find_gate_second_top, find the gaussian gates of the second top for a vector xx
#' @param xx, vector of numbers 
#' @param lower_gate_prop, propotions for lower gate
#' @param upper_gate_prop, propotions for upper gate
#' @return list of lower and upper gates 

find_gate_second_top <- function(xx, lower_gate_prop, upper_gate_prop){
  dens <- density(xx)
  ts_y<-ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  bunn1 <- dens$x[tp$pits][1]
  xx[xx < bunn1] <- NA
  dens <- density(xx[!is.na(xx)])  
  lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop])
  upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop])
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}


#' find_gaussian_gates_second_top, find the gaussian gates of the second top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset) 
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)  
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_second_top <- function(data, channel, lower_gate_percent, upper_gate_percent){
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
    res <-  find_gate_second_top(xx = data[[i]][, column], lower_gate_prop = lower_gate_prop, upper_gate_prop = upper_gate_prop)
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates))
}




#' find_split_first_second_top, find the first bottom
#' @param xx, vector of numbers 
#' @return vector of splits

find_split <- function(xx){
  dens <- density(xx)
  ts_y <- ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  bunn1 <- dens$x[tp$pits][1]
  return(bunn1)
}

#' find_gaussian_gates_first_top, find the gaussian gates of the first top for all subsets
#' @param data, data 
#' @param channel, which channel to plot
#' @return vector of splits

find_split_first_second_top <- function(data, channel){
  column <- which(colnames(data[[1]]) == channel)
  splits <- rep(NA, length(data))
  for(i in 1:length(data)){
    splits[[i]] <-  find_split(data[[i]][, column])
  }
  return(splits)
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




