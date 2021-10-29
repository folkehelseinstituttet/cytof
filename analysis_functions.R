#' prosent_senario
#' @param posneg, result per senario from channel gating 
#' @param channels, which channels to select
#' @param values, the values for the channels in the same order, 0 = neg, 1 = pos (or low), 2 = high, 3 = pos (high and low), -1 = neg og low  
#' @return vector with prosent that have the asked combination.


prosent_senario <- function(posneg, channels, values){
  if(length(values) == 1){
    values <- rep(values, length(channels))
  }
  res <- rep(NA, length(posneg))
  for(i in 1:length(posneg)){
    xx <- rep(1, nrow(posneg[[i]]))
    for(j in 1:length(channels)){
      column_j <- which(colnames(posneg[[i]]) == channels[j]) 
      if(!(values[j] == 3 | values[j] == -1)){
        if(values[j] %in% c(0,1,2)){
          x <- posneg[[i]][, column_j] == values[j]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: -1 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 3 = c(1,2")
        }
      } else{
        if(values[j] == 3){
          x <- posneg[[i]][, column_j] %in% c(1,2)
        } else {
          x <- posneg[[i]][, column_j] %in% c(0,1)
        }
      }  
      xx <- x & xx
    }
    res[i] <- table(xx)["TRUE"]/length(xx) * 100
  }
  return(res)
}


#' prosent_per_channel
#' @param posneg, result per senario from channel gating 
#' @param channels, which channels to select
#' @param values, the values for the channels in the same order, 0 = neg, 1 = pos (or low), 2 = high, 3 = pos (high and low)  
#' @return vector with prosent that have the asked combination.
prosent_per_channel <- function(posneg, channels = "all", values = 0){
  if(channels[1] == "all"){
    channels <- colnames(posneg[[1]])
  }
  if(length(values) == 1){
    values <- rep(values, length(channels))
  }
  
  mat <- matrix(NA, ncol = length(channels), nrow = length(posneg))
  colnames(mat) <- channels
  rownames(mat) <- names(posneg)
  
  
  for(i in 1:length(channels)){
    column_i <- which(colnames(posneg[[1]]) == channels[i]) 
    
    for(j in 1:length(posneg)){
      
      if(!(values[i] == 3 | values[i] == -1)){
        if(values[i] %in% c(0,1,2)){
          
          xx <- posneg[[j]][, column_i] == values[i]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: -1 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 3 = c(1,2")
        }
      } else{
        if(values[j] == 3){
          xx <- posneg[[j]][, column_i] %in% c(1,2)
        } else {
          xx <- posneg[[j]][, column_i] %in% c(0,1)
        }
      } 
      mat[j,i] <- table(xx)["TRUE"]/length(xx) * 100
    }
  }
  return(mat)
}



#' prosent_per_channel_senario
#' @param posneg, result per senario from channel gating 
#' @param channels, which channels to select
#' @param values, the values for the channels in the same order, 0 = neg, 1 = pos (or low), 2 = high, 3 = pos (high and low)  
#' @return vector with prosent that have the asked combination.
prosent_per_channel_senario <- function(posneg, channels = "all", values = 0, senario_channels, senario_values){
  if(channels[1] == "all"){
    channels <- colnames(posneg[[1]])
  }
  if(length(values) == 1){
    values <- rep(values, length(channels))
  }
  
  mat <- matrix(NA, ncol = length(channels), nrow = length(posneg))
  colnames(mat) <- channels
  rownames(mat) <- names(posneg)
  
  
  for(j in 1:length(posneg)){
    d <- posneg[[j]]
    
    xx <- rep(1, nrow(d))
    for(senario_j in 1:length(senario_channels)){
      column_j <- which(colnames(data) == senario_channels[senario_j]) 
      if(!(senario_values[senario_j] == 3 | senario_values[senario_j] == -1)){
        if(senario_values[senario_j] %in% c(0,1,2)){
          x <- d[, column_j] == senario_values[senario_j]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: -1 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 3 = c(1,2")
        }
      } else{
        if(senario_values[senario_j] == 3){
          x <- d[, column_j] %in% c(1,2)
        } else {
          x <- d[, column_j] %in% c(0,1)
        }
      }  
      xx <- x & xx
    }
    d <- d[xx,]
    
    
    for(i in 1:length(channels)){
      column_i <- which(colnames(d) == channels[i]) 
      
      
      if(!(values[i] == 3 | values[i] == -1)){
        if(values[i] %in% c(0,1,2)){
          
          xx <- d[, column_i] == values[i]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: -1 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 3 = c(1,2")
        }
      } else{
        if(values[j] == 3){
          xx <- d[, column_i] %in% c(1,2)
        } else {
          xx <- d[, column_i] %in% c(0,1)
        }
      } 
      mat[j,i] <- table(xx)["TRUE"]/length(xx) * 100
    }
  }
  return(mat)
}
