#' prosent_senario
#' @param posneg, result per senario from channel gating 
#' @param channels, which channels to select
#' @param values, the values for the channels in the same order, 0 = neg, 1 = pos (or low), 2 = high, 12 = pos (high and low), 10 = neg og low  
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
      if(!(values[j] == 12 | values[j] == 10)){
        if(values[j] %in% c(0,1,2)){
          x <- posneg[[i]][, column_j] == values[j]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(values[j] == 12){
          x <- posneg[[i]][, column_j] %in% c(1,2)
        } else {
          if(values[j] == 10){
            x <- posneg[[i]][, column_j] %in% c(0,1)
          }
        }
      }  
      xx <- x & xx
    }
    xTrue <- max(c(table(xx)["TRUE"],0), na.rm = T)
    res[i] <- xTrue/length(xx) * 100
  }
  return(res)
}


any_value <- function(x, value){
  res <- 0
  if(value == 1){
    res <- any(x == 1)
  } else {
    if(value == 0){
      res <- any(x == 0)
    } else  {
      if(value == 12){
        res <- any(x %in% c(1,2))
      } else {
        if(value == 10){
          res <- any(x %in% c(0, 1))
        }
      }
    }
  }
  return(res)
}


#' prosent_senario_with_atleat_one_of_some_channels
#' @param posneg, result per senario from channel gating 
#' @param channels, which channels to select
#' @param values, the values for the channels in the same order, 0 = neg, 1 = pos (or low), 2 = high, 12 = pos (high and low), 10 = neg og low  
#' @return vector with prosent that have the asked combination.


prosent_senario_with_atleat_one_of_some_channels <- function(posneg, channels, values, atleast_one_of, value_atleast_one_of = 1){
  if(length(values) == 1){
    values <- rep(values, length(channels))
  }
  res <- rep(NA, length(posneg))
  for(i in 1:length(posneg)){
    columns <-  which(colnames(posneg[[i]]) %in% atleast_one_of)
      xx <- apply(posneg[[i]][, columns], 1, any_value, value = value_atleast_one_of)
    
    for(j in 1:length(channels)){
      column_j <- which(colnames(posneg[[i]]) == channels[j]) 
      if(!(values[j] == 12 | values[j] == 10)){
        if(values[j] %in% c(0,1,2)){
          x <- posneg[[i]][, column_j] == values[j]
        } else {
          browser()
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(values[j] == 12){
          x <- posneg[[i]][, column_j] %in% c(1,2)
        } else {
          if(values[j] == 10){
            x <- posneg[[i]][, column_j] %in% c(0,1)
          }
        }
      }  
      xx <- x & xx
    }
    xTrue <- max(c(table(xx)["TRUE"],0), na.rm = T)
    res[i] <- xTrue/length(xx) * 100
  }
  return(res)
}



#' prosent_per_channel
#' @param posneg, result per senario from channel gating 
#' @param channels, which channels to select
#' @param values, the values for the channels in the same order, 0 = neg, 1 = pos (or low), 2 = high, 12 = pos (high and low), 10 = neg (neg and low)  
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
  #browser()
  for(i in 1:length(channels)){
    column_i <- which(colnames(posneg[[1]]) == channels[i]) 
    
    for(j in 1:length(posneg)){
      
      if(!(values[i] == 12 | values[i] == 10)){
        if(values[i] %in% c(0,1,2)){
          
          xx <- posneg[[j]][, column_i] == values[i]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(values[i] == 12){
          xx <- posneg[[j]][, column_i] %in% c(1,2)
          
        } else {
          if(values[i] == 10){
            xx <- posneg[[j]][, column_i] %in% c(0,1)
          }
        }
      } 
      mat[j,i] <- table(xx)["TRUE"]/length(xx) * 100
      
      #   print(j)
    }
    
  }
  
  return(mat)
}



senario_kolonne <- function(data, senario_channels, senario_values){
  if(length(senario_values) == 1){
    senario_values <- rep(senario_values, length(senario_channels))
  }
  d <- data
  xx <- rep(1, nrow(d))
  for(senario_j in 1:length(senario_channels)){
    column_j <- which(colnames(data) == senario_channels[senario_j]) 
    if(!(senario_values[senario_j] == 12 | senario_values[senario_j] == 10)){
      if(senario_values[senario_j] %in% c(0,1,2)){
        x <- d[, column_j] == senario_values[senario_j]
      } else {
        print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
      }
    } else{
      if(senario_values[senario_j] == 12){
        x <- d[, column_j] %in% c(1,2)
      } else {
        if(senario_values[senario_j] == 10){
          x <- d[, column_j] %in% c(0,1)
        }
      }
    }  
    xx <- x & xx
  }
  return(xx)
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
      column_j <- which(colnames(posneg) == senario_channels[senario_j]) 
      if(!(senario_values[senario_j] == 12 | senario_values[senario_j] == 10)){
        if(senario_values[senario_j] %in% c(0,1,2)){
          x <- d[, column_j] == senario_values[senario_j]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(senario_values[senario_j] == 12){
          x <- d[, column_j] %in% c(1,2)
        } else {
          if(senario_values[senario_j] == 10){
            x <- d[, column_j] %in% c(0,1)
          }
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




cells_to_select_from <- function(file_names, posNeg, senario_channels, senario_values){
  result <- NULL
  if(length(file_names) != length(posNeg)){
    stop("Need same length of filnames and posNeg")
  }
  
  for(i in 1:length(file_names)){
    
    d <- posNeg[[i]]
    xx <- rep(TRUE, nrow(d))
    # browser()    
    for(senario_j in 1:length(senario_channels)){
      column_j <- which(colnames(posNeg[[i]]) == senario_channels[senario_j]) 
      if(!(senario_values[senario_j] == 12 | senario_values[senario_j] == 10)){
        if(senario_values[senario_j] %in% c(0,1,2)){
          x <- d[, column_j] == senario_values[senario_j]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(senario_values[senario_j] == 12){
          x <- d[, column_j] %in% c(1,2)
        } else {
          if(senario_values[senario_j] == 10){
            x <- d[, column_j] %in% c(0,1)
          }
        }
      }  
      xx <- x & xx
    }
    result[[i]] <- which(xx)
  }
  
  return(result)
}


# # må lages
# random_events_from_selected_cells <- function(mulige_celler, marker = NA, n = 25000){
#   result <- NULL
#  # browser()
#   for(i in 1:length(mulige_celler)){
#    if(is.na(marker)){
#       result[[i]] <- sample(mulige_celler[[i]], size = min(n, length(mulige_celler[[i]])), replace = FALSE)
#     } else {
#       result[[i]] <- ssample(mulige_celler[[i]][,marker], size = min(n, sum(mulige_celler[[i]][,marker])), replace = FALSE)
#     }
#   }
#   
#   
#   return(result)
# }




median_per_cluster_marker<- function(data, kluster){
  result <- as.data.frame(matrix(NA, ncol = ncol(data), nrow = length(unique(kluster))))
  colnames(result) <- colnames(data)
  rownames(result) <- sort(unique(kluster))
  for(i in rownames(result)){
    result[i,] <- apply(data[as.character(kluster) %in% i,], 2, median)
  }
  return(result)
}





q_per_cluster_marker<- function(data, kluster, probs){
  if(probs > 1){
    probs <- probs/100
  }
  result <- as.data.frame(matrix(NA, ncol = ncol(data), nrow = length(unique(kluster))))
  colnames(result) <- colnames(data)
  rownames(result) <- sort(unique(kluster))
  for(i in rownames(result)){
    result[i,] <- apply(data[as.character(kluster) %in% i,], 2, quantile, probs = probs)
  }
  return(result)
}



extra_column_posneg <- function(posNeg, column_name, senario_channels, senario_values){
  for(i in 1:length(posNeg)){
    kol <- ncol(posNeg[[i]])
    posNeg[[i]][,kol + 1] <- rep(TRUE, nrow(posNeg[[i]]))
    for(senario_j in 1:length(senario_channels)){
      column_j <- which(colnames(posNeg[[i]]) == senario_channels[senario_j]) 
      if(!(senario_values[senario_j] == 12 | senario_values[senario_j] == 10)){
        if(senario_values[senario_j] %in% c(0,1,2)){
          x <- posNeg[[i]][, column_j] == senario_values[senario_j]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(senario_values[senario_j] == 12){
          x <- posNeg[[i]][, column_j] %in% c(1,2)
        } else {
          if(senario_values[senario_j] == 10){
            x <- posNeg[[i]][, column_j] %in% c(0,1)
          }
        }
      }  
      posNeg[[i]][,kol + 1] <- x &  posNeg[[i]][,kol + 1]
    }
    colnames(posNeg[[i]])[kol + 1] <- column_name
  }
  return(posNeg)
}
