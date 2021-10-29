#' arc_sinh_transform_selected_channels transform data with arc_sinh
#' @params fcs_data a set of fcs files
#' @params channels transforms only those channels listed
#' @params cofactor can be changed from 5 that is default
#' @params scale, if scaling wanted this have to be true, default false
#' @params scaling, interval to scale each parameter to, only used if scale = true
#' @return scaled data for those channels chosen

arc_sinh_transform_selected_channels <- function(fcs_data, channels, cofactor = 5, scale = F,  scaling = c(-5, 12000)){
  new_data <- NULL 
  number_of_files <- length(fcs_data)
  
  for (i in 1:number_of_files){
    new_data[[i]] <- as.data.frame(asinh(flowCore::exprs(fcs_data[[i]][,channels])/cofactor))
    if(scale){
      for (j in 1:length(channels)){
        new_data[[i]][,j] <- rescale(new_data[[i]][,j],to = scaling)
      }
    }
  }
  # Add Time
  for (i in 1:number_of_files){
    new_data[[i]]$Time <- flowCore::exprs(fcs_data[[i]][,"Time"])
  }
  return(new_data)
}




#' log_transform_selected_channels transform data with log
#' @params fcs_data a set of fcs files
#' @params channels transforms only those channels listed
#' @params scale, if scaling wanted this have to be true, default false
#' @params scaling, interval to scale each parameter to, only used if scale = true
#' @return scaled data for those channels chosen

log_transform_selected_channels <- function(fcs_data, channels, scale = F,  scaling = c(-5, 12000)){
  new_data <- NULL 
  number_of_files <- length(fcs_data)
  
  for (i in 1:number_of_files){
    new_data[[i]] <- as.data.frame(log(flowCore::exprs(fcs_data[[i]][,channels]) + 1))
    if(scale){
      for (j in 1:length(channels)){
        new_data[[i]][,j] <- rescale(new_data[[i]][,j],to = scaling)
      }
    }
  }
  # Add Time
  for (i in 1:number_of_files){
    new_data[[i]]$Time <- flowCore::exprs(fcs_data[[i]][,"Time"])
  }
  return(new_data)
}

