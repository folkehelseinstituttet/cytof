#' Read all fcs files in folder
#' 
#' @param path a path to the folder with fcs files
#' @return fcs_data a set of all data produced. 

read_data_from_folder <- function(data_path){
  fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
  files_to_open <- basename(fcs_files)
  files_to_open <- files_to_open[grepl(".fcs", files_to_open)]
  setwd(dirname(fcs_files[1]))
  file_names <- gsub(".fcs", "", files_to_open)
  # Read the files into a flowset
  fcs_data <- flowCore::read.flowSet(files_to_open, transformation = FALSE,
                                     truncate_max_range = FALSE)
  return(list(fcs_data = fcs_data, file_names = file_names))
}


#' Read some fcs files in folder
#' 
#' @param path a path to the folder with fcs files
#' @return fcs_data a set of all data produced. 

read_some_data_from_folder <- function(data_path, file_number = 1:4){
  fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
  files_to_open <- basename(fcs_files)
  files_to_open <- files_to_open[grepl(".fcs", files_to_open)]
  files_to_open <- files_to_open[file_number]
  setwd(dirname(fcs_files[1]))
  file_names <- gsub(".fcs", "", files_to_open)
  # Read the files into a flowset
  fcs_data <- flowCore::read.flowSet(files_to_open, transformation = FALSE,
                                     truncate_max_range = FALSE)
  return(list(fcs_data = fcs_data, file_names = file_names))
}


#' Read specific fcs files in folder
#' 
#' @param path a path to the folder with fcs files
#' @return fcs_data a set of all data produced. 

read_specific_data_from_folder <- function(data_path, files_to_open){
  setwd(data_path)
  file_names <- gsub(".fcs", "", files_to_open)
  # Read the files into a flowset
  fcs_data <- flowCore::read.flowSet(files_to_open, transformation = FALSE,
                                     truncate_max_range = FALSE)
  return(list(fcs_data = fcs_data, file_names = file_names))
}




#' get all parameters in fcs_data with info
#' 
#' @param x one file in the fcs_data, e.g. fcs_data[[1]]
#' @return the parameters in fcs_data.

get_params_fcs_data <- function(x = fcs_data[[1]]){
  params <- flowCore::pData(flowCore::parameters(x))
  return(params)
}




list_to_matrix <- function(data, file_names = NA){
  n <- length(data)
  if(is.na(file_names[1])){
    file_names <- 1:n
  }
  
  if(!length(file_names) == n){
    print("There has to be equal amount of file_names and datasets")
    stop()
  }
  
  mat <- data[[1]]
  mat$dataset <- file_names[1]
  if(n > 1){
    for(i in 2:n){
      mat0 <- data[[i]]
      mat0$dataset <- file_names[i]
      mat <- rbind(mat, mat0)
    }
  }
  
  mat <- as.data.frame(mat)
  return(mat)
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

#' list_to_matrix_selected_events
#' @param data, data 
#' @param kept_events, list of vectors of true/false for alle events in each subdataset. Those events that are true will be kept
#' @param file_names, filnavn til kolonne, hvis na nummerert filnavn
#' @param channels, kanaler som skal være med
#' @param archSin, default true
#' @param cofactor, default 5
#' @param scale, default F
#' @param scaling, default c(-5,12000), only used when scale = T
#' @param cytofData, default True
#' @return new dataset with only those events that we want to keep. 

list_to_matrix_selected_events <- function(data, kept_events, file_names = NA, channels, archSin = T, 
                                           cofactor = 5, scale = F,  scaling = c(-5, 12000), cytofData = T){
  #  browser()
  n <- length(data)
  if(is.na(file_names[1])){
    file_names <- 1:n
  }
  
  if(!length(file_names) == n){
    print("There has to be equal amount of file_names and datasets")
    stop()
  }
  
  if(cytofData == F){
    mat <- as.data.frame(data[[1]][kept_events[[1]], channels])
    mat$dataset <- file_names[1]
    if(n > 1){
      for(i in 2:n){
        mat0 <- as.data.frame(data[[i]][kept_events[[i]], channels])
        mat0$dataset <- file_names[i]
        mat <- rbind(mat, mat0)
      }
    }
  } else {
    
    
    if(archSin == T){
      if(scale == F){
        mat <- as.data.frame(asinh(flowCore::exprs(data[[1]][kept_events[[1]], channels])/cofactor))
        mat$dataset <- file_names[1]
        mat$Time <- flowCore::exprs(fcs_data[[1]][kept_events[[1]],"Time"])
        if(n > 1){
          for(i in 2:n){
            mat0 <- as.data.frame(asinh(flowCore::exprs(data[[i]][kept_events[[i]], channels])/cofactor))
            mat0$dataset <- file_names[i]
            mat0$Time <- flowCore::exprs(fcs_data[[i]][kept_events[[i]],"Time"])
            mat <- rbind(mat, mat0)
          }
        }
      } else {  #rescale(new_data[[i]][,j],to = scaling)
        mat <- as.data.frame(rescale(asinh(flowCore::exprs(data[[1]][kept_events[[1]], channels])/cofactor),to = scaling))
        mat$dataset <- file_names[1]
        mat$Time <- flowCore::exprs(fcs_data[[1]][,"Time"])
        if(n > 1){
          for(i in 2:n){
            mat0 <- as.data.frame(rescale(asinh(flowCore::exprs(data[[i]][kept_events[[i]], channels])/cofactor),to = scaling))
            mat0$dataset <- file_names[i]
            mat0$Time <- flowCore::exprs(fcs_data[[i]][,"Time"])
            mat <- rbind(mat, mat0)
          }
        }
      }
      
    } else {
      mat <- as.data.frame(flowCore::exprs(data[[1]][kept_events[[1]], channels]))
      mat$dataset <- file_names[1]
      mat$Time <- flowCore::exprs(fcs_data[[1]][,"Time"])
      if(n > 1){
        for(i in 2:n){
          mat0 <- as.data.frame(flowCore::exprs(data[[i]][kept_events[[i]], channels]))
          mat0$dataset <- file_names[i]
          mat0$Time <- flowCore::exprs(fcs_data[[i]][,"Time"])
          mat <- rbind(mat, mat0)
        }
      }
    }
    
    
  }
  mat <- as.data.frame(mat)
  return(mat)
  
}
