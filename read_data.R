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


#' get all parameters in fcs_data with info
#' 
#' @param x one file in the fcs_data, e.g. fcs_data[[1]]
#' @return the parameters in fcs_data.

get_params_fcs_data <- function(x = fcs_data[[1]]){
  params <- flowCore::pData(flowCore::parameters(x))
  return(params)
}
