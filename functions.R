#' make_project_folders, genetrate the directory structure for a project, do not overwrite excisting folders.
#' @param path, the path to where you want to start or continue a project 
#' @return vector with prosent that have the asked combination.

make_project_folders <- function(path){
  org::initialize_project(
    home = path,
    create_folders = TRUE, 
    
    raw_data_path = fs::path(path, "raw_data"),
    R_path = fs::path(path, "R"),
    script_path = fs::path(path, "script"),
    clean_up_path = fs::path(path, "clean_up"),
    clean_data_path = fs::path(path, "clean_up", "clean_data"),
    clean_data_info_path = fs::path(path, "clean_up",  "clean_data_info"),
    clean_data_fig_density_path = fs::path(path, "clean_up", "clean_data_fig_density"),
    clean_data_fig_signal_path = fs::path(path, "clean_up", "clean_data_fig_signal"),
    marker_gating_path = fs::path(path, "clean_up", "marker_gating"),  
    marker_gating_fig_density_path = fs::path(path, "clean_up", "marker_gating",  "fig_density"),  
    marker_gating_fig_signal_path = fs::path(path, "clean_up", "marker_gating",  "fig_signal"),  
    marker_gating_results_path = fs::path(path, "clean_up", "marker_gating",  "results"),
    clean_data_flowSOM_results_path = fs::path(path, "clean_up", "flowSOM_results"),
    meta_data_path = fs::path(path, "meta_data")
  )
 paths <- list()
 paths$raw_data_path <- fs::path(path, "raw_data")
 paths$R_path <- fs::path(path, "R")
 paths$script_path <- fs::path(path, "script")
 paths$clean_up_path <- fs::path(path, "clean_up")
 paths$clean_data_path <- fs::path(path, "clean_up", "clean_data")
 paths$clean_data_fig_density_path <- fs::path(path, "clean_up", "clean_data_fig_density")
 paths$clean_data_fig_signal_path <- fs::path(path, "clean_up", "clean_data_fig_signal")
 paths$clean_data_info_path <- fs::path(path, "clean_up",  "clean_data_info")
 paths$clean_data_posNeg_path <- fs::path(path, "clean_up",  "clean_data_posNeg")
 paths$marker_gating_path <- fs::path(path, "clean_up", "marker_gating")
 paths$marker_gating_fig_density_path <- fs::path(path, "clean_up", "marker_gating",  "fig_density")
 paths$marker_gating_fig_signal_path <- fs::path(path, "clean_up", "marker_gating",  "fig_signal")
 paths$marker_gating_results_path <- fs::path(path, "clean_up", "marker_gating",  "results")
 paths$clean_data_flowSOM_results_path <- fs::path(path, "clean_up", "flowSOM_results")
 paths$meta_data_path <- fs::path(path, "meta_data")
 
 return(paths)
}




#' any_value, used in function prosent_senario
#' @param x, vector of markeres where atleast one should have the value, value
#' @param value, the value that atleast one markere should have

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


#' prosent_senario, calculate the percentages of cells that follows this senario for each sample.
#' @param posneg, result per marker from gating 
#' @param markers, which markers to select
#' @param values, the values for the markers in the same order, 0 = neg, 1 = pos (or low), 2 = high, 12 = pos (high and low), 10 = neg og low  
#' @param atleast_one_of, vector of markers where atleast one of the markers should have the value "value_atleast_one_of"
#' @param value_atleast_one_of, default 12 which means that atleast one of the markers in atleast_one_of should be 1 or 2
#' @return vector with prosent that have the asked combination.


prosent_senario <- function(posneg, markers, values, atleast_one_of = NA, value_atleast_one_of = 12){
  if(length(values) == 1){
    values <- rep(values, length(markers))
  }
  res <- rep(NA, length(posneg))
  for(i in 1:length(posneg)){
    xx <- rep(1, nrow(posneg[[i]]))
    for(j in 1:length(markers)){
      column_j <- which(colnames(posneg[[i]]) == markers[j]) 
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
    
    if(!is.na(atleast_one_of)){
      columns <-  which(colnames(posneg[[i]]) %in% atleast_one_of)
      x <- apply(posneg[[i]][, columns], 1, any_value, value = value_atleast_one_of)
      xx <- x & xx
    }
    
    
    
    xTrue <- max(c(table(xx)["TRUE"],0), na.rm = T)
    res[i] <- xTrue/length(xx) * 100
  }
  return(res)
}


#' number_of_positive_events
#' @param posNeg, result per marker from gating 
#' @param marker, which marker to calculate number of events for 
#' @return vector with number of events per file
number_of_positive_events <- function(posNeg, marker){
  number_of_events <- rep(NA, length(posNeg))
  for(i in 1:length(posNeg)){
    tab <- table(posNeg[[i]][,marker])
    number_of_events[i] <- max(tab["TRUE"], tab["1"], tab["1"] + tab["2"], 0, na.rm = T)
  }
  return(number_of_events)
}

#' prosent_per_marker, calculate the percentages of cells that follows each of the selected markers for each sample.
#' @param posneg, result per marker from gating 
#' @param markers, which markers to calculate perentages for. 
#' @param values, the values for the markers in the same order, 0 = neg, 1 = pos (or low), 2 = high, 12 = pos (high and low), 10 = neg (neg and low)  
#' @return matrix with prosent that have each of the selected markers for each sample.

prosent_per_marker <- function(posneg, markers = "all", values = 0){
  if(markers[1] == "all"){
    markers <- colnames(posneg[[1]])
  }
  if(length(values) == 1){
    values <- rep(values, length(markers))
  }
  
  mat <- matrix(NA, ncol = length(markers), nrow = length(posneg))
  colnames(mat) <- markers
  rownames(mat) <- names(posneg)
  #browser()
  for(i in 1:length(markers)){
    column_i <- which(colnames(posneg[[1]]) == markers[i]) 
    
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


#' extra_column_posneg, introduce an extra column to each matrix in posneg
#' @param posneg, result per marker from gating 
#' @param column_name, name on the extra column generated
#' @param markers, which markers to combine to obtain this column 
#' @param values, the values for the markers in the same order, 0 = neg, 1 = pos (or low), 2 = high, 12 = pos (high and low), 10 = neg (neg and low)  
#' @return matrix with prosent that have each of the selected markers for each sample.


extra_column_posneg <- function(posNeg, column_name, markers, markers_values){
    for(i in 1:length(posNeg)){
    kol <- ncol(posNeg[[i]])
    posNeg[[i]][,kol + 1] <- rep(TRUE, nrow(posNeg[[i]]))
    for(markersj in 1:length(markers)){
      column_j <- which(colnames(posNeg[[i]]) == markers[markersj]) 
      if(!(markers_values[markersj] == 12 | markers_values[markersj] == 10)){
        if(markers_values[markersj] %in% c(0,1,2)){
          x <- posNeg[[i]][, column_j] == markers_values[markersj]
        } else {
          print("ugyldig valg av verdi, gyldige verdier er: 10 = c(0,1), 0 = 0, 1 = 1, 2 = 2, 12 = c(1,2)")
        }
      } else{
        if(markers_values[markersj] == 12){
          x <- posNeg[[i]][, column_j] %in% c(1,2)
        } else {
          if(markers_values[markersj] == 10){
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



#' q_per_cluster_marker, calucate the q percentage for all signal in one cluster
#' @param data, cytof data
#' @param kluster, result from clustering which is a vector of 1 to number of clusters for each sample in dataset  
#' @param probs, the quantile used, for median 0.5
#' @return the signal value for this quantile for each cluster


q_per_cluster_marker<- function(data, cluster, probs){
  if(probs > 1){
    probs <- probs/100
  }
  result <- as.data.frame(matrix(NA, ncol = ncol(data), nrow = length(unique(cluster))))
  colnames(result) <- colnames(data)
  rownames(result) <- sort(unique(cluster))
  for(i in rownames(result)){
    result[i,] <- apply(data[as.character(cluster) %in% i,], 2, quantile, probs = probs)
  }
  return(result)
}



#' median_per_cluster_marker, calucate the 0.5 percentage for all signal in one cluster
#' @param data, cytof data
#' @param kluster, result from clustering which is a vector of 1 to number of clusters for each sample in dataset  
#' @return the signal value for the median for each cluster

median_per_cluster_marker <- function(data, cluster){
  q_per_cluster_marker(data = data, cluster = cluster, probs = 0.5)
}





#' col50, give 50 different  colors 

col50 <- c(
  "1" = "deepskyblue2",
  "2" = "#E31A1C", # red
  "3" = "green4",
  "4" = "#6A3D9A", # purple
  "5" = "#FF7F00", # orange
  "6" = "black", 
  "7" = "pink",
  "8" = "cyan4",
  "9" = "darkgoldenrod4",
  "10" = "blueviolet", 
  "11" = "darkorange4", 
  "12" = "coral", 
  "13" = "khaki2",
  "14" = "green3", 
  "15" = "orchid1", 
  "16" = "deeppink1", 
  "17" = "blue1", 
  "18" = "steelblue4",
  "19" = "darkturquoise", 
  "20" = "green1", 
  "21" = "yellow4", 
  "22" = "yellow3",
  "23" = "#FDBF6F", # lt orange
  "24" = "brown", 
  "25" = "gold1",
  "26" = "cyan2",
  "27" = "darkMagenta",
  "28" = "darkorange",
  "29" = "deeppink3", 
  "30" = "#CAB2D6", # lt purple
  "31" = "cornflowerblue",
  "32" = "darkslategrey", 
  "33" = "darkseagreen1", 
  "34" = "darkolivegreen4", 
  "35" = "azure3", 
  "36" = "brown", 
  "37" = "chocolate", 
  "38" = "darkgoldenrod2", 
  "39" = "chartreuse4", 
  "40" = "cadetblue4", 
  "41" = "darkorchid1", 
  "42" = "deeppink4",
  "43" = "palegreen2",
  "44" = "cornsilk4",
  "45" = "#FB9A99", # lt pink
  "46" = "darkolivegreen1", 
  "47" = "deepskyblue4",
  "48" = "dodgerblue2", 
  "49" = "gray70",
  "50" = "coral3"
  
  
)




pch25 <- c(
  "1" = 0, 
  "2" = 1,
  "3" = 2,
  "4" = 6,
  "5" = 0,
  "6" = 1,
  "7" = 2,
  "8" = 6,
  "9" = 0,
  "10" = 1,
  "11" = 2,
  "12" = 6,
  "13" = 0,
  "14" = 1,
  "15" = 2,
  "16" = 6,
  "17" = 0,
  "18" = 1,
  "19" = 2,
  "20" = 6,
  "21" = 0,
  "22" = 1,
  "23" = 2,
  "24" = 6,
  "25" = 0,
  "26" = 1
)


#' colfunc, give 4 different colors for time-signal plots. 
colfunc <- colorRampPalette(c("black", "purple4", "red", "yellow"))



#' number_of_events, calculate the number of cells/events in each sample 
#' @param data, list observations in all fcs files 
#' @param file_names, default NA will give the name 1,2,3 etc.
#' @return vector with number of events in each sample

number_of_events <- function(data, file_names = NA){
  events <- NULL
  number_of_files <- length(data)
  if(is.na(file_names[1])){
    file_names <- 1:number_of_files
  }
  file_names <- as.character(file_names)
  
  
  if(number_of_files > 1){
    for (i in 1:number_of_files){
      events[i] <- nrow(data[[i]])
    }
  }
  
  if(number_of_files > 1){
    names(events) <- file_names
  }
  return(events)
}


#' random_events, give list of random events for each sample (if numb_events given) or a vector of random events to use (if vector_with_filenames given)
#' @param numb_events, vector with number of events in each sample, use either vector_with_filenames or numb_events 
#' @param vector_with_filenames, vector of filename for each observation, use either either vector_with_filenames or numb_events 
#' @param n, number of events from wanted for each samples. Default equal 10000. 
#' If n greater than number of observation in a sample then n for that sample will be equal to number of observations 
#' @return list of vectors with position for the random events for each sub dataset
#' 

random_events <- function(numb_events = NULL, vector_with_filenames = NULL, n = 10000, posNeg = NULL, selected_events = NULL, values = 12){
  rand_events <- NULL
  if(!is.null(numb_events)){
    number_of_files <- length(numb_events)
    if(is.null(selected_events)){
      for (i in 1:number_of_files){
       rand_events[i] <- list(sort(sample(1:numb_events[i], min(numb_events[i], n))))
      } 
    } else {
      if(selected_events == "all"){
        for (i in 1:number_of_files){
          rand_events[i] <- list(sort(sample(1:numb_events[i], min(numb_events[i], n))))
        } 
      } else {
        if(is.null(posNeg)){
          print("a posNeg list with matrixes of the selected events has to be included")
        } else {
          rand_events <- random_events_from_selected_events(posNeg = posNeg, marker = selected_events, n = n, values = values)
        }
      }
    }
  } else {
    if(!is.null(vector_with_filenames)){
      datasets <- unique(vector_with_filenames)
        for(datasets_i in datasets){
          events <- which(vector_with_filenames %in% datasets_i)
          rand_events <- c(rand_events, sort(sample(events, n)))
        }
    }
  }
  return(rand_events)
}





#hm, har pr?vd ? sl? sammen med random_events..
#' #' random_events_vector give vector of random events for the whole dataset
#' #' @param datasetvector, vector of filename for each observation
#' #' @param n, number of events from wanted for each samples. Default equal 10000.
#' #' @return list of vectors with position for the random events for each sub dataset
#' #' 
#' 
#' random_events_vector <- function(datasetvector, n = 10000){
#'   datasets <- unique(datasetvector)
#'   rand_events <- NULL
#'   for(datasets_i in datasets){
#'     events <- which(datasetvector %in% datasets_i)
#'     rand_events <- c(rand_events, sort(sample(events, n)))
#'   }
#'   return(rand_events)
#' }
#' 
#' 
# HM Nytt navn????#' random_events_from_selected_events give list of random events for each sample
#' @param numb_events, vector with number of events in each sample
#' @param n, number of events from wanted for each samples. Default equal 10000.
#' If n greater than number of observation in a file then n for that file will be equal to number of observations
#' @return list of vectors with position for the random events for each sub dataset
#'

random_events_from_selected_events <- function(posNeg, marker, n = 10000, values = 12){
  number_of_files <- length(posNeg)
  rand_events <- NULL
  if(marker %in% colnames(posNeg[[1]])){
  for (i in 1:number_of_files){
    if(values == 12){
      possible <- which(posNeg[[i]][,marker] %in% c(1,2))
    } else {
      if(values == 1){
        possible <- which(posNeg[[i]][,marker] %in% c(1))
      } else {
        if(values == 2){
          possible <- which(posNeg[[i]][,marker] %in% c(2))
        } else {
          if(values == 0){
            possible <- which(posNeg[[i]][,marker] %in% c(0))
          }
          
        }
        
      }
      
    }
    rand_events[i] <- list(sort(sample(possible, min(length(possible), n))))
  }
  return(rand_events)
  } else {
    print(paste0(marker, " is not a column_name in posNeg[[1]]"))
  }
  
}



#' time_signal_plot, plot x = time and y = signal of marker in random_events for each sample
#' @param data, transformed data 
#' @param random_events, list of which events to plot for each sample
#' @param marker, which marker to plot
#' @param plot_title, vector with title for each plot, typical file names
#' @param prop_after_this_gating, proportion of events remaining after this gating
#' @param prop_final_event, propotion of events remaining in total
#' @param lower_gate, vector with values for lower gate or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate or NA (no upper gating)
#' @param time_div, value to divide time with, to get better values on x-axis, default 60*1000 which gives time in min_utes.
#' @param ylim
#' @return a list of time signal plots, one for each sample. 

time_signal_plot <- function(data, random_events, marker,  plot_title, 
                             prop_after_this_gating = NA, prop_final_event = NA, 
                             lower_gate = NA, upper_gate = NA, time_div = 60 * 1000, ylim = NA){
  marker <- ggplot2::sym(marker)
  plot_list <- list()
  for (i in 1:length(plot_title)){
    max_time <- max(data[[i]][random_events[[i]],"Time"]/time_div)
    gg <- ggplot2::ggplot(data[[i]][random_events[[i]],], ggplot2::aes(x=Time/time_div, y=!!marker)) +
      #scale on x axis 
      ggplot2::scale_x_continuous(breaks=seq(0,round(max_time ,1),round(max_time /2,1))) + 
      # Plot all points
      ggplot2::geom_point(shape=".",alpha=0.5) + 
      # Fill with transparent colour fill using density stats
      # ndensity scales each graph to its own min/max values
      ggplot2::stat_density2d(geom="raster", ggplot2::aes(fill=..ndensity.., alpha = ..ndensity..), 
                              contour = FALSE) +
      # Produces a colour scale based on the colours in the colfunc list
      ggplot2::scale_fill_gradientn(colours=colfunc(128)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Time (min)")+
      ggplot2::coord_cartesian(expand=FALSE)
    # Draw gate
    if(!is.na(lower_gate[1]) ){
      gate_line <- data.frame(y0 = lower_gate[i], xmax = max_time )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = 0, xend = xmax, y = y0, yend = y0), 
                                       color = "black") 
    }
    if(!is.na(upper_gate[1])){
      gate_line <- data.frame(y1 = upper_gate[i], xmax = max_time )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = 0, xend = xmax, y = y1, yend = y1), 
                                       color = "black") 
    }
    #title
    if(is.na(prop_after_this_gating) & is.na(prop_final_event)){
      title <- plot_title[i]
    } else {
      if(!is.na(prop_after_this_gating) & is.na(prop_final_event)){    
        title <- paste0(round(prop_after_this_gating[i]*100,1)," %,     ", plot_title[i])
      } else{
        if(is.na(prop_after_this_gating) & !is.na(prop_final_event)){
          title <- paste0(round(prop_final_event[i]*100,1), " % of total,    ", plot_title[i])
        } else {
          title <- paste0(round(prop_after_this_gating[i]*100,1)," % (", round(prop_final_event[i]*100,1), " % of total),    ", plot_title[i])
        }
      }
    }
    
    gg <- gg + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(size=8))
    
    if(!is.na(ylim)[1]){
      gg <- gg + ggplot2::coord_cartesian(ylim = ylim, default = TRUE) 
    }
    
    
    plot_list[[i]] <- gg
  }
  return(plot_list)
}




#######################


#' density_plot, plot density for all events in each sample
#' @param data, transformed data 
#' @param marker, which marker to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param lower_gate, vector with values for lower gate or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate or NA (no upper gating)
#' @param xlim, xlim default NA.
#' @return density plots

density_plot <- function(data, marker, plot_title = NA, lower_gate = NA, upper_gate = NA, xlim = NA, main_title = "", included_files = NA, max_event_used = NA, minimum_signal = NA){
  #browser()
  number_of_files <- length(data)  
  possible_i <- 1:number_of_files 
  
  if(!is.na(included_files[1])){
    possible_i <- included_files
  }
  
  if(is.na(plot_title[1])){
    plot_title <- as.character(possible_i)
  }
  n_used <- max_event_used
  plot_title_nr <- 1:length(possible_i)
  column <- which(colnames(data[[possible_i[1]]]) == marker)
  if(is.na(max_event_used)){
    df <- data.frame(Values = data[[possible_i[1]]][,column], Sample = rep(plot_title[possible_i[1]], nrow(data[[possible_i[1]]])))
    for(i in possible_i[2:length(possible_i)]){
      df <- rbind(df, data.frame(Values = data[[i]][,column], Sample = rep(plot_title[i], nrow(data[[i]]))))
    }
  } else {
    n_used <- min(c(max_event_used, nrow(data[[possible_i[1]]])))
    tamed <- sample(nrow(data[[possible_i[1]]]), n_used)
    df <- data.frame(Values = data[[possible_i[1]]][tamed,column], Sample = rep(plot_title[possible_i[1]], length(tamed)))
    for(i in possible_i[2:length(possible_i)]){
      n_used <- min(c(max_event_used, nrow(data[[possible_i[i]]])))
      tamed <- sample(nrow(data[[possible_i[i]]]), n_used)
      df <- rbind(df, data.frame(Values = data[[i]][tamed,column], Sample = rep(plot_title[i], length(tamed))))
    }
  }
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = Values, y = Sample, col = Sample, fill = Sample, group = Sample))+
    ggjoy::geom_joy(scale = 2, alpha=0.5, show.legend= F) +
    ggplot2::scale_y_discrete(expand=c(0.01, 0)) +
    ggplot2::scale_x_continuous(expand=c(0.01, 0)) +
    ggplot2::ggtitle(marker) +
    ggplot2::ggtitle(main_title)+
    ggjoy::theme_joy() 
  
  if(!is.na(lower_gate[1]) ){
    gate_line <- data.frame(Sample = plot_title, x0 = lower_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x0, xend = x0, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "black") 
    # gg <- gg + ggplot2::geom_vline(data = gate_line, xintercept = x0, col = Sample)
  }
  if(!is.na(upper_gate[1])){
    gate_line <- data.frame(Sample = plot_title,  x1 = upper_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x1, xend = x1, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "black") 
    #  gg <- gg + ggplot2::geom_vline(data = gate_line, ggplot2::aes(xintercept = x1, color = Sample))
  }
  
  if(!is.na(xlim)[1]){
    gg <- gg + ggplot2::coord_cartesian(xlim = xlim) 
  }
  
  return(gg)
}






#' density_plot_without_neg, plot density for positive events in each sample
#' @param data, transformed data 
#' @param marker, which marker to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param lower_gate, vector with values for lower gate or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate or NA (no upper gating)
#' @param xlim, xlim default NA.
#' @return density plots

density_plot_without_neg <- function(data, marker, plot_title = NA, lower_gate = NA, upper_gate = NA, xlim = NA, main_title = "",included_files = NA){
  
  if(length(lower_gate) == 1){
    lower_gate <- rep(lower_gate, length(data))
  }
  
  
  
  number_of_files <- length(data)
  
  
  if(is.na(included_files[1])){
    possible_i <- 1:number_of_files 
  } else {
    possible_i <- included_files
  }
  
  if(is.na(plot_title[1])){
    plot_title <- as.character(possible_i)
  }
  
  
  plot_title_nr <- 1:possible_i
  column <- which(colnames(data[[1]]) == marker)
  xx <- data[[possible_i[1]]][,column]
  xx <- xx[xx > lower_gate[possible_i[1]]]
  df <- data.frame(Values = xx, Sample = rep(plot_title[possible_i[1]], length(xx)))
  for(i in possible_i[2:length(possible_i)]){
    xx <- data[[i]][,column]
    xx <- xx[xx > lower_gate[i]]
    df <- rbind(df, data.frame(Values = xx, Sample = rep(plot_title[i], length(xx))))
  }
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = Values, y = Sample, col = Sample, fill = Sample, group = Sample))+
    ggjoy::geom_joy(scale = 2, alpha=0.5, show.legend= F) +
    ggplot2::scale_y_discrete(expand=c(0.01, 0)) +
    ggplot2::ggtitle(main_title)+
    ggplot2::scale_x_continuous(expand=c(0.01, 0)) +
    ggplot2::ggtitle(marker) +
    ggjoy::theme_joy()
  
  if(!is.na(lower_gate[1]) ){
    gate_line <- data.frame(Sample = plot_title, x0 = lower_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x0, xend = x0, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "blue", size = 1.3) 
    # gg <- gg + ggplot2::geom_vline(data = gate_line, xintercept = x0, col = Sample)
  }
  if(!is.na(upper_gate[1])){
    gate_line <- data.frame(Sample = plot_title,  x1 = upper_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x1, xend = x1, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "blue", size = 1.3) 
    #  gg <- gg + ggplot2::geom_vline(data = gate_line, ggplot2::aes(xintercept = x1, color = Sample))
  }
  
  if(!is.na(xlim)[1]){
    gg <- gg + ggplot2::coord_cartesian(xlim = xlim) 
  }
  return(gg)
}



density_plot_selected_cells <- function(data, marker, include, mark, positiv = TRUE, plot_title = NA, lower_gate = NA, upper_gate = NA, xlim = NA, main_title = "", included_files = NA){
  
  if(main_title == ""){
    main_title <- marker
  }
  number_of_files <- length(data)
  if(is.na(included_files[1])){
    possible_i <- 1:number_of_files 
  } else {
    possible_i <- included_files
  }
  
  if(is.na(plot_title[1])){
    plot_title <- as.character(possible_i)
  }
  
  
  plot_title_nr <- 1:possible_i
  
  
  column <- which(colnames(data[[1]]) == marker)
  df <- data.frame(Values = data[[possible_i[1]]][,column], Sample = rep(plot_title[1], nrow(data[[possible_i[1]]])))
  for(i in possible_i[2:length(possible_i)]){
    if(positiv){
      xx <- data[[i]][include[[i]][,mark],column]
    } else {
      xx <- data[[i]][!include[[i]][,mark],column]
    }
    # xx <- xx[xx > lower_gate]
    df <- rbind(df, data.frame(Values = xx, Sample = rep(plot_title[i], length(xx))))
  }
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = Values, y = Sample, col = Sample, fill = Sample))+
    ggjoy::geom_joy(scale = 2, alpha=0.5, show.legend= F) +
    ggplot2::scale_y_discrete(expand=c(0.01, 0)) +
    ggplot2::ggtitle(main_title) +
    ggplot2::scale_x_continuous(expand=c(0.01, 0)) +
    ggjoy::theme_joy()
  
  if(!is.na(lower_gate[1]) ){
    gate_line <- data.frame(Sample = plot_title, x0 = lower_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x0, xend = x0, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "black") 
    # gg <- gg + ggplot2::geom_vline(data = gate_line, xintercept = x0, col = Sample)
  }
  if(!is.na(upper_gate[1])){
    gate_line <- data.frame(Sample = plot_title,  x1 = upper_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x1, xend = x1, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "black") 
    #  gg <- gg + ggplot2::geom_vline(data = gate_line, ggplot2::aes(xintercept = x1, color = Sample))
  }
  
  if(!is.na(xlim)[1]){
    gg <- gg + ggplot2::coord_cartesian(xlim = xlim) 
  }
  
  return(gg)
}






#' signal_signal_just_plot, scatterplot of two different signals
#' @param data, transformed data 
#' @param marker1, which marker to plot
#' @param marker2, which marker to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param xlim, xlim default NA.
#' @param ylim, ylim default NA.
#' @return scatterplots of two differnt signals per file. 


signal_signal_just_plot <- function(data, random_events, marker1, marker2, 
                                    plot_title = NA, xlim = NA, ylim = NA){
  marker1 <- ggplot2::sym(marker1)
  marker2 <- ggplot2::sym(marker2)
  columnVar1 <- which(colnames(data[[1]]) == marker1)
  columnVar2 <- which(colnames(data[[1]]) == marker2)
  
  
  
  if(is.na(plot_title[1])){
    plot_title <- paste0("file ", 1:length(data))
  }
  
  plotList <- list()
  
  for (i in 1:length(data)){
    var1_i <- data[[i]][, columnVar1]
    if(is.na(xlim[1])){
      xlim <- c(min(var1_i), max(var1_i))
    }
    
    var2_i <- data[[i]][, columnVar2]
    if(is.na(ylim[1])){
      ylim <- c(min(var2_i), max(var2_i))
    }
    
    gg <- ggplot2::ggplot(data[[i]][random_events[[i]],], ggplot2::aes(x=!!marker1, y=!!marker2)) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +      
      # Plot all points
      ggplot2::geom_point(shape=".",alpha=0.5)+
      # Fill with transparent colour fill using density stats
      # ndensity scales each graph to its own min/max values
      ggplot2::stat_density2d(geom="raster", ggplot2::aes(fill=..ndensity.., alpha = ..ndensity..), contour = FALSE) +
      # Produces a colour scale based on the colours in the colfunc list
      ggplot2::scale_fill_gradientn(colours=colfunc(128)) +
      ggplot2::theme(legend.position = "none")  
    #   stat_ellipse(level = 0.8)
    plotList[[i]] <- gg
  }
  return(list(plotList = plotList))
  
  
}



#' signal_signal_plot, scatterplot of two different signals
#' @param data, transformed data 
#' @param marker1, which marker to plot
#' @param marker2, which marker to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param xlim, xlim default NA.
#' @param ylim, ylim default NA.
#' @return scatterplots of two differnt signals per file. 


signal_signal_plot <- function(data, random_events, marker1, marker2, xname = marker1, yname = marker2, 
                               xlow = NA, ylow = NA, xhigh = NA, yhigh = NA, 
                               plot_title = NA, xlim = NA, ylim = NA, title_size = 10, contour = FALSE){
  marker1 <- ggplot2::sym(marker1)
  marker2 <- ggplot2::sym(marker2)
  columnVar1 <- which(colnames(data[[1]]) == marker1)
  columnVar2 <- which(colnames(data[[1]]) == marker2)
  
  
  
  if(is.na(plot_title[1])){
    plot_title <- paste0("file ", 1:length(data))
  }
  
  plotList <- list()
  
  for (i in 1:length(data)){
    var1_i <- data[[i]][random_events[[i]], columnVar1]
    if(is.na(xlim[1])){
      xlim <- c(min(var1_i), max(var1_i))
    }
    
    var2_i <- data[[i]][random_events[[i]], columnVar2]
    if(is.na(ylim[1])){
      ylim <- c(min(var2_i), max(var2_i))
    }
    #browser()
    
    gg <- ggplot2::ggplot(data[[i]][random_events[[i]],], ggplot2::aes(x=!!marker1, y=!!marker2)) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +      
      # Plot all points
      ggplot2::geom_point(shape=".",alpha=0.5)
    # Fill with transparent colour fill using density stats
    # ndensity scales each graph to its own min/max values
    if(contour == F){  
      gg <- gg +  ggplot2::stat_density2d(geom="raster", ggplot2::aes(fill=..ndensity.., alpha = ..ndensity..), contour = FALSE) + 
        ggplot2::scale_fill_gradientn(colours=colfunc(128))
      # Produces a colour scale based on the colours in the colfunc list
    } else  {
      gg <- gg + ggplot2::stat_density2d(ggplot2::aes(fill = stat(level)), geom = "polygon", colour="white")
    }
    gg <- gg  +
      ggplot2::theme(legend.position = "none")  +
      ggplot2::xlab(xname) +
      ggplot2::ylab(yname) +
      ggplot2::ggtitle(plot_title[i]) + 
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size, face = "bold"))
    
    if(!is.na(xlow[1]) ){
      gate_line <- data.frame(x0 = xlow[i], ymax = ylim[2], ymin = ylim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = x0, xend = x0, y = ymin, yend = ymax), 
                                       color = "blue", size = 1.3) 
    }
    if(!is.na(xhigh[1]) ){
      gate_line <- data.frame(x0 = xhigh[i], ymax = ylim[2], ymin = ylim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = x0, xend = x0, y = ymin, yend = ymax), 
                                       color = "blue", size = 1.3) 
    }    
    if(!is.na(ylow[1]) ){
      gate_line <- data.frame(y0 = ylow[i], xmax = xlim[2], xmin = xlim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = xmin, xend = xmax, y = y0, yend = y0), 
                                       color = "blue", size = 1.3) 
    }
    if(!is.na(yhigh[1]) ){
      gate_line <- data.frame(y0 = yhigh[i], xmax = xlim[2], xmin = xlim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = xmin, xend = xmax, y = y0, yend = y0), 
                                       color = "blue", size = 1.3) 
    }    
    #   stat_ellipse(level = 0.8)
    plotList[[i]] <- gg
  }
  return(list(plotList = plotList))
  
  
}





#' signal_signal_plot_selected_cells, scatterplot of two different signals
#' @param data, transformed data 
#' @param marker1, which marker to plot
#' @param marker2, which marker to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param xlim, xlim default NA.
#' @param ylim, ylim default NA.
#' @param include, list with matrices over which cells to include'
#' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' @param mark, which column from include to use
#' @return scatterplots of two differnt signals per file. 


signal_signal_plot_selected_cells <- function(data, number_random_events = 10000,  marker1, marker2, xname = marker1, yname = marker2, 
                                              xlow = NA, ylow = NA, xhigh = NA, yhigh = NA, 
                                              plot_title = NA, xlim = NA, ylim = NA, title_size = 10, include, mark, positiv = TRUE){
  marker1 <- ggplot2::sym(marker1)
  marker2 <- ggplot2::sym(marker2)
  columnVar1 <- which(colnames(data[[1]]) == marker1)
  columnVar2 <- which(colnames(data[[1]]) == marker2)
  
  
  
  if(is.na(plot_title[1])){
    plot_title <- paste0("file ", 1:length(data))
  }
  
  plotList <- list()
  
  for (i in 1:length(data)){
    if(positiv){
      var1_i <- data[[i]][include[[i]][,mark], columnVar1]
      var2_i <- data[[i]][include[[i]][,mark], columnVar2]
      data_i <- data[[i]][include[[i]][,mark],]
    } else {
      var1_i <- data[[i]][!include[[i]][,mark], columnVar1]
      var2_i <- data[[i]][!include[[i]][,mark], columnVar2]
      data_i <- data[[i]][!include[[i]][,mark],]
    }
    
    if(nrow(data_i) > number_random_events){
      data_i <- data_i[sample(1:nrow(data_i), number_random_events),]
    }
    
    if(is.na(xlim[1])){
      xlim <- c(min(var1_i), max(var1_i))
    }
    
    var2_i <- data[[i]][include[[i]][,mark], columnVar2]
    if(is.na(ylim[1])){
      ylim <- c(min(var2_i), max(var2_i))
    }
    
    gg <- ggplot2::ggplot(data_i, ggplot2::aes(x=!!marker1, y=!!marker2)) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +      
      # Plot all points
      ggplot2::geom_point(shape=".",alpha=0.5)
    # Fill with transparent colour fill using density stats
    # ndensity scales each graph to its own min/max values
    if(contour == F){  
      gg <- gg +  ggplot2::stat_density2d(geom="raster", ggplot2::aes(fill=..ndensity.., alpha = ..ndensity..), contour = FALSE) 
      # Produces a colour scale based on the colours in the colfunc list
    } else  {
      gg <- gg + ggplot2::stat_density2d(aes(fill = stat(level)), geom = "polygon")
    }
    
    gg <- gg +  ggplot2::scale_fill_gradientn(colours=colfunc(128)) +
      ggplot2::theme(legend.position = "none")  +
      ggplot2::xlab(xname) +
      ggplot2::ylab(yname) +
      ggplot2::ggtitle(plot_title[i]) + 
      ggplot2::theme(plot.title = ggplot2::element_text(size = title_size, face = "bold"))
    
    if(!is.na(xlow[1]) ){
      gate_line <- data.frame(x0 = xlow[i], ymax = ylim[2], ymin = ylim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = x0, xend = x0, y = ymin, yend = ymax), 
                                       color = "blue", size = 1.3) 
    }
    if(!is.na(xhigh[1]) ){
      gate_line <- data.frame(x0 = xhigh[i], ymax = ylim[2], ymin = ylim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = x0, xend = x0, y = ymin, yend = ymax), 
                                       color = "blue", size = 1.3) 
    }    
    if(!is.na(ylow[1]) ){
      gate_line <- data.frame(y0 = ylow[i], xmax = xlim[2], xmin = xlim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = xmin, xend = xmax, y = y0, yend = y0), 
                                       color = "blue", size = 1.3) 
    }
    if(!is.na(yhigh[1]) ){
      gate_line <- data.frame(y0 = yhigh[i], xmax = xlim[2], xmin = xlim[1] )
      gg <- gg + ggplot2::geom_segment(data = gate_line, 
                                       ggplot2::aes(x = xmin, xend = xmax, y = y0, yend = y0), 
                                       color = "blue", size = 1.3) 
    }    
    #   stat_ellipse(level = 0.8)
    plotList[[i]] <- gg
  }
  return(list(plotList = plotList))
  
  
}





#' density_plot_per_cluster, plot density for all markers in each cluster
#' @param data, transformed data 
#' @param cluster_per_cell, vector of clusters
#' @param rand_events, which events to plot
#' @param plot_cluster, default NA give alle clusters, else vector of clusters to plot
#' @param nrow_plot, tells how many rows of plots to export, default = 4,
#' @param strip_text_size, give size of stripe text, defaut = 10, 
#' @param legend_text_size, give size of legend text, defaut = 10, 
#' @param axis_text_size, give size of axis text, defaut = 8
#' @return density plots for clusters against rest of data. 



density_plot_per_cluster <- function(data, cluster_per_cell, rand_events = NA, plot_cluster = NA, nrow_plot = 4, strip_text_size = 10, legend_text_size = 10, axis_text_size = 8){
  if(!is.na(rand_events[1])){
    data <- data[rand_events,]
    cluster_per_cell <- cluster_per_cell[rand_events]
  }  
  
  if(is.na(plot_cluster[1])){
    unique_cluster_per_cell <- sort(unique(cluster_per_cell))
  } else {
    unique_cluster_per_cell <- plot_cluster
  }
  
  n_cluster_per_cell <- length(unique_cluster_per_cell)
  
  for(i in 1:n_cluster_per_cell){
    d2 <- as.data.frame(cbind(data, cluster_per_cell %in% unique_cluster_per_cell[i]))
    colnames(d2) <- c(colnames(data),  "cluster")
    d2$id <- 1:nrow(d2)
    d3 <- reshape2::melt(d2, id.vars = c("id", "cluster"))
    d3$cluster[d3$cluster == TRUE] <- paste("cluster", unique_cluster_per_cell[i])
    d3$cluster[d3$cluster == FALSE] <- "the rest"
    g <- ggplot2::ggplot(d3, ggplot2::aes(x = value, color = cluster)) + 
      ggplot2::geom_density(size = 2) +
      
      ggplot2::facet_wrap(~ variable, nrow = nrow_plot, scales = "free") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), 
                     strip.text = ggplot2::element_text(size = strip_text_size), 
                     axis.text = ggplot2::element_text(size = axis_text_size),
                     legend.text=ggplot2::element_text(size = legend_text_size),
                     legend.title = ggplot2::element_blank()) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = 1)) #+
    print(paste("cluster", unique_cluster_per_cell[i]))
    print(g)
  }
  
}


#' barplot_per_sample, barplot of events in each cluster, per sample
#' @param file_names_per_cell, transformed data 
#' @param cluster_per_cell, vector of clusters
#' @param rand_events, which events to plot
#' @return density plots for clusters against rest of data. 

barplot_per_sample <- function(file_names_per_cell, cluster_per_cell, rand_events = NA, mar = c(5.1, 4.1, 4.1, 8.1)){
  if(!is.na(rand_events[1])){
    file_names_per_cell <- file_names_per_cell[rand_events]
    cluster_per_cell <- cluster_per_cell[rand_events]
  }  
  
  n_clusters <- length(unique(cluster_per_cell))
  tab <- table( cluster_per_cell, file_names_per_cell) 
  
  par(mar=mar, xpd=TRUE)
  barplot(tab, col = col25, las = 2)
  legend("topright", inset=c(-0.1,0), col= c(col40, col40)[n_clusters:1], legend = names(c(col50, col50)[n_clusters:1]), pch =15)
}


violin_per_sample <- function(data_mat, x, y, colour = NA, med = NA, q10 = NA, q90 = NA, main = ""){
  # browser()
  xcol <- which(colnames(data_mat) == x)
  ycol <- which(colnames(data_mat) == y)
  colcol <- which(colnames(data_mat) == colour)
  d2 <- data_mat[, c(xcol, ycol, colcol)]
  colnames(d2) <- c("x", "y", "Status")
  g <-  ggplot(d2, aes( x = x, y = y, fill = Status)) + 
    geom_violin() +
    labs(title = main)
  if(!is.na(med)){
    g <- g + geom_hline(aes(yintercept = med), lwd = 1.5)  
  }
  if(!is.na(q10)){
    g <- g + geom_hline(aes(yintercept = q10), lty = 5, lwd = 1.5)  
  }
  if(!is.na(q90)){
    g <- g + geom_hline(aes(yintercept = q90), lty = 5, lwd = 1.5)  
  }
  
  
  return(g)
}


#' Read fcs files in folder
#' 
#' @param data_path a path to the folder with fcs files
#' @param files_to_open list of files to open, if NULL all files will be used
#' @return fcs_data a set of all data produced. 

read_data_from_folder <- function(data_path, files_to_open = NULL){
  if(is.null(files_to_open)){
    fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
    fcs_files <- fcs_files[grepl(".fcs", fcs_files)]
    files_to_open <- basename(fcs_files)
  } else {
     if(files_to_open == "all"){
        fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
        fcs_files <- fcs_files[grepl(".fcs", fcs_files)]
        files_to_open <- basename(fcs_files)
     } else {
       files_to_open <- files_to_open[grepl(".fcs", files_to_open)]
       setwd(dirname(fcs_files[1]))
       file_names <- gsub(".fcs", "", files_to_open)
     }
  }
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


#' transform_selected_markers transform data with arc_sinh or log, might include more
#' @params fcs_data,  a set of fcs files
#' @params markers, transforms only those markers listed
#' @params method, possible values are "arc_sinh" and "log", "arc_sinh" is default 
#' @params cofactor can be changed from 5 that is default
#' @params scale, if scaling wanted this have to be true, default false
#' @params scaling, interval to scale each parameter to, only used if scale = true
#' @return scaled data for those markers chosen

transform_selected_markers <- function(fcs_data, markers, method = "arc_sinh", cofactor = 5, scale = F,  scaling = c(-5, 12000)){
  new_data <- NULL 
  number_of_files <- length(fcs_data)
  
  for (i in 1:number_of_files){
    if(method == "arc_sinh"){
      new_data[[i]] <- as.data.frame(asinh(flowCore::exprs(fcs_data[[i]][,markers])/cofactor))
    } else {
      if(method == "log"){
        new_data[[i]] <- as.data.frame(log(flowCore::exprs(fcs_data[[i]][,markers]) + 1))
      } else {
        print("method must be given as either arc_sinh or log")
      }
    }
    if(scale){
      for (j in 1:length(markers)){
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


#' max_marker find max value for the marker
#' @params fcs_data,  a set of fcs files
#' @params markers, transforms only those markers listed
max_marker <- function(data, marker){
  marker_max <- max(data[[1]][,marker])
  for(i in 1:n_files){
    marker_max <- max(marker_max, max(data[[i]][,marker]))
  }
  return(marker_max)
}



#' min_marker find min value for the marker
#' @params fcs_data,  a set of fcs files
#' @params markers, transforms only those markers listed
min_marker <- function(data, marker){
  marker_min <- min(data[[1]][,marker])
  for(i in 1:n_files){
    marker_min <- min(marker_min, min(data[[i]][,marker]))
  }
  return(marker_min)
}



##her



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
#' @param markers, kanaler som skal være med
#' @param arcsinh, default true
#' @param cofactor, default 5
#' @param scale, default F
#' @param scaling, default c(-5,12000), only used when scale = T
#' @param cytofData, default True
#' @return new dataset with only those events that we want to keep. 

list_to_matrix_selected_events <- function(data, kept_events, file_names = NA, markers, arcsinh = T, 
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
    mat <- as.data.frame(data[[1]][kept_events[[1]], markers])
    mat$dataset <- file_names[1]
    if(n > 1){
      for(i in 2:n){
        mat0 <- as.data.frame(data[[i]][kept_events[[i]], markers])
        mat0$dataset <- file_names[i]
        mat <- rbind(mat, mat0)
      }
    }
  } else {
    
    
    if(arcsinh == T){
      if(scale == F){
        mat <- as.data.frame(asinh(flowCore::exprs(data[[1]][kept_events[[1]], markers])/cofactor))
        mat$dataset <- file_names[1]
        mat$Time <- flowCore::exprs(fcs_data[[1]][kept_events[[1]],"Time"])
        if(n > 1){
          for(i in 2:n){
            mat0 <- as.data.frame(asinh(flowCore::exprs(data[[i]][kept_events[[i]], markers])/cofactor))
            mat0$dataset <- file_names[i]
            mat0$Time <- flowCore::exprs(fcs_data[[i]][kept_events[[i]],"Time"])
            mat <- rbind(mat, mat0)
          }
        }
      } else {  #rescale(new_data[[i]][,j],to = scaling)
        mat <- as.data.frame(rescale(asinh(flowCore::exprs(data[[1]][kept_events[[1]], markers])/cofactor),to = scaling))
        mat$dataset <- file_names[1]
        mat$Time <- flowCore::exprs(fcs_data[[1]][,"Time"])
        if(n > 1){
          for(i in 2:n){
            mat0 <- as.data.frame(rescale(asinh(flowCore::exprs(data[[i]][kept_events[[i]], markers])/cofactor),to = scaling))
            mat0$dataset <- file_names[i]
            mat0$Time <- flowCore::exprs(fcs_data[[i]][,"Time"])
            mat <- rbind(mat, mat0)
          }
        }
      }
      
    } else {
      mat <- as.data.frame(flowCore::exprs(data[[1]][kept_events[[1]], markers]))
      mat$dataset <- file_names[1]
      mat$Time <- flowCore::exprs(fcs_data[[1]][,"Time"])
      if(n > 1){
        for(i in 2:n){
          mat0 <- as.data.frame(flowCore::exprs(data[[i]][kept_events[[i]], markers]))
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


#' 
#' 
#' #' find_gate_lower_noise_plus_high_low_selected_cells find position for gating for each subdataset
#' #' @param data fcs_data sets
#' #' @param marker which marker used for gating
#' #' @return value of gating for the given marker in each subdataset. 
#' 
#' find_gate_lower_noise_plus_high_low_selected_cells <- function(data, marker, include, mark, positiv = TRUE){
#'   column <- which(colnames(data[[1]]) == marker)
#'   results <- rep(NA, length(data))
#'   for(i in 1:length(data)){
#'     if(positiv){
#'       res <-  find_gate_lower_noise_per_file_plus_high_low(xx = data[[i]][include[[i]][,mark], column])
#'     } else {
#'       res <-  find_gate_lower_noise_per_file_plus_high_low(xx = data[[i]][!include[[i]][,mark], column])
#'     }
#'     lower_gates[i] <- res[[1]]
#'     upper_gates[i] <- res[[2]]
#'   }
#'   return(list(lower_gates = lower_gates, upper_gates = upper_gates))  
#' }
#' 
#' 
#' 
#' 
#' #' find_gate_lower_noise_plus_high_low find position for gating for each subdataset
#' #' @param data fcs_data sets
#' #' @param marker which marker used for gating
#' #' @return value of gating for the given marker in each subdataset. 
#' 
#' find_gate_lower_noise_plus_high_low <- function(data, marker){
#'   column <- which(colnames(data[[1]]) == marker)
#'   results <- rep(NA, length(data))
#'   for(i in 1:length(data)){
#'     results[i] <- find_gate_lower_noise_per_file_plus_high_low(data[[i]][, column])
#'   }
#'   return(results)
#' }
#' 
#' 
#' #' find_gate_lower_noise_per_file_plus_high_low function
#' #' @param xx vector of values
#' #' @return the value that correspond to the percentage in the density plot.
#' 
#' find_gate_lower_noise_per_file_plus_high_low <- function(xx){
#'   dens <- density(xx)
#'   #  ts_y<-ts(smooth(dens$y))
#'   #  tp <- pastecs::turnpoints(ts_y)
#'   #  top1 <- dens$x[tp$peaks][1]
#'   top1 <- dens$x[dens$y == max(dens$y)]
#'   firstGate <- top1 + (top1 - dens$x[1])
#'   
#'   xxx <- xx[xx > firstGate]
#'   dens <- density(xxx)
#'   browser()
#'   
#'   return(value)
#' }
#' 
#' 
#' 

#' find_gate_lower_noise find position for gating for each subdataset
#' @param data fcs_data sets
#' @param marker which marker used for gating
#' @return value of gating for the given marker in each subdataset.

find_gate_lower_noise <- function(data, marker){
  column <- which(colnames(data[[1]]) == marker)
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



#' 
#' #' find_gate_gaussian_first_top_per_file function
#' #' @param xx vector of values
#' #' @param upper_perc how many percentage to trow away.
#' #' @return the value that correspond to the percentage in the density plot.
#' find_gate_gaussian_first_top_per_file <- function(xx, perc_included = 0.9995){
#'   fit <- density(xx)
#'   x.new <- rnorm(10000, sample(xx, size = 10000, replace = TRUE), fit$bw)
#'   #plot(fit)
#'   #lines(density(x.new), col = "red")
#'   fit <- VGAM::vglm(x.new ~ 1, 
#'                     VGAM::mix2normal(eq.sd = F), iphi=0.5, imu= 0, isd1=2, imu2=5, 
#'                     isd2=1)
#'   #  fit2 <- vglm(x.new ~ 1, uninormal(), lmean = 0, lsd = 1)
#'   
#'   # Calculated parameters
#'   pars <- as.vector(coef(fit))
#'   w <- VGAM::logitlink(pars[1], inverse=TRUE, )
#'   m1 <- pars[2]
#'   sd1 <- exp(pars[3])
#'   m2 <- pars[4]
#'   sd2 <- exp(pars[5])
#'   gate <- m1 + qnorm(perc_included) * sd1
#'   return(gate)
#' }
#' 
#' 
#' #' find_gate_gaussian_first_top find position for gating for each subdataset
#' #' @param data fcs_data sets
#' #' @param marker which marker used for gating
#' #' @param perc_included how many percentages that are assumed to be noise
#' #' @return value of gating for the given marker in each subdataset. 
#' 
#' find_gate_gaussian_first_top <- function(data, marker, perc_included = 0.9995){
#'   column <- which(colnames(data[[1]]) == marker)
#'   results <- rep(NA, length(data))
#'   for(i in 1:length(data)){
#'     results[i] <- find_gate_gaussian_first_top_per_file(data[[i]][, column], perc_included = perc_included)
#'   }
#'   return(results)
#' }
#' 
#' 
#' 
#' 
#' events_to_keep, find which events to keep in each subdataset based on lower and/or upper gate.
#' @param data, transformed data
#' @param marker, which marker to plot
#' @param lower_gate, vector with values for lower gate, a number (same lower gate for all subset) or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate, a number (same upper gate for all subset)  or NA (no upper gating)
#' @return list of vectors of true/false for each events (true means keep), one vector for each subdataset.


events_to_keep <- function(data, marker, lower_gate = NA, upper_gate = NA){
  number_of_files <- length(data)
  column <- which(colnames(data[[1]]) == marker)
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

#' 
#' #' percent_to_keep_this_gating
#' #' @param data, data 
#' #' @param kept_events, list of vectors of true/false for alle events in each subdataset. Those events that are true will be kept
#' #' @param file_names, default NA will give the name 1,2,3 etc.
#' #' @return percentage kept this gating
#' #' 
#' percent_to_keep_this_gating <- function(kept_events, file_names = NA){
#'   number_of_files <- length(kept_events)
#'   if(is.na(file_names[1])){
#'     file_names <- 1:number_of_files
#'   }
#'   file_names <- as.character(file_names)
#'   percent <- rep(NA, number_of_files)
#'   for (i in 1:number_of_files){
#'     percent[i] <- table(kept_events[[i]])["TRUE"]/length(kept_events[[i]])
#'   }
#'   names(percent) <- file_names
#'   return(percent)
#' }
#' 
#' 
#' #' update_data_based_on_events_to_keep
#' #' @param data, data 
#' #' @param kept_events, list of vectors of true/false for alle events in each subdataset. Those events that are true will be kept
#' #' @return new dataset with only those events that we want to keep. 
#' 
#' update_data_based_on_events_to_keep <- function(data, kept_events){
#'   number_of_files <- length(data)
#'   for (i in 1:number_of_files){
#'     data[[i]] <- data[[i]][kept_events[[i]],]
#'   }
#'   return(data)
#' }
#' 
#' 
#' 
#' 
#' #' find_gate_second_top, find the gaussian gates of the second top for a vector xx
#' #' @param xx, vector of numbers 
#' #' @param lower_gate_prop_height, propotions for lower gate
#' #' @param upper_gate_prop_height, propotions for upper gate
#' #' @return list of lower and upper gates 
#' 
#' # find_gate_second_top <- function(xx, lower_gate_prop_height, upper_gate_prop_height, perc_included){
#' #   dens <- density(xx)
#' #   ts_y<-ts(smooth(dens$y))
#' #   tp <- pastecs::turnpoints(ts_y)
#' #   bunn1 <- dens$x[tp$pits][1]
#' #   xx[xx < bunn1] <- NA
#' #   dens <- density(xx[!is.na(xx)])  
#' #   lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop_height])
#' #   upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop_height])
#' #   return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
#' # }
#' 
#' 
#' 
#' 
#' # 
#' # 
#' # find_gate_second_top <- function(xx, lower_gate_prop_height, upper_gate_prop_height, perc_included, main_top_to_left = F){
#' # #  browser()
#' #   dens <- density(xx)
#' #   ts_y<-ts(smooth(dens$y))
#' #   tp <- pastecs::turnpoints(ts_y)
#' #   bunn1 <- dens$x[tp$pits][1]
#' #   xx[xx < bunn1] <- NA
#' #   dens <- density(xx[!is.na(xx)])
#' #   if(!is.na(lower_gate_prop_height[1])){
#' #     lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop_height])
#' #     upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop_height])
#' #   } else{
#' # 
#' #       xx <- xx[!is.na(xx)]
#' #       fit <- density(xx)
#' #       x.new <- rnorm(10000, sample(xx, size = 10000, replace = TRUE), fit$bw)
#' #        #plot(fit)
#' #          #lines(density(x.new), col = "red")
#' #          fit <- VGAM::vglm(x.new ~ 1,
#' #                           VGAM::mix2normal(eq.sd = F), iphi=0.5, imu= 0, isd1=2, imu2=5,
#' #                           isd2=1)
#' #        #  fit2 <- vglm(x.new ~ 1, uninormal(), lmean = 0, lsd = 1)
#' # 
#' #            # Calculated parameters
#' #            pars <- as.vector(coef(fit))
#' #            w <- VGAM::logitlink(pars[1], inverse=TRUE)
#' #            m1 <- pars[2]
#' #            sd1 <- exp(pars[3])
#' #            m2 <- pars[4]
#' #            sd2 <- exp(pars[5])
#' #            if(main_top_to_left == TRUE){
#' #              lower_gate <- m1 - qnorm(perc_included) * sd1
#' #              upper_gate <- m1 + qnorm(perc_included) * sd1
#' #            } else {
#' #              lower_gate <- m2 - qnorm(perc_included) * sd2
#' #              upper_gate <- m2 + qnorm(perc_included) * sd2
#' # 
#' #            }
#' # 
#' #   }
#' #   return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
#' # }
#' # 
#' # 
#' # 
#' # 
#' # 
#' # 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 






#' 
#' 
#' 
#' #' find_gaussian_gates_second_top_selected_cells, find the split between the first and secound top for all subsets
#' #' @param data, data
#' #' @param data, data
#' #' @param marker, which marker to plot
#' #' @param include, list with matrices over which cells to include'
#' #' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' #' @param mark, which column from include to use
#' #' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset)
#' #' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)
#' #' @return list of vectors with lower and upper gates for each subset.
#' 
#' find_gaussian_gates_second_top_top_selected_cells <- function(data, marker, lower_gate_percent = NA, upper_gate_percent = NA, perc_included, main_top_to_left ,  minimum = 0, include, mark, positiv = TRUE){
#'   # browser()
#'   column <- which(colnames(data[[1]]) == marker)
#'   if(!is.na(lower_gate_percent[[1]])){
#'     if(lower_gate_percent >= 1){
#'       lower_gate_prop_height <- lower_gate_percent/100
#'     } else {
#'       lower_gate_prop_height <- lower_gate_percent
#'     }
#'     if(upper_gate_percent >= 1){
#'       upper_gate_prop_height <- upper_gate_percent/100
#'     } else {
#'       upper_gate_prop_height <- upper_gate_percent
#'     }
#'   } else {
#'     lower_gate_prop_height <- NA
#'     upper_gate_prop_height <- NA
#'   }
#'   lower_gates <- rep(NA, length(data))
#'   upper_gates <- rep(NA, length(data))
#'   for(i in 1:length(data)){
#'     if(positiv){
#'       res <-  find_gate_second_top(xx = data[[i]][include[[i]][,mark], column], lower_gate_prop_height = lower_gate_prop_height, upper_gate_prop_height = upper_gate_prop_height, perc_included = perc_included, main_top_to_left = main_top_to_left, minimum = minimum)
#'     } else {
#'       res <-  find_gate_second_top(xx = data[[i]][!include[[i]][,mark], column], lower_gate_prop_height = lower_gate_prop_height, upper_gate_prop_height = upper_gate_prop_height, perc_included = perc_included, main_top_to_left = main_top_to_left, minimum = minimum)
#'     }
#'     lower_gates[i] <- res[[1]]
#'     upper_gates[i] <- res[[2]]
#'   }
#'   return(list(lower_gate = lower_gates, upper_gate = upper_gates))
#' }
#' 


#' 
#' 
#' 
#' #' find_split, find the first bottom
#' #' @param xx, vector of numbers 
#' #' @return vector of splits
#' 
#' find_split <- function(xx, minimum){
#'   dens <- density(xx)
#'   ts_y <- ts(smooth(dens$y))
#'   tp <- pastecs::turnpoints(ts_y)
#'   bunn1 <- min(dens$x[tp$pits][dens$x[tp$pits] > minimum])
#'   return(bunn1)
#' }
#' 
#' #' find_split_first_second_top, find the split between the first and secound top for all subsets
#' #' @param data, data 
#' #' @param marker, which marker to plot
#' #' @return vector of splits
#' 
#' find_split_first_second_top <- function(data, marker, minimum = 0){
#'   column <- which(colnames(data[[1]]) == marker)
#'   splits <- rep(NA, length(data))
#'   for(i in 1:length(data)){
#'     splits[[i]] <-  find_split(data[[i]][, column], minimum)
#'   }
#'   return(splits)
#' }
#' 
#' 
#' 
#' 
#' 
#' 
#' #' find_split_first_second_top_selected_cells, find the split between the first and secound top for all subsets
#' #' @param data, data 
#' #' @param marker, which marker to plot
#' #' @param include, list with matrices over which cells to include'
#' #' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' #' @param mark, which column from include to use
#' #' @return vector of splits
#' 
#' find_split_first_second_top_selected_cells <- function(data, marker, minimum = 0, include, mark, positiv = TRUE){
#'   column <- which(colnames(data[[1]]) == marker)
#'   splits <- rep(NA, length(data))
#'   for(i in 1:length(data)){
#'     if(positiv){
#'       splits[[i]] <-  find_split(data[[i]][include[[i]][,mark], column], minimum)
#'     } else {
#'       splits[[i]] <-  find_split(data[[i]][!include[[i]][,mark], column], minimum)
#'     }
#'   }
#'   return(splits)
#' }
#' 
#' 
#' 
#' #' find_split_first_second_top, find the split between the first and secound top for all subsets
#' #' @param data, data 
#' #' @param marker, which marker to plot
#' #' @return vector of splits
#' 
#' find_split_neg_low_high <- function(data, marker, neg = 0.05, minLowHigh = 0.1){
#'   column <- which(colnames(data[[1]]) == marker)
#'   splits <- rep(NA, length(data))
#'   neg <- rep(NA, length(data))
#'   lower_gates <- splits
#'   upper_gates <- splits
#'   for(i in 1:length(data)){
#'     xx <- data[[i]][,column]
#'     neg[[i]] <- find_split(xx, minimum = 0)
#'     #     neg[[i]] <- find_gate_first_top(xx, lower_gate_prop_height = negProp, upper_gate_prop_height = negProp)$upper_gate
#'     splits[[i]] <-  find_split(xx[xx > neg[[i]]], minimum = minLowHigh)
#'     
#'   }
#'   
#'   return(list(neg_splits = neg, low_high_splits = splits))
#'   
#'   
#'   
#' }
#' 
#' 
#' 
#' 
#' #' find_split_first_second_top, find the split between the first and secound top for all subsets
#' #' @param data, data 
#' #' @param marker, which marker to plot
#' #' @param include, list with matrices over which cells to include'
#' #' @param positiv, default equal TRUE, include those cells that are positiv for mark in include list
#' #' @param mark, which column from include to use
#' 
#' #' @return vector of splits
#' 
#' find_split_neg_low_high_selected_cells <- function(data, marker, include, mark, positiv = TRUE, neg = 0.05, minLowHigh = 0.1){
#'   column <- which(colnames(data[[1]]) == marker)
#'   splits <- rep(NA, length(data))
#'   neg <- rep(NA, length(data))
#'   lower_gates <- splits
#'   upper_gates <- splits
#'   for(i in 1:length(data)){
#'     if(positiv){
#'       xx <- data[[i]][include[[i]][,mark],column]
#'     } else{
#'       xx <- data[[i]][!include[[i]][,mark],column]
#'     }
#'     neg[[i]] <- find_split(xx, minimum = 0)
#'     splits[[i]] <-  find_split(xx[xx > neg[[i]]], minimum = minLowHigh)
#'   }
#'   return(list(neg_splits = neg, low_high_splits = splits))
#' }









#' find_gaussian_gates find the gaussian gates 
#' @param data, data
#' @param marker, which marker to plot
#' @param method, "first top", "second top" or "highest top", "upper noise"
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset)
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)
#' @param minimum, minimum values for where to look for top, except highest top which is minimum of upper gate.
#' @param perc_included,
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates <- function(data, marker, method, lower_gate_percent = NA, upper_gate_percent = NA, minimum = NA,  perc_included = NA){
  if(method == "first top"){
    if(is.na(perc_included)){
      find_gaussian_gates_first_top(data = data, marker = marker, lower_gate_percent = lower_gate_percent, upper_gate_percent = upper_gate_percent, min_upper_gate = minimum)
    } else {
      find_gaussian_gates_first_top(data = data, marker = marker, lower_gate_percent = NA, upper_gate_percent = NA, perc_included = perc_included, main_top_to_left = T, min_upper_gate = minimum)
    }
  } else {
    if(method == "second top"){
      if(is.na(perc_included)){
        find_gaussian_gates_second_top(data = data, marker = marker, lower_gate_percent = lower_gate_percent, upper_gate_percent = upper_gate_percent, minimum = minimum)
      } else {
        find_gaussian_gates_second_top(data = data, marker = marker, lower_gate_percent = NA, upper_gate_percent = NA, perc_included = perc_included, main_top_to_left = F, minimum = minimum)
      }
    } else {
      if(method == "highest top"){
        find_gaussian_gates_highest_top(data = data, marker = marker, lower_gate_percent = lower_gate_percent, upper_gate_percent = upper_gate_percent, min_upper_gate = minimum)
      } else {
        if(method == "upper noise"){
          if(is.na(upper_gate_percent)){
            upper_gate_percent <- lower_gate_percent
          }
          find_gate_perc_height_upper_noise(data = data, marker =  marker, upper_gate_prop_height = upper_gate_percent/100)
        } else {
          if(method == "lower noise"){
            find_gate_lower_noise(data = data, marker = marker)
          } else {
            print("You have to use the argument top; top could be first top, second top, highest top, upper noise or lower noise")
          }
        }
      }
    }
  }
}
  


#' find_gate_first_top, find the gaussian gates of the first top for a vector xx
#' @param xx, vector of numbers
#' @param lower_gate_prop_height, propotions for lower gate
#' @param upper_gate_prop_height, propotions for upper gate
#' @return list of lower and upper gates

find_gate_first_top <- function(xx, lower_gate_prop_height, upper_gate_prop_height, min_upper_gate = NA){
  dens <- density(xx)
  ts_y <- ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  bunn1 <- dens$x[tp$pits][1]
  xx[xx > bunn1] <- NA
  dens <- density(xx[!is.na(xx)])
  lower_gate <- max(min(dens$x[dens$y > max(dens$y) * lower_gate_prop_height]))
  upper_gate <- max(c(dens$x[dens$y > max(dens$y) * upper_gate_prop_height], min_upper_gate), na.rm = T)
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}




#' find_gaussian_gates_first_top, find the gaussian gates of the first top for all subsets
#' @param data, data
#' @param marker, which marker to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset)
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_first_top <- function(data, marker, lower_gate_percent, upper_gate_percent, min_upper_gate = NA){
  column <- which(colnames(data[[1]]) == marker)
  if(lower_gate_percent > 1){
    lower_gate_prop_height <- lower_gate_percent/100
  } else {
    lower_gate_prop_height <- lower_gate_percent
  }
  if(upper_gate_percent > 1){
    upper_gate_prop_height <- upper_gate_percent/100
  } else {
    upper_gate_prop_height <- upper_gate_percent
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
  for(i in 1:length(data)){
    res <-  find_gate_first_top(data[[i]][, column], lower_gate_prop_height = lower_gate_prop_height, upper_gate_prop_height = upper_gate_prop_height, min_upper_gate = min_upper_gate)
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates ))
}




#' find_gaussian_gates_highest_top, find the gaussian gates of the highest top for all subsets
#' @param data, data
#' @param marker, which marker to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset)
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)
#' @param min_upper_gate, minimum value for upper gate
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_highest_top <- function(data, marker, lower_gate_percent, upper_gate_percent, min_upper_gate = NA){
  column <- which(colnames(data[[1]]) == marker)
  if(lower_gate_percent > 1){
    lower_gate_prop_height <- lower_gate_percent/100
  } else {
    lower_gate_prop_height <- lower_gate_percent
  }
  if(upper_gate_percent > 1){
    upper_gate_prop_height <- upper_gate_percent/100
  } else {
    upper_gate_prop_height <- upper_gate_percent
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
  for(i in 1:length(data)){
    res <-  find_gate_highest_top(data[[i]][, column], lower_gate_prop_height = lower_gate_prop_height, upper_gate_prop_height = upper_gate_prop_height, min_upper_gate = min_upper_gate)
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates ))
}




#' find_gaussian_gates_second_top, find the gaussian gates of the second top for all subsets
#' @param data, data
#' @param marker, which marker to plot
#' @param lower_gate_percent, vector with percentage for lower gate, a number (same percentage for all subset)
#' @param upper_gate_percent,  vector with percentage for uppe gate, a number (same percentage for all subset)
#' @return list of vectors with lower and upper gates for each subset.

find_gaussian_gates_second_top <- function(data, marker, lower_gate_percent = NA, upper_gate_percent = NA, perc_included = NA, main_top_to_left = F , minimum = NA){
  
  column <- which(colnames(data[[1]]) == marker)
  if(!is.na(lower_gate_percent[[1]])){
    if(lower_gate_percent >= 1){
      lower_gate_prop_height <- lower_gate_percent/100
    } else {
      lower_gate_prop_height <- lower_gate_percent
    }
    if(upper_gate_percent >= 1){
      upper_gate_prop_height <- upper_gate_percent/100
    } else {
      upper_gate_prop_height <- upper_gate_percent
    }
  } else {
    lower_gate_prop_height <- NA
    upper_gate_prop_height <- NA
  }
  lower_gates <- rep(NA, length(data))
  upper_gates <- rep(NA, length(data))
  for(i in 1:length(data)){
    res <-  find_gate_second_top(xx = data[[i]][, column], lower_gate_prop_height = lower_gate_prop_height, upper_gate_prop_height = upper_gate_prop_height, perc_included = perc_included, main_top_to_left = main_top_to_left, minimum = minimum)
    
    
    lower_gates[i] <- res[[1]]
    upper_gates[i] <- res[[2]]
  }
  return(list(lower_gates = lower_gates, upper_gates = upper_gates))
}


#' find_gate_highest_top, find the gaussian gates of the main top for a vector xx
#' @param xx, vector of numbers
#' @param lower_gate_prop_height, propotions for lower gate
#' @param upper_gate_prop_height, propotions for upper gate
#' @return list of lower and upper gates

find_gate_highest_top <- function(xx, lower_gate_prop_height, upper_gate_prop_height, min_upper_gate = NA){
 #browser()
  dens <- density(xx)
  ts_y <- ts(smooth(dens$y))
  top_position <- which(dens$y == max(dens$y))[1]
  possible <- which(dens$y < max(dens$y) * lower_gate_prop_height) 
  lower_gate <- dens$x[max(possible[possible < top_position])]
  possible <- which(dens$y < max(dens$y) * upper_gate_prop_height) 
  upper_gate <- max(c(dens$x[min(possible[possible > top_position])], min_upper_gate), na.rm = T)
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}




#' find_gate_second_top, find the gaussian gates of the secund top for a vector xx
#' @param xx, vector of numbers
#' @param lower_gate_prop_height, propotions for lower gate; NA if perc_included used
#' @param upper_gate_prop_height, propotions for upper gate; NA if perc_included used
#' @param perc_included, percentage of gaussian top included
#' @param main_top_to_left = F, only needed when perc_included used, second top is equal to main_top_to_left
#' @param minimum = NA
#' @return list of lower and upper gates

find_gate_second_top <- function(xx, lower_gate_prop_height, upper_gate_prop_height, perc_included, main_top_to_left = F, minimum = NA){
  dens <- density(xx)
  ts_y<-ts(smooth(dens$y))
  tp <- pastecs::turnpoints(ts_y)
  if(!is.na(minimum)){
    mini <- which(dens$x > minimum)[1]
  } else {
    mini <- 0
  }
  bunn_n <- min(max(which(tp$pits)[1], mini),  length(dens$y) - 1, na.rm = T)
  top_n <-  bunn_n + min(which(dens$y[bunn_n:length(dens$y)] == max(dens$y[bunn_n:length(dens$y)]))[1],  length(bunn_n:length(dens$y)), na.rm = T) - 1
  h <- dens$y[top_n] - dens$y[bunn_n]
  if(!is.na(lower_gate_prop_height[1])){
    cutoff_h_lower <- dens$y[bunn_n] + lower_gate_prop_height*h
    cutoff_n_lower <- bunn_n + which(dens$y[bunn_n:length(dens$y)] > cutoff_h_lower)[1] - 1
    lower_gate <- min(dens$x[cutoff_n_lower])
    upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop_height])
  } else{

    temp <- vgam_gates(xx,  perc_included = perc_included, main_top_to_left = main_top_to_left)
    lower_gate <- temp$lower_gate
    upper_gate <- temp$upper_gate
    
  }
  return(list(lower_gate = lower_gate, upper_gate = upper_gate ))
}


vgam_gates <- function(xx, perc_included, main_top_to_left){
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
    temp <- find_gate_highest_top(xx, lower_gate_prop_height = 0.06, upper_gate_prop_height = 0.06)
    lower_gate <- temp$lower_gate
    upper_gate <- temp$upper_gate
  }
  
  return(list(lower_gate = lower_gate, upper_gate = upper_gate))
}


gaus_gates <- function(dens,  lower_gate_prop_height, upper_gate_prop_height){
  lower_gate <- min(dens$x[dens$y > max(dens$y) * lower_gate_prop_height])
  upper_gate <- max(dens$x[dens$y > max(dens$y) * upper_gate_prop_height])
  return(list(lower_gate = lower_gate, upper_gate = upper_gate))
  
}






#' find_gate_perc_upper_noise find position for gating for each subdataset
#' @param data fcs_data sets
#' @param marker which marker used for gating
#' @param upper_gate_prop_height how many percentages that are assumed to be noise
#' @return value of gating for the given marker in each subdataset.

find_gate_perc_height_upper_noise <- function(data, marker, upper_gate_prop_height = 0.001){
  column <- which(colnames(data[[1]]) == marker)
  results <- rep(NA, length(data))
  for(i in 1:length(data)){
    results[i] <- find_gate_upper_noise_per_file(data[[i]][, column], upper_gate_prop_height = upper_gate_prop_height)
  }
  return(results)
}


#' find_gate_upper_noise_per_file function
#' @param xx vector of values
#' @param upper_gate_prop_height how many percentage to trow away.
#' @return the value that correspond to the percentage in the density plot.

find_gate_upper_noise_per_file <- function(xx, upper_gate_prop_height){
  dens <- density(xx)
  posible <- dens$y > max(dens$y) * upper_gate_prop_height
  xposible <- which(posible[2:length(posible)] - posible[1:(length(posible)-1)] == -1)[1]
  value <- dens$x[xposible]
  return(value)
}




#' marker_plot
#' @param path path to where the result files from flowSOM analysis is
#' @param k which k to make plot for
#' @param seed which seed to make plot for
#' @param highlight_cluster, number of cluster or vector of clusters that should get an other color than grey
#' @param gates, optional can be included for non, all or some markers
#' @param order_marker_shortname, if a specific order of the figure per marker is wanted.
#' @return value of gating for the given marker in each subdataset.
marker_plot <- function(path, k, seed, highlight_cluster = NULL, selectedEvents = "all", gates = NULL, order_marker_shortname = NULL){
  
  q5 <- read.csv2(fs::path(path, paste0("q5_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  q10 <- read.csv2(fs::path(path, paste0("q10_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  q25 <- read.csv2(fs::path(path, paste0("q25_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  q50 <- read.csv2(fs::path(path, paste0("medians_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  q75 <- read.csv2(fs::path(path, paste0("q75_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  q90 <- read.csv2(fs::path(path, paste0("q90_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  q95 <- read.csv2(fs::path(path, paste0("q95_per_cluster_k_", k, "_seed", seed, selectedEvents,".csv")), stringsAsFactors = FALSE) 
  
  colnames(q5)[1] <- "Cluster"
  colnames(q10)[1] <- "Cluster"
  colnames(q25)[1] <- "Cluster"
  colnames(q50)[1] <- "Cluster"
  colnames(q75)[1] <- "Cluster"
  colnames(q90)[1] <- "Cluster"
  colnames(q95)[1] <- "Cluster"

  names <- read.csv2(fs::path(paths$clean_data_info_path, "marker_names_included_manual_shortnames.csv"), stringsAsFactors = FALSE)

  q5long <- data.table::melt(data.table::setDT(q5),id.vars = "Cluster", variable.names = "marker")
  colnames(q5long)[3] <- "q5"
  q10long <- data.table::melt(data.table::setDT(q10),id.vars = "Cluster", variable.names = "marker")
  colnames(q10long)[3] <- "q10"
  q25long <- data.table::melt(data.table::setDT(q25),id.vars = "Cluster", variable.names = "marker")
  colnames(q25long)[3] <- "q25"
  q50long <- data.table::melt(data.table::setDT(q50),id.vars = "Cluster", variable.names = "marker")
  colnames(q50long)[3] <- "q50"
  q75long <- data.table::melt(data.table::setDT(q75),id.vars = "Cluster", variable.names = "marker")
  colnames(q75long)[3] <- "q75"
  q90long <- data.table::melt(data.table::setDT(q90),id.vars = "Cluster", variable.names = "marker")
  colnames(q90long)[3] <- "q90"
  q95long <- data.table::melt(data.table::setDT(q95),id.vars = "Cluster", variable.names = "marker")
  colnames(q95long)[3] <- "q95"
  
  d <-merge(q5long, q10long)
  d <-merge(d, q25long)
  d <-merge(d, q50long)
  d <-merge(d, q75long)
  d <- merge(d, q90long)
  d <- merge(d, q95long)
  d <- merge(d, names, by.x= "variable", by.y = "marker_name")
  if(!is.null(gates)){
    d <- merge(d, gates, by.x = "marker_short_name", by.y = "X")
  }
  d$X <- NULL
  
  if(!is.null(order_marker_shortname[1])){
    if(!is.na(order_marker_shortname[1])){
      if(all(d$marker_short_name %in%  order_marker_shortname )){
        d$marker_short_name <- factor(d$marker_short_name, levels = order_marker_shortname)
      } else {
        d$marker_short_name <- factor(d$marker_short_name)
        print("could not change order of markers in plot due to some markers not included in list (or misspelled)")
      }
    } else {
      d$marker_short_name <- factor(d$marker_short_name)
    }
  } else {
    d$marker_short_name <- factor(d$marker_short_name)
  }

  d$col <- "0"
  if(!is.null(highlight_cluster[1])){
    if(!is.na(highlight_cluster[1])){
        for(i in highlight_cluster){
          d$col[d$Cluster == i] <- i
        }
    }
  }
  

  unique_col <- unique(d$col)
  
  our_colors <- c("grey", c(col50, col50))
  names(our_colors) <- 0:100
  used_colors <- our_colors[as.character(unique_col)]
  
  g <- ggplot(d, aes(x = q50, y = Cluster, xmin = q10, xmax = q90, col = col)) + 
    geom_point(size = 2) + 
    geom_errorbar() +
    geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q5, xmax = q95, col = col))+ 
    geom_errorbar(data = d, aes(x = q50, y = Cluster, xmin = q25, xmax = q75, col = col))+ 
    facet_wrap(vars(marker_short_name), ncol = 11, scale = "free_x") + 
    scale_color_manual(values = used_colors)
 
  if(!is.null(gates)){ 
    g <- g + geom_vline(aes(xintercept  = low), linetype = "dashed", size = 1.1) + 
    geom_vline(aes(xintercept  = high), linetype = "dashed", size = 1.1) 
  }
  return(g)
}



#' get_cluster_from_quantiles
#' @param lower_upper_limites, matrix with the columns: marker_short_name, lower, upper, either lower upper or both should be specified for all markers.
#' @param data, dataset list of matrixes that are transformed in the same manner as lower and upper.
#' @param selected_markers, either "all" or vector of marker_short_names to beused
#' @param correlation, 
#' @return list of rest_prop and
get_cluster_from_quantiles <- function(lower_upper_limites, data, selected_markers  = "all"){
 # browser()
  if(selected_markers [1] == "all"){
    markers <- lower_upper_limites$marker_short_name 
  } else {
    markers <- selected_markers 
  }
  
  lower_upper_limites <- lower_upper_limites[lower_upper_limites$marker_short_name %in% markers,]
  temp <- lower_upper_limites[,c("lower", "upper")]
  both_NA <- apply(is.na(temp), 1, sum) == 2
  lower_upper_limites <-lower_upper_limites[!both_NA,]

  res_prop <- rep(NA, length(data))
  result <- list()
  for(j in 1:length(data)){
    mat <-  as.data.frame(matrix(NA, ncol = nrow(lower_upper_limites), nrow =  nrow(data[[j]])))
    colnames(mat) <- lower_upper_limites$marker_short_name
    result[[j]] <- mat
  }
  
  for(marker_i in as.character(lower_upper_limites$marker_short_name)){
    marker_name_i <- as.character(lower_upper_limites[lower_upper_limites$marker_short_name == marker_i, "marker_name"]) 
    lower_i <- lower_upper_limites[lower_upper_limites$marker_short_name == marker_i, "lower"]
    upper_i <- lower_upper_limites[lower_upper_limites$marker_short_name == marker_i, "upper"]
    for(j in 1:length(data)){
      if(!is.na(lower_i) & !is.na(upper_i)){
        result[[j]][,marker_i] <- data[[j]][, marker_name_i] >= lower_i & data[[j]][, marker_name_i] <= upper_i
      } else {
        if(is.na(upper_i)){
          result[[j]][,marker_i] <- data[[j]][, marker_name_i] >= lower_i 
        } else {
          result[[j]][,marker_i] <- data[[j]][, marker_name_i] <= upper_i
        }
      }
    }
  }
  
 cor_per_markers <- NULL
  
  for(j in 1:length(data)){
    temp <- apply(result[[j]], 1, sum) == nrow(lower_upper_limites)
    if(!is.na( table(temp)["TRUE"])){
      res_prop[j] <- table(temp)["TRUE"]/length(temp)
      cor_per_markers <- rbind(cor_per_markers, cor(temp,result[[j]]))
    } else {
      res_prop[j] <- 0
    }
    median_cor_per_markers <- apply(cor_per_markers, 2, median)
  }
    
 return(list(res_prop = res_prop, median_cor_per_markers = median_cor_per_markers))
  
}




#' get_marker_name_from_marker_short_name
#' @param marker_info tabel with both marker_name and marker_short name, can be found in ..\clean_up\clean_data_info
#' @param marker_short_name shorter marker name
#' @return makes lots of files in the result folder
get_marker_name_from_marker_short_name <- function(marker_info, marker_short_name){
  marker_name <- rep(NA, length(marker_short_name))
  for(i in 1:length(marker_short_name)){
    marker_name[i] <- marker_info$marker_name[marker_info$marker_short_name == marker_short_name[i]]
  }
  return(marker_name)
}



#' run_flowSOM
#' @param path path to where the result files from flowSOM analysis is
#' @param k which k to make plot for
#' @param seed which seed to make plot for
#' @param included_files = file_names, # this will include all files in the clean data folder, can be changed to a vector of filenames to use.  
#' @param n_per_file number of events to include per file
#' @param included_markers, use marker_name, might use function get_marker_name_from_marker_short_name(marker_info, marker_short_name) 
#' @param transformation = "arc_sinh"
#' @param scaling_flowSOM = TRUE
#' @param k_s number of clusters, k_s has to be a number or a vector of numbers.
#' @param xdim xdim * ydim gives the number of nodes that FlowSOM will combined to k_s clusters.
#' @param ydim
#' @param resultpath paths$clean_data_flowSOM_results_path, can be changed if new folder is wanted for different analysis
#' @param seed the seed ensure that the same events are chosen, and hence the same result are uptained next time the exactly same analysis are done
#' @param selectedEvents
#' @param posNeg
#' @param make_heatmap = TRUE
#' @param heatmap_cluster_column = FALSE
#' @return makes lots of files in the result folder

run_flowSOM <- function(fcs_data, file_names, included_files = "all", n_per_file, included_markers = "all", transformation = "arc_sinh", scaling_flowSOM = TRUE, k_s, xdim = 10, ydim = 10, resultpath, seed, selectedEvents = "all", posNeg = NULL, make_heatmap = TRUE, heatmap_cluster_column = FALSE){
  if(!file.exists(resultpath)){
    dir.create(file.path(resultpath))
  }
  set.seed(seed)
  if(is.null(selectedEvents)){
    number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
    random_events_data <- random_events(number_of_events_data, n = n_per_file)
    ext_name <- "all"
   # print("NULL")
  } else {
     if(selectedEvents == "all"){
        number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
        random_events_data <- random_events(number_of_events_data, n = n_per_file)
        ext_name <- "all"
      #  print("all")
     } else {
        random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = selectedEvents, n = n_per_file)
        ext_name <- selectedEvents
      #  print(selectedEvents)
     }
  }
  
  
  if(transformation != "arc_sinh"){
    print("the data analysis will  be done on arc sinh transformeddata.")
  }
  
  if(included_markers[1] == "all"){
    included_markers <- colnames(fcs_data[[1]])
  }
  

  
  arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data, 
                                                     kept_events = random_events_data, 
                                                     file_names = file_names, 
                                                     markers = included_markers,
                                                     arcsinh = T, 
                                                     cofactor = 5, 
                                                     scale = F
  )
  
  colnames(arcSindataMatrix)[1:length(included_markers)] <- included_markers
  
  
  if(included_files[1] == "all"){
    included_files <- file_names
  }
  
  data_to_analyse <- arcSindataMatrix[arcSindataMatrix$dataset %in% included_files, ]
  
  
  set.seed(seed) # to ensure same plot every time
  out <- FlowSOM::ReadInput(as.matrix(data_to_analyse[,included_markers]), transform = F, scale = scaling_flowSOM)
  out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(data_to_analyse[,included_markers])), xdim = xdim, ydim = ydim)
  out <- FlowSOM::BuildMST(out)
  cluster_FlowSOM_pre <- out$map$mapping[, 1]
  
  for(k in k_s){
    set.seed(seed)
    out_k <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)
    cluster_FlowSOM_k <- out_k[cluster_FlowSOM_pre]
    cluster_FlowSOM_k_factor <- factor(cluster_FlowSOM_k, levels = 1:k)
    q5_k <- q_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k, probs = 0.05)
    write.csv2(q5_k, fs::path(resultpath, paste0("q5_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    q10_k <- q_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k, probs = 0.1)
    write.csv2(q10_k, fs::path(resultpath, paste0("q10_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    q25_k <- q_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k, probs = 0.25)
    write.csv2(q25_k, fs::path(resultpath, paste0("q25_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    q75_k <- q_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k, probs = 0.75)
    write.csv2(q75_k, fs::path(resultpath, paste0("q75_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    q90_k <- q_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k, probs = 0.9)
    write.csv2(q90_k, fs::path(resultpath, paste0("q90_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    q95_k <- q_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k, probs = 0.95)
    write.csv2(q95_k, fs::path(resultpath, paste0("q95_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    medians_k <- median_per_cluster_marker(data = data_to_analyse[,included_markers], cluster = cluster_FlowSOM_k)
    write.csv2(medians_k, fs::path(resultpath, paste0("medians_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    
    if(make_heatmap == TRUE){
      tiff(fs::path(resultpath, paste0("heatmap_median_k_", k, "_cluster_seed", seed, ext_name, ".tiff")), width = 1000, height = 800)
      print(Heatmap(as.matrix(medians_k[,included_markers], cluster_columns = heatmap_cluster_column)))
      dev.off()
    }
    
    total_number_of_events_in_cluster <- table(cluster_FlowSOM_k_factor)
    
    number_of_events_per_cluster <- cbind(total_number_of_events_in_cluster, total_number_of_events_in_cluster/nrow(data_to_analyse))
    colnames(number_of_events_per_cluster) <- c("number_of_events", "percent")
    write.csv2(number_of_events_per_cluster, fs::path(resultpath, paste0("number_of_events_per_cluster_k_", k, "_seed", seed, ext_name, ".csv")))
    
    perSample <- table(data_to_analyse$dataset, cluster_FlowSOM_k_factor)
    colnames(perSample) <- paste0("cluster_", 1:k)
    
    write.csv2(perSample, fs::path(resultpath, paste0("number_of_events_per_cluster_and_file_k_", k, "_seed", seed, ext_name, ".csv")))
    
  }
}
