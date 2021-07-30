#' col25, give 25 different  colors 

col25 <- c(
  "1" = "dodgerblue2", 
  "2" = "#E31A1C", # red
  "3" = "green4",
  "4" = "#6A3D9A", # purple
  "5" = "#FF7F00", # orange
  "6" = "black", 
  "7" = "skyblue2", 
  "8" = "#FB9A99", # lt pink
  "9" = "palegreen2",
  "10" = "#CAB2D6", # lt purple
  "11" = "#FDBF6F", # lt orange
  "12" = "gray70", 
  "13" = "khaki2",
  "14" = "maroon", 
  "15" = "orchid1", 
  "16" = "deeppink1", 
  "17" = "blue1", 
  "18" = "steelblue4",
  "19" = "darkturquoise", 
  "20" = "green1", 
  "21" = "yellow4", 
  "22" = "yellow3",
  "23" = "darkorange4", 
  "24" = "brown", 
  "25" = "gold1"
)


#' colfunc, give 4 different colors for time-signal plots. 
colfunc <- colorRampPalette(c("black", "purple4", "red", "yellow"))



#' number of cells in each file
#' @param data, list observations in all fcs files
#' @return vector with number of cells in each subdataset

number_of_cells <- function(data){
  number_of_files <- length(data)
  cells <- NULL
  for (i in 1:number_of_files){
    cells[i] <- nrow(data[[i]])
  }
  return(cells)
}


#' random_cells give list of random cells for each subdataset
#' @param numb_cells, vector with number of cells in each subdataset
#' @param n, number of cells from wanted for each subdatasets. Default equal 10000.
#' @return list of vectors with position for the random cells for each sub dataset
#' 

random_cells <- function(numb_cells, n = 10000){
  number_of_files <- length(numb_cells)
  rand_cells <- NULL
  for (i in 1:number_of_files){
    rand_cells[i] <- list(sort(sample(1:numb_cells[i], min(numb_cells[i], n))))
  }
  return(rand_cells)
}

#' random_cells_vector give vector of random cells for the whole dataset
#' @param 
#' @param n, number of cells from wanted for each subdatasets. Default equal 10000.
#' @return list of vectors with position for the random cells for each sub dataset
#' 

random_cells_vector <- function(datasetvector, n = 10000){
  datasets <- unique(datasetvector)
  rand_cells <- NULL
  for(datasets_i in datasets){
    cells <- which(datasetvector %in% datasets_i)
    rand_cells <- c(rand_cells, sort(sample(cells, n)))
  }
  return(rand_cells)
}




#' time_signal_plot, plot x = time and y = signal of channel in random_cells for each subdataset
#' @param data, transformed data 
#' @param random_cells, list of which cells to plot for each subdataset
#' @param channel, which channel to plot
#' @param plot_title, vector with title for each plot, typical file names
#' @param prop_after_this_gating, proportion of cells remaining after this gating
#' @param prop_final_event, propotion of cells remaining in total
#' @param lower_gate, vector with values for lower gate or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate or NA (no upper gating)
#' @param time_div, value to divide time with, to get better values on x-axis, default 60*1000 which gives time in min_utes.
#' @param ylim
#' @return a list of time signal plots, one for each subdataset. 

time_signal_plot <- function(data, random_cells, channel,  plot_title, 
                             prop_after_this_gating = NA, prop_final_event = NA, 
                             lower_gate = NA, upper_gate = NA, time_div = 60 * 1000, ylim = NA){
  channel <- ggplot2::sym(channel)
  plot_list <- list()
  for (i in 1:length(plot_title)){
    max_time <- max(data[[i]][random_cells[[i]],"Time"]/time_div)
    gg <- ggplot2::ggplot(data[[i]][random_cells[[i]],], ggplot2::aes(x=Time/time_div, y=!!channel)) +
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
      gg <- gg + ggplot2::coord_cartesian(ylim = ylim) 
    }
    
    
    plot_list[[i]] <- gg
  }
  return(plot_list)
}



#' density_plot, plot density for all cells in each subdataset
#' @param data, transformed data 
#' @param channel, which channel to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param lower_gate, vector with values for lower gate or NA (no lower gating)
#' @param upper_gate,  vector with values for uppe gate or NA (no upper gating)
#' @param xlim, xlim default NA.
#' @return density plots

density_plot <- function(data, channel, plot_title = NA, lower_gate = NA, upper_gate = NA, xlim = NA){
  number_of_files <- length(data)
  if(is.na(plot_title[1])){
    plot_title <- as.character(1:number_of_files)
  }
  column <- which(colnames(data[[1]]) == channel)
  df <- data.frame(Values = data[[1]][,column], Sample = rep(plot_title[1], nrow(data[[1]])))
  for(i in 2:number_of_files){
    df <- rbind(df, data.frame(Values = data[[i]][,column], Sample = rep(plot_title[i], nrow(data[[i]]))))
  }
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = Values, y = Sample, col = Sample, fill = Sample))+
    ggjoy::geom_joy(scale = 2, alpha=0.5, show.legend= F) +
    ggplot2::scale_y_discrete(expand=c(0.01, 0)) +
    ggplot2::scale_x_continuous(expand=c(0.01, 0)) +
    ggplot2::ggtitle(channel) +
    ggjoy::theme_joy()
  
  if(!is.na(lower_gate[1]) ){
    gate_line <- data.frame(Sample = plot_title, x0 = lower_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x0, xend = x0, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "black") 
  }
  if(!is.na(upper_gate[1])){
    gate_line <- data.frame(Sample = plot_title,  x1 = upper_gate)
    gg <- gg + ggplot2::geom_segment(data = gate_line, ggplot2::aes(x = x1, xend = x1, y = as.numeric(Sample), yend = as.numeric(Sample) + 0.9), color = "black") 
  }
  
  if(!is.na(xlim)[1]){
    gg <- gg + ggplot2::coord_cartesian(xlim = xlim) 
  }
  
  return(gg)
}



#' signal_signal_just_plot, scatterplot of two different signals
#' @param data, transformed data 
#' @param channel1, which channel to plot
#' @param channel2, which channel to plot
#' @param plot_title, vector with title for each plot, default NA where the plots are numbered 1, 2, 3, etc.
#' @param xlim, xlim default NA.
#' @param ylim, ylim default NA.
#' @return scatterplots of two differnt signals per file. 


signal_signal_just_plot <- function(data, random_cells, channel1, channel2, 
                                    plot_title = NA, xlim = NA, ylim = NA){
  channel1 <- ggplot2::sym(channel1)
  channel2 <- ggplot2::sym(channel2)
  columnVar1 <- which(colnames(data[[1]]) == channel1)
  columnVar2 <- which(colnames(data[[1]]) == channel2)
  
  
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
    
    gg <- ggplot2::ggplot(data[[i]][random_cells[[i]],], ggplot2::aes(x=!!channel1, y=!!channel2)) +
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






#' density_plot_per_cluster, plot density for all markers in each cluster
#' @param data, transformed data 
#' @param clusters, vector of clusters
#' @param rand_cells, which cells to plot
#' @param plot_cluster, default NA give alle clusters, else vector of clusters to plot
#' @param nrow_plot, tells how many rows of plots to export, default = 4,
#' @param strip_text_size, give size of stripe text, defaut = 10, 
#' @param legend_text_size, give size of legend text, defaut = 10, 
#' @param axis_text_size, give size of axis text, defaut = 8
#' @return density plots for clusters against rest of data. 



density_plot_per_cluster <- function(data, cluster_per_cell, rand_cells = NA, plot_cluster = NA, nrow_plot = 4, strip_text_size = 10, legend_text_size = 10, axis_text_size = 8){
  if(!is.na(rand_cells[1])){
    data <- data[rand_cells,]
    cluster_per_cell <- cluster_per_cell[rand_cells]
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


#' barplot_per_sample, barplot of cells in each cluster, per sample
#' @param file_names_per_cell, transformed data 
#' @param cluster_per_cell, vector of clusters
#' @param rand_cells, which cells to plot
#' @return density plots for clusters against rest of data. 

barplot_per_sample <- function(file_names_per_cell, cluster_per_cell, rand_cells = NA){
  if(!is.na(rand_cells[1])){
    file_names_per_cell <- file_names_per_cell[rand_cells]
    cluster_per_cell <- cluster_per_cell[rand_cells]
  }  
  
  n_clusters <- length(unique(cluster_per_cell))
  tab <- table( cluster_per_cell, file_names_per_cell) 

  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  barplot(tab, col = col25, las = 2)
  legend("topright", inset=c(-0.1,0), col= col25[n_clusters:1], legend = names(col25[n_clusters:1]), pch =15)
}
