# inspiration:
# https://www.fluidigm.com/binaries/content/documents/fluidigm/search/hippo%3Aresultset/approach-to-bivariate-analysis-of-data-acquired-using-the-maxpar-direct-immune-profiling-assay-technical-note/fluidigm%3Afile


#read https://www.vsh.com/publication/PathsetterCleanup.pdf


# Example
# 
# 
ptm <- proc.time()

source(fs::path(paths$script_path,  "functions.R"))

fcs_files <- fs::path(paths$raw_data_path, rownames(file.info(list.files(paths$raw_data_path))))
fcs_files


files_to_open <- basename(fcs_files)
files_to_open <- files_to_open[grepl(".fcs", files_to_open)]
setwd(dirname(fcs_files[1]))
file_names <- gsub(".fcs", "", files_to_open)
file_names_in_plot <- 1:length(file_names)

n_files <- length(file_names)
filenumbers <- 1:n_files

percent_lost_each_gating <- as.data.frame(matrix(NA, ncol = 9, nrow = n_files))
percent_lost_from_full_dataset <-   as.data.frame(matrix(NA, ncol = 9, nrow = n_files))


colnames(percent_lost_each_gating) <- c("Ce140Di", "Residual", "Center", "Offset", "Width",
                                        "Event_length", "Pt194Di", "Ir191Di", "Ir193Di")

colnames(percent_lost_from_full_dataset) <- colnames(percent_lost_each_gating)

rownames(percent_lost_each_gating) <- file_names
rownames(percent_lost_from_full_dataset) <- rownames(percent_lost_each_gating)




  # read all files in data_path into one dataset fcs_data
  fcs_data_with_info <- read_data_from_folder(paths$raw_data_path)
  fcs_data <- fcs_data_with_info$fcs_data
  file_names <- factor(fcs_data_with_info$file_names, levels = fcs_data_with_info$file_names)
  rm(fcs_data_with_info)
  
  # get the parameters of fcs_data and store them in params_fcs.
  params_fcs <- get_params_fcs_data(fcs_data[[1]])
  
  
  #***************************************************
  #gating on Beads ----
  #***************************************************
  
  
  bead_markers <- as.character(params_fcs$name[grepl("140|151|153|165|175", params_fcs$name)])
  
  beads_data <-  transform_selected_markers(fcs_data = fcs_data, markers = bead_markers, method = "arc_sinh", cofactor = 5)
  
  number_of_events_raw_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_raw_data, n = 10000)
  
  
  #find upper_gate for noise gating based on Ce140Di
  upper_gate_Ce140Di <- find_gaussian_gates(data = beads_data, marker = "Ce140Di", method = "upper noise", upper_gate_percent = 0.01)
   
  upper_gate_Ce140Di <- rep(4, length( upper_gate_Ce140Di))
  
  
  ymax <- max_marker(data = beads_data, marker = "Ce140Di")
  ymin <- min_marker(data = beads_data, marker = "Ce140Di")
  
  time_signal_plots <- time_signal_plot(data = beads_data, random_events = random_events_for_plotting, 
                                        marker = "Ce140Di",  plot_title = file_names_in_plot, upper_gate = upper_gate_Ce140Di, ylim = c(ymin, ymax))
  
 
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig1_bead_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()

  density_plots <- density_plot(data = beads_data, marker = "Ce140Di", plot_title = file_names_in_plot, upper_gate = upper_gate_Ce140Di, max_event_used = 25000)
 
  
  tiff(fs::path(paths$clean_data_fig_density_path, "fig1_bead_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()

  
  
  ## use the upper_gates found for gating. Either on one of the beads or why not all. 
  #"Ce140Di"
  events_to_keep_after_gating <- events_to_keep(data = beads_data, marker = "Ce140Di",  upper_gate = upper_gate_Ce140Di)
   
  #overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_gating)
  beads_data <- update_data_based_on_events_to_keep(data = beads_data, kept_events = events_to_keep_after_gating)
  number_of_events_after_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_after_gating)
  percent_lost_each_gating[as.character(file_names),"Ce140Di"] <- number_of_events_after_gating/number_of_events_raw_data * 100 
  percent_lost_from_full_dataset[as.character(file_names),"Ce140Di"] <- number_of_events_after_gating/number_of_events_raw_data * 100 
  
  
  number_of_events_after_beads_gating <- number_of_events_after_gating
  #finished working with beads_data, remove to keep space in memory
  rm(beads_data)
  
  
  
  
  #************************************************
  #clean_up_data
  #************************************************
  clean_up_markers <- as.character(params_fcs$name[grep("Center|Offset|Width|Residual|Event|Ir191|Ir193|Pt195Di|Pt194Di", params_fcs$name)])  #Pt194Di og 195 tilsvarer Cis
  as.character(params_fcs$desc[grep("CD3|CD45", params_fcs$desc)])  #CD3 og CD45 are here chosen when plotting against a real marker.
  extra_markers <- as.character(params_fcs$name[grep("CD3|CD45", params_fcs$desc)[1:2]])  #keep only CD3 and CD45 (not CD45RA which also are found by grep as the third element on my list)
  CD3 <- extra_markers[2] #important to cheek that you got the correct marker, might have to be changed
  CD45 <- extra_markers[1]
  
  number_of_events_before_clean_up_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  clean_up_data <-  transform_selected_markers(fcs_data = fcs_data, markers = c(clean_up_markers, extra_markers), method = "arc_sinh", cofactor = 5)
  
  
 
  

  
  
  #************************************************
  #gating on Residual+ ----
  #************************************************
  number_of_events_before_residual_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_residual_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  residual_gates <- find_gaussian_gates(data = clean_up_data, marker = "Residual", method = "second top", lower_gate_percent = 25, upper_gate_percent = 25)
  density_plots <- density_plot(data = clean_up_data, "Residual", plot_title = file_names_in_plot, lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate, max_event_used = 25000)

  
  ymax <- max_marker(data = clean_up_data, marker = "Residual")
  ymin <- min_marker(data = clean_up_data, marker = "Residual")
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Residual", plot_title = file_names_in_plot,  lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate, ylim = c(ymin, ymax))
  
  
  # if you want to save all plots for later evaluation:
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig2_residual_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()

  tiff(fs::path(paths$clean_data_fig_density_path, "fig2_residual_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()

  events_to_keep_after_gating <- events_to_keep(data = clean_up_data, marker = "Residual", lower_gate = residual_gates$lower_gates, 
                                                upper_gate = residual_gates$upper_gate)
    
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_gating)
  number_of_events_after_residual_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Residual"] <- number_of_events_after_residual_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Residual"] <- number_of_events_after_residual_gating/number_of_events_before_residual_gating * 100 #percent remaining from bead gating
  
  
  
  #************************************************
  #gating on Center+ ----
  #************************************************
  number_of_events_before_center_gating <- number_of_events_after_residual_gating
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_center_gating)
  
  #update lower_gate_percent, upper_gate_percent
  center_gates <- find_gaussian_gates(data = clean_up_data, marker = "Center", method = "second top", lower_gate_percent = 25, upper_gate_percent = 25)
  density_plots <- density_plot(data = clean_up_data, "Center", plot_title = file_names_in_plot, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate, max_event_used = 25000)

  
  ymax <- max_marker(data = clean_up_data, marker = "Center")
  ymin <- min_marker(data = clean_up_data, marker = "Center")
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Center", plot_title = file_names_in_plot, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate, ylim = c(ymin, ymax))
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  # if you want to save all plots for later evaluation:
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig3_center_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  

  tiff(fs::path(paths$clean_data_fig_density_path, "fig3_center_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  events_to_keep_after_center_gating <- events_to_keep(data = clean_up_data, marker = "Center", lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)

  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_center_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_center_gating)
  number_of_events_after_center_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Center"] <- number_of_events_after_center_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Center"] <-number_of_events_after_center_gating/number_of_events_before_center_gating * 100 
  
  
  
  #************************************************
  #gating on Offset+ ----
  #************************************************
  number_of_events_before_offset_gating <-  number_of_events_after_center_gating
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_offset_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  offset_gates <- find_gaussian_gates(data = clean_up_data, marker = "Offset", method = "second top", lower_gate_percent = 25, upper_gate_percent = 25)
 
  #MANUALLY CHANGING SOME GATES
  # offset_gates$lower_gates[grep("T1_FHI003", file_names)] <- offset_gates$lower_gates[grep("T2_FHI003", file_names)] 
  # offset_gates$upper_gates[grep("T1_FHI003", file_names)] <- offset_gates$upper_gates[grep("T2_FHI003", file_names)] 
  # offset_gates$lower_gates[grep("T1_FHI016", file_names)] <- offset_gates$lower_gates[grep("T1_FHI81", file_names)] 
  # offset_gates$upper_gates[grep("T1_FHI016", file_names)] <- offset_gates$upper_gates[grep("T1_FHI81", file_names)] 
  # 
  
  
  density_plots <- density_plot(data = clean_up_data, "Offset", plot_title = file_names_in_plot, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate, max_event_used = 25000)

  tiff(fs::path(paths$clean_data_fig_density_path, "fig4_offset_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  events_to_keep_after_offset_gating <- events_to_keep(data = clean_up_data, marker = "Offset",  lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
  
  ymax <- max_marker(data = clean_up_data, marker = "Offset")
  ymin <- min_marker(data = clean_up_data, marker = "Offset")
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Offset", plot_title = file_names_in_plot, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate, ylim = c(ymin, ymax))
 
  # if you want to save all plots for later evaluation:
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig4_offset_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
   
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_offset_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_offset_gating)
  number_of_events_after_offset_gating <- number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Offset"] <- number_of_events_after_offset_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Offset"] <- number_of_events_after_offset_gating/number_of_events_before_offset_gating * 100 #percent remaining from  center gating
  
  #************************************************
  #gating on Width+ ----
  #************************************************
  number_of_events_before_width_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_width_gating)
  
  
  
  #update lower_gate_percent, upper_gate_percent
  width_gates <- find_gaussian_gates(data = clean_up_data, marker = "Width", method = "highest top", lower_gate_percent = 20, upper_gate_percent = 20)
  
  density_plots <- density_plot(data = clean_up_data, "Width", plot_title = file_names_in_plot, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate, max_event_used = 25000)
   
  
  tiff(fs::path(paths$clean_data_fig_density_path, "fig5_width_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  
  ymax <- max_marker(data = clean_up_data, marker = "Width")
  ymin <- min_marker(data = clean_up_data, marker = "Width")
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Width", plot_title = file_names_in_plot, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate, ylim = c(ymin, ymax))
 
  
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig5_width_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  events_to_keep_after_width_gating <- events_to_keep(data = clean_up_data, marker = "Width", lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_width_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_width_gating)
  number_of_events_after_width_gating <- number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Width"] <- number_of_events_after_width_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Width"] <- number_of_events_after_width_gating/number_of_events_before_width_gating * 100 
  
  
  
  
  #************************************************
  #gating on Event+ ----
  #************************************************
  number_of_events_before_event_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_event_gating)
  
  #update lower_gate_percent, upper_gate_percent
  EventGates <- find_gaussian_gates(data = clean_up_data, marker = "Event_length", method = "second top", lower_gate_percent = 20, upper_gate_percent = 20)
  density_plots <- density_plot(data = clean_up_data, "Event_length", plot_title = file_names_in_plot, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate, max_event_used = 25000)
  
  
  tiff(fs::path(paths$clean_data_fig_density_path, "fig6_event_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  
  ymax <- max_marker(data = clean_up_data, marker = "Event_length")
  ymin <- min_marker(data = clean_up_data, marker = "Event_length")
  
  time_signal_plots <- time_signal_plot(data =  clean_up_data, random_events = random_events_for_plotting, marker = "Event_length", plot_title = file_names_in_plot, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate, ylim = c(ymin, ymax))
  
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig6_event_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  
  
  events_to_keep_after_event_gating <- events_to_keep(data = clean_up_data, marker = "Event_length",  upper_gate = EventGates$upper_gate)
    
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_event_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_event_gating)
  number_of_events_after_event_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Event_length"] <- number_of_events_after_event_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Event_length"] <- number_of_events_after_event_gating/number_of_events_before_event_gating * 100 
  
  
  #
  #HER
  #
  #************************************************
  #gating on Live Dead----  Cis, 
  #************************************************
  number_of_events_before_cis_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_cis_gating)
  
  
  
  #update lower_gate_percent, upper_gate_percent
  cis_gates <- find_gaussian_gates(data = clean_up_data, marker = "Pt194Di", method = "second top", lower_gate_percent = 2, upper_gate_percent = 40)
  
  
  density_plots <- density_plot(data = clean_up_data, "Pt194Di", plot_title = file_names_in_plot, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate, max_event_used = 25000)
   
  tiff(fs::path(paths$clean_data_fig_density_path, "fig7_cis_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  
  ymax <- max_marker(data = clean_up_data, marker = "Pt194Di")
  ymin <- min_marker(data = clean_up_data, marker = "Pt194Di")
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Pt194Di", plot_title = file_names_in_plot, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate, ylim = c(ymin, ymax)) #time_signal_plots[1] # to see first plot
 
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig7_cis_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  
  
  events_to_keep_after_cis_gating <- events_to_keep(data = clean_up_data, marker = "Pt194Di",  lower_gate = cis_gates$lower_gate,  upper_gate = cis_gates$upper_gate)
   
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_cis_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_cis_gating)
  number_of_events_after_cis_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Pt194Di"] <- number_of_events_after_cis_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Pt194Di"] <- number_of_events_after_cis_gating/number_of_events_before_cis_gating * 100 
  
  
  
  #************************************************
  #gating on DNA1, Ir191Di+ ----
  #************************************************
  number_of_events_before_Ir191Di_gating <-  number_of_events(fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_Ir191Di_gating)
  
  #update lower_gate_percent, upper_gate_percent
 
#  Ir191di_gates <- find_gaussian_gates_second_top(data = clean_up_data, marker = "Ir191Di", lower_gate_percent = NA, upper_gate_percent = NA, perc_included = 0.99, main_top_to_left = F)
  Ir191di_gates <- find_gaussian_gates(data = clean_up_data, marker = "Ir191Di", method = "second top", lower_gate_percent = 25, upper_gate_percent = 25)
  
  
  
  density_plots <- density_plot(data = clean_up_data, "Ir191Di", plot_title = file_names_in_plot, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate, max_event_used = 25000)
   
  tiff(fs::path(paths$clean_data_fig_density_path, "fig8_Ir191_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  
  ymax <- max_marker(data = clean_up_data, marker = "Ir191Di")
  ymin <- min_marker(data = clean_up_data, marker = "Ir191Di")
  

  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Ir191Di", plot_title = file_names_in_plot, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate, ylim = c(ymin, ymax))
    
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig8_Ir191_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  
  
  events_to_keep_after_Ir191Di_gating <- events_to_keep(data = clean_up_data, marker = "Ir191Di",  lower_gate = Ir191di_gates$lower_gate,  upper_gate = Ir191di_gates$upper_gate)
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_Ir191Di_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_Ir191Di_gating)
  number_of_events_after_Ir191Di_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Ir191Di"] <- number_of_events_after_Ir191Di_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Ir191Di"] <- number_of_events_after_Ir191Di_gating/number_of_events_before_Ir191Di_gating * 100 
  
  
  
  #************************************************
  #gating on DNA1, Ir193Di+ ----
  #************************************************
  number_of_events_before_Ir193Di_gating <-  number_of_events(fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(numb_events = number_of_events_before_Ir193Di_gating)
  
  #update lower_gate_percent, upper_gate_percent
  Ir193di_gates <- find_gaussian_gates(data = clean_up_data, marker = "Ir193Di", method = "highest top", lower_gate_percent = 25, upper_gate_percent = 25)
  
  
  density_plots <- density_plot(data = clean_up_data, "Ir193Di", plot_title = file_names_in_plot, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate, max_event_used = 25000)
  
  
  tiff(fs::path(paths$clean_data_fig_density_path, "fig9_Ir193_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  
  
  
  ymax <- max_marker(data = clean_up_data, marker = "Ir193Di")
  ymin <- min_marker(data = clean_up_data, marker = "Ir193Di")
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, marker = "Ir193Di", plot_title = file_names_in_plot, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate, ylim = c(ymin, ymax))
    
  tiff(fs::path(paths$clean_data_fig_signal_path, "Signal_fig9_Ir193_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  events_to_keep_after_Ir193Di_gating <- events_to_keep(data = clean_up_data, marker = "Ir193Di",  lower_gate = Ir193di_gates$lower_gate,  upper_gate = Ir193di_gates$upper_gate)
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_Ir193Di_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_Ir193Di_gating)
  number_of_events_after_Ir193Di_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Ir193Di"] <- number_of_events_after_Ir193Di_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Ir193Di"] <- number_of_events_after_Ir193Di_gating/number_of_events_before_Ir193Di_gating * 100 
  
  
  #************************************************
  #save fcs_data ----
  #************************************************
  
  flowCore::write.flowSet(fcs_data, outdir = paths$clean_data_path, filename = as.character(file_names))
  


write.csv2(percent_lost_each_gating, fs::path(paths$clean_data_info_path, "percent_kept_each_gating.csv"))
write.csv2(percent_lost_from_full_dataset, fs::path(paths$clean_data_info_path, "percent_kept_from_full_dataset.csv"))


markers <- params_fcs$desc[grepl("_",params_fcs$desc)]
markers <- markers[!grepl("_DNA", markers)]
markers <- markers[!grepl("_Cis", markers)]
markers <- markers[!grepl("_Beads", markers)]

markers_name <- params_fcs$name[params_fcs$desc %in% markers]
print("number of markers included:")
print(length(markers_name))

d <- cbind( as.character(markers_name),  as.character(markers),  as.character(markers))
colnames(d) <- c("markers_name", "marker", "markers_short_name")

for(i in 1:nrow(d)){
  if(grepl("_", as.character(d[i,3]))){
    d[i,3] <- strsplit(as.character(d[i,3]), "_")[[1]][2]
  }
}
write.csv2(d, fs::path(paths$clean_data_info_path, "marker_names_included_manual_shortnames.csv"))


proc.time() - ptm