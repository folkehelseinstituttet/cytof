# inspiration:
# https://www.fluidigm.com/binaries/content/documents/fluidigm/search/hippo%3Aresultset/approach-to-bivariate-analysis-of-data-acquired-using-the-maxpar-direct-immune-profiling-assay-technical-note/fluidigm%3Afile



# Example
# 
# 
ptm <- proc.time()
# path to where the fcs files are stored F:\Forskningsprosjekter\PDB 2794 - Immune responses aga_\Forskningsfiler\JOBO\CyTOF\Datafiles\Panel 1 all files
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Datafiles", "Panel 1 all files")
outDataPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Gating", "Gating fra R_FINAL_CISnedre")
scriptPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "gating")
out_result <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "gating", "gating_results_FINAL_CISnedre")


outFigPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "gating", "gating_results_FINAL_CISnedre", "Fig")

source(fs::path(scriptPath, "gating_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))

fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
fcs_files

tvungetLavereCISgate <- 0.5 #sett inn NA hvis du heller vil bruke gaussian gate for CIS. evt annet tall..

files_to_open <- basename(fcs_files)
files_to_open <- files_to_open[grepl(".fcs", files_to_open)]
setwd(dirname(fcs_files[1]))
file_names <- gsub(".fcs", "", files_to_open)
# Read the files into a flowset


n_files <- length(file_names)
filenumber <- 1:n_files

percent_loss_each_gating <- as.data.frame(matrix(NA, ncol = 9, nrow = n_files))
percent_loss_from_full_dataset <-   as.data.frame(matrix(NA, ncol = 9, nrow = n_files))

colnames(percent_loss_each_gating) <- c("Ce140Di", "Residual", "Center", "Offset", "Width",
                                        "Event_length", "Pt194Di", "Ir191Di", "Ir193Di")

colnames(percent_loss_from_full_dataset) <- c("Ce140Di", "Residual", "Center", "Offset", "Width",
                                              "Event_length", "Pt194Di", "Ir191Di", "Ir193Di")

rownames(percent_loss_each_gating) <- file_names
rownames(percent_loss_from_full_dataset) <- file_names

for(ii in 1:floor(n_files/6)){
  
  i <- filenumber[((ii-1) * 6 + (1:6))]
  i <- i[!is.na(i)]
  # read all files in data_path into one dataset fcs_data
  fcs_data_with_info <- read_some_data_from_folder(data_path, file_number = i)
  fcs_data <- fcs_data_with_info$fcs_data
  file_names <- factor(fcs_data_with_info$file_names, levels = fcs_data_with_info$file_names)
  rm(fcs_data_with_info)
  
  # get the parameters of fcs_data and store them in params.
  params <- get_params_fcs_data(fcs_data[[1]])
  
  
  #***************************************************
  #gating on Beads ----
  #***************************************************
  
  
  bead_channels <- as.character(params$name[grepl("140|151|153|165|175", params$name)])
  
  beads_data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = bead_channels)
  
  number_of_events_raw_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_raw_data, n = 10000)
  
  
  #find upper_gate for noise gating based on Ce140Di
  upper_gate_Ce140Di <- find_gate_perc_height_upper_noise(data = beads_data, channel =  "Ce140Di", upper_perc_height = 0.001)
  #time_signal_plots <- time_signal_plot(data = beads_data, random_events = random_events_for_plotting, 
  #                                      channel = "Ce140Di",  plot_title = file_names, upper_gate = upper_gate_Ce140Di)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  density_plots <- density_plot(data = beads_data, channel = "Ce140Di", plot_title = file_names, upper_gate = upper_gate_Ce140Di)
  #density_plots # to see the plots
  
  
  tiff(fs::path(outFigPath, paste0("fig1_bead_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  
  
  
  ## use the upper_gates found for gating. Either on one of the beads or why not all. 
  #"Ce140Di"
  events_to_keep_after_gating <- events_to_keep(data = beads_data, channel = "Ce140Di",  upper_gate = upper_gate_Ce140Di)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_gating, file_names = file_names)
  
  
  
  #overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_gating)
  beads_data <- update_data_based_on_events_to_keep(data = beads_data, kept_events = events_to_keep_after_gating)
  number_of_events_after_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_after_gating)
  percent_loss_each_gating[as.character(file_names),"Ce140Di"] <- number_of_events_after_gating/number_of_events_raw_data * 100 
  percent_loss_from_full_dataset[as.character(file_names),"Ce140Di"] <- number_of_events_after_gating/number_of_events_raw_data * 100 
  
  
  number_of_events_after_beads_gating <- number_of_events_after_gating
  #finished working with beads_data, remove to keep space in memory
  rm(beads_data)
  
  
  
  
  #************************************************
  #clean_up_data
  #************************************************
  clean_up_channels <- as.character(params$name[grep("Center|Offset|Width|Residual|Event|Ir191|Ir193|Pt195Di|Pt194Di", params$name)])  #Pt194Di og 195 tilsvarer Cis
  number_of_events_before_clean_up_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  clean_up_data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = clean_up_channels)
  
  #************************************************
  #gating on Residual+ ----
  #************************************************
  number_of_events_before_residual_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_residual_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  residual_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Residual", lower_gate_percent = 10, upper_gate_percent = 10)
  density_plots <- density_plot(data = clean_up_data, "Residual", plot_title = file_names, lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate)
  #density_plots
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Residual", plot_title = file_names,  lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  tiff(fs::path(outFigPath, paste0("fig2_residual_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  
  events_to_keep_after_gating <- events_to_keep(data = clean_up_data, channel = "Residual", lower_gate = residual_gates$lower_gates, 
                                                upper_gate = residual_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_gating, file_names = file_names)
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_gating)
  number_of_events_after_residual_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Residual"] <- number_of_events_after_residual_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Residual"] <- number_of_events_after_residual_gating/number_of_events_after_beads_gating * 100 #percent remaining from bead gating
  
  
  
  #************************************************
  #gating on Center+ ----
  #************************************************
  number_of_events_before_center_gating <- number_of_events_after_residual_gating
  random_events_for_plotting <- random_events(number_of_events_before_center_gating)
  
  #update lower_gate_percent, upper_gate_percent
  center_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Center", lower_gate_percent = 10, upper_gate_percent = 10)
  density_plots <- density_plot(data = clean_up_data, "Center", plot_title = file_names, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
  #density_plots
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Center", plot_title = file_names, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  tiff(fs::path(outFigPath, paste0("fig3_event_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  events_to_keep_after_center_gating <- events_to_keep(data = clean_up_data, channel = "Center", lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_center_gating, file_names = file_names)
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_center_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_center_gating)
  number_of_events_after_center_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Center"] <- number_of_events_after_center_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Center"] <-number_of_events_after_center_gating/number_of_events_after_residual_gating * 100 #percent remaining from residual gating
  
  
  
  #************************************************
  #gating on Offset+ ----
  #************************************************
  number_of_events_before_offset_gating <-  number_of_events_after_center_gating
  random_events_for_plotting <- random_events(number_of_events_before_offset_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  offset_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Offset", lower_gate_percent = 7, upper_gate_percent = 7)
  density_plots <- density_plot(data = clean_up_data, "Offset", plot_title = file_names, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
  #density_plots
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Offset", plot_title = file_names, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  tiff(fs::path(outFigPath, paste0("fig4_offset_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  events_to_keep_after_offset_gating <- events_to_keep(data = clean_up_data, channel = "Offset",  lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_offset_gating, file_names = file_names)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_offset_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_offset_gating)
  number_of_events_after_offset_gating <- number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Offset"] <- number_of_events_after_offset_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Offset"] <- number_of_events_after_offset_gating/number_of_events_after_center_gating * 100 #percent remaining from  center gating
  
  #************************************************
  #gating on Width+ ----
  #************************************************
  number_of_events_before_width_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_width_gating)
  
  
  
  #update lower_gate_percent, upper_gate_percent
  width_gates <- find_gaussian_gates_highest_top(data = clean_up_data, channel = "Width", lower_gate_percent = 2.5, upper_gate_percent = 4)
  
  density_plots <- density_plot(data = clean_up_data, "Width", plot_title = file_names, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
  #density_plots
  
  
  tiff(fs::path(outFigPath, paste0("fig5_width_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Width", plot_title = file_names, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  events_to_keep_after_width_gating <- events_to_keep(data = clean_up_data, channel = "Width", lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_width_gating, file_names = file_names)
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_width_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_width_gating)
  number_of_events_after_width_gating <- number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Width"] <- number_of_events_after_width_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Width"] <- number_of_events_after_width_gating/number_of_events_after_offset_gating * 100 #percent remaining from offset gating
  
  
  
  
  #************************************************
  #gating on Event+ ----
  #************************************************
  number_of_events_before_event_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_event_gating)
  
  #update lower_gate_percent, upper_gate_percent
  EventGates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Event_length", lower_gate_percent = 20, upper_gate_percent = 20)
  density_plots <- density_plot(data = clean_up_data, "Event_length", plot_title = file_names, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate)
  #density_plots
  
  
  tiff(fs::path(outFigPath, paste0("fig6_event_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  #time_signal_plots <- time_signal_plot(data =  clean_up_data, random_events = random_events_for_plotting, channel = "Event_length", plot_title = file_names, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  events_to_keep_after_event_gating <- events_to_keep(data = clean_up_data, channel = "Event_length",  upper_gate = EventGates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_event_gating, file_names = file_names)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_event_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_event_gating)
  number_of_events_after_event_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Event_length"] <- number_of_events_after_event_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Event_length"] <- number_of_events_after_event_gating/number_of_events_after_width_gating * 100 #percent remaining from width gating
  
  
  
  
  #************************************************
  #gating on Live Dead----  Cis, 
  #************************************************
  number_of_events_before_cis_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_cis_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  cis_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Pt194Di", lower_gate_percent = 5, upper_gate_percent = 30)
  if(!is.na(tvungetLavereCISgate)){
    cis_gates$lower_gates <- rep(tvungetLavereCISgate, length(cis_gates$lower_gates))
  }
  
  density_plots <- density_plot(data = clean_up_data, "Pt194Di", plot_title = file_names, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate)
  #density_plots
  
  
  tiff(fs::path(outFigPath, paste0("fig7_cis_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Pt194Di", plot_title = file_names, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Pt194Di", plot_title = file_names, lower_gate = NA, upper_gate = NA)
  #time_signal_plots # to see all plots
  events_to_keep_after_cis_gating <- events_to_keep(data = clean_up_data, channel = "Pt194Di",  lower_gate = cis_gates$lower_gate,  upper_gate = cis_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_cis_gating, file_names = file_names)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_cis_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_cis_gating)
  number_of_events_after_cis_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Pt194Di"] <- number_of_events_after_cis_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Pt194Di"] <- number_of_events_after_cis_gating/number_of_events_after_event_gating * 100 #percent remaining from event gating
  
  
  
  #************************************************
  #gating on DNA1, Ir191Di+ ----
  #************************************************
  number_of_events_before_Ir191Di_gating <-  number_of_events(fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_Ir191Di_gating)
  
  #update lower_gate_percent, upper_gate_percent
  # Ir191di_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Ir191Di", lower_gate_percent = NA, upper_gate_percent = NA, perc_included = 0.99995, main_top_to_left = F)
  # 
  # density_plots <- density_plot(data = clean_up_data, "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate)
  #density_plots
  
  
  Ir191di_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Ir191Di", lower_gate_percent = NA, upper_gate_percent = NA, perc_included = 0.99, main_top_to_left = F)
  
  density_plots <- density_plot(data = clean_up_data, "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate)
  
  
  tiff(fs::path(outFigPath, paste0("fig8_Ir191_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  events_to_keep_after_Ir191Di_gating <- events_to_keep(data = clean_up_data, channel = "Ir191Di",  lower_gate = Ir191di_gates$lower_gate,  upper_gate = Ir191di_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_Ir191Di_gating, file_names = file_names)
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_Ir191Di_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_Ir191Di_gating)
  number_of_events_after_Ir191Di_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Ir191Di"] <- number_of_events_after_Ir191Di_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Ir191Di"] <- number_of_events_after_Ir191Di_gating/number_of_events_after_cis_gating * 100 #percent remaining from cis gating
  
  
  
  #************************************************
  #gating on DNA1, Ir193Di+ ----
  #************************************************
  number_of_events_before_Ir193Di_gating <-  number_of_events(fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_Ir193Di_gating)
  
  #update lower_gate_percent, upper_gate_percent
  Ir193di_gates <- find_gaussian_gates_highest_top(data = clean_up_data, channel = "Ir193Di", lower_gate_percent = 15, upper_gate_percent = 5)
  
  
  density_plots <- density_plot(data = clean_up_data, "Ir193Di", plot_title = file_names, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate)
  
  
  tiff(fs::path(outFigPath, paste0("fig9_Ir193_gating", ii, ".tiff")))
  print(density_plots)
  dev.off()#density_plots
  #time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Ir193Di", plot_title = file_names, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  events_to_keep_after_Ir193Di_gating <- events_to_keep(data = clean_up_data, channel = "Ir193Di",  lower_gate = Ir193di_gates$lower_gate,  upper_gate = Ir193di_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_Ir193Di_gating, file_names = file_names)
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_Ir193Di_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_Ir193Di_gating)
  number_of_events_after_Ir193Di_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_loss_from_full_dataset[as.character(file_names),"Ir193Di"] <- number_of_events_after_Ir193Di_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_loss_each_gating[as.character(file_names),"Ir193Di"] <- number_of_events_after_Ir193Di_gating/number_of_events_after_cis_gating * 100 #percent remaining from cis gating
  
  
  #************************************************
  #save fcs_data ----
  #************************************************
  
  flowCore::write.flowSet(fcs_data, outdir = outDataPath, filename = as.character(file_names))
  print(".")
}

write.csv2(percent_loss_each_gating, fs::path(out_result, "percent_kept_each_gatingnr2.csv"))
write.csv2(percent_loss_from_full_dataset, fs::path(out_result, "percent_kept_from_full_datasetnr2.csv"))

proc.time() - ptm