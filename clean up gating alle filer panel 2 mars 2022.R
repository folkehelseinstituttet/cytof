# inspiration:
# https://www.fluidigm.com/binaries/content/documents/fluidigm/search/hippo%3Aresultset/approach-to-bivariate-analysis-of-data-acquired-using-the-maxpar-direct-immune-profiling-assay-technical-note/fluidigm%3Afile



# Example
# 
# 
ptm <- proc.time()
# path to where the fcs files are stored F:\Forskningsprosjekter\PDB 2794 - Immune responses aga_\Forskningsfiler\JOBO\CyTOF\Datafiles\Panel 1 all files
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Datafiles", "Panel 2", "Panel 2 all files")
outDataPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "clean data")
scriptPath <- fs::path("H:", "git", "cytof")
out_result <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022")
outFigSignalPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "FigSignal")
outFigDensityPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "FigDensity")

source(fs::path(scriptPath, "gating_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))
source(fs::path(scriptPath, "analysis_functions.R"))


fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
fcs_files


files_to_open <- basename(fcs_files)
files_to_open <- files_to_open[grepl(".fcs", files_to_open)]
setwd(dirname(fcs_files[1]))
file_names <- gsub(".fcs", "", files_to_open)
# Read the files into a flowset


n_files <- length(file_names)
filenumbers <- 1:n_files

percent_lost_each_gating <- as.data.frame(matrix(NA, ncol = 9, nrow = n_files))
percent_lost_from_full_dataset <-   as.data.frame(matrix(NA, ncol = 9, nrow = n_files))


colnames(percent_lost_each_gating) <- c("Ce140Di", "Residual", "Center", "Offset", "Width",
                                        "Event_length", "Pt194Di", "Ir191Di", "Ir193Di")

colnames(percent_lost_from_full_dataset) <- colnames(percent_lost_each_gating)

rownames(percent_lost_each_gating) <- file_names
rownames(percent_lost_from_full_dataset) <- rownames(percent_lost_each_gating)

plotSignal <- function(plot_list){
g <- gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]],
             plot_list[[3]], plot_list[[4]], plot_list[[5]],
             plot_list[[6]],  plot_list[[7]], plot_list[[8]],
             plot_list[[9]], plot_list[[10]], plot_list[[11]],
             plot_list[[12]], plot_list[[13]], plot_list[[14]],
             plot_list[[15]], plot_list[[16]], plot_list[[17]],
             plot_list[[18]], plot_list[[19]], plot_list[[20]],
             plot_list[[21]], plot_list[[22]], plot_list[[23]],
             plot_list[[24]], plot_list[[25]], plot_list[[26]],
             plot_list[[27]], plot_list[[28]], plot_list[[29]],
             plot_list[[30]], plot_list[[31]], plot_list[[32]],
             plot_list[[33]], plot_list[[34]], plot_list[[35]],
             plot_list[[36]], plot_list[[37]], plot_list[[38]],
             plot_list[[39]], plot_list[[40]], plot_list[[41]],
             plot_list[[42]], plot_list[[43]], plot_list[[44]],
             plot_list[[45]], plot_list[[46]], plot_list[[47]],
             plot_list[[48]], plot_list[[49]], plot_list[[50]],
             plot_list[[51]], plot_list[[52]], plot_list[[53]],
             plot_list[[54]], plot_list[[55]], plot_list[[56]],
             plot_list[[57]], plot_list[[58]], plot_list[[59]],
             plot_list[[60]], plot_list[[61]], plot_list[[62]],
             plot_list[[63]], plot_list[[64]], plot_list[[65]],
             plot_list[[66]], plot_list[[67]],  ncol = 12, nrow = 6)
   return(g)
}


  # read all files in data_path into one dataset fcs_data
  fcs_data_with_info <- read_some_data_from_folder(data_path, file_number = filenumbers)
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
  time_signal_plots <- time_signal_plot(data = beads_data, random_events = random_events_for_plotting, 
                                        channel = "Ce140Di",  plot_title = file_names, upper_gate = upper_gate_Ce140Di)
  
  # remove # for those lines that you want to use. 
  
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  # if you want to save all plots for later evaluation:
  tiff(fs::path(outFigSignalPath, "Signal_fig1_bead_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()

  density_plots <- density_plot(data = beads_data, channel = "Ce140Di", plot_title = file_names, upper_gate = upper_gate_Ce140Di, maksCellsUsed = 25000)
  #density_plots # to see the plots
  
  
  tiff(fs::path(outFigDensityPath, "fig1_bead_gating.tiff"), width = 1200, height = 2000)
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
  percent_lost_each_gating[as.character(file_names),"Ce140Di"] <- number_of_events_after_gating/number_of_events_raw_data * 100 
  percent_lost_from_full_dataset[as.character(file_names),"Ce140Di"] <- number_of_events_after_gating/number_of_events_raw_data * 100 
  
  
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
  residual_gates$upper_gates[residual_gates$upper_gates < 1] <- mean(residual_gates$upper_gates[residual_gates$upper_gates > 1])
  density_plots <- density_plot(data = clean_up_data, "Residual", plot_title = file_names, lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Residual", plot_title = file_names,  lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  # if you want to save all plots for later evaluation:
  tiff(fs::path(outFigSignalPath, "Signal_fig2_residual_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()

  tiff(fs::path(outFigDensityPath, "fig2_residual_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()

  events_to_keep_after_gating <- events_to_keep(data = clean_up_data, channel = "Residual", lower_gate = residual_gates$lower_gates, 
                                                upper_gate = residual_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_gating, file_names = file_names)
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_gating)
  number_of_events_after_residual_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Residual"] <- number_of_events_after_residual_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Residual"] <- number_of_events_after_residual_gating/number_of_events_before_residual_gating * 100 #percent remaining from bead gating
  
  
  
  #************************************************
  #gating on Center+ ----
  #************************************************
  number_of_events_before_center_gating <- number_of_events_after_residual_gating
  random_events_for_plotting <- random_events(number_of_events_before_center_gating)
  
  #update lower_gate_percent, upper_gate_percent
  center_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Center", lower_gate_percent = 10, upper_gate_percent = 10)
  density_plots <- density_plot(data = clean_up_data, "Center", plot_title = file_names, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Center", plot_title = file_names, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  # if you want to save all plots for later evaluation:
  tiff(fs::path(outFigSignalPath, "Signal_fig3_center_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  

  tiff(fs::path(outFigDensityPath, "fig3_center_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  events_to_keep_after_center_gating <- events_to_keep(data = clean_up_data, channel = "Center", lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)

  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_center_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_center_gating)
  number_of_events_after_center_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Center"] <- number_of_events_after_center_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Center"] <-number_of_events_after_center_gating/number_of_events_before_center_gating * 100 
  
  
  
  #************************************************
  #gating on Offset+ ----
  #************************************************
  number_of_events_before_offset_gating <-  number_of_events_after_center_gating
  random_events_for_plotting <- random_events(number_of_events_before_offset_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  offset_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Offset", lower_gate_percent = 25, upper_gate_percent = 25)
 
  
  
  
  density_plots <- density_plot(data = clean_up_data, "Offset", plot_title = file_names, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  tiff(fs::path(outFigDensityPath, "fig4_offset_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  events_to_keep_after_offset_gating <- events_to_keep(data = clean_up_data, channel = "Offset",  lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
  
  
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Offset", plot_title = file_names, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  # if you want to save all plots for later evaluation:
  tiff(fs::path(outFigSignalPath, "Signal_fig4_offset_gating.tiff"), width = 1800, height = 1200)
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
  random_events_for_plotting <- random_events(number_of_events_before_width_gating)
  
  
  
  #update lower_gate_percent, upper_gate_percent
  width_gates <- find_gaussian_gates_highest_top(data = clean_up_data, channel = "Width", lower_gate_percent = 7, upper_gate_percent = 7)
  
  density_plots <- density_plot(data = clean_up_data, "Width", plot_title = file_names, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  
  
  tiff(fs::path(outFigDensityPath, "fig5_width_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Width", plot_title = file_names, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  tiff(fs::path(outFigSignalPath, "Signal_fig5_width_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  events_to_keep_after_width_gating <- events_to_keep(data = clean_up_data, channel = "Width", lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_width_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_width_gating)
  number_of_events_after_width_gating <- number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Width"] <- number_of_events_after_width_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Width"] <- number_of_events_after_width_gating/number_of_events_before_width_gating * 100 
  
  
  
  
  #************************************************
  #gating on Event+ ----
  #************************************************
  number_of_events_before_event_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_event_gating)
  
  #update lower_gate_percent, upper_gate_percent
  EventGates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Event_length", lower_gate_percent = 20, upper_gate_percent = 20)
  EventGates$lower_gates[is.na(EventGates$lower_gates)] <- mean(EventGates$lower_gates[!is.na(EventGates$lower_gates)])
  
  density_plots <- density_plot(data = clean_up_data, "Event_length", plot_title = file_names, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  
  
  tiff(fs::path(outFigDensityPath, "fig6_event_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  time_signal_plots <- time_signal_plot(data =  clean_up_data, random_events = random_events_for_plotting, channel = "Event_length", plot_title = file_names, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  tiff(fs::path(outFigSignalPath, "Signal_fig6_event_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  
  
  events_to_keep_after_event_gating <- events_to_keep(data = clean_up_data, channel = "Event_length",  upper_gate = EventGates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_event_gating, file_names = file_names)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_event_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_event_gating)
  number_of_events_after_event_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Event_length"] <- number_of_events_after_event_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Event_length"] <- number_of_events_after_event_gating/number_of_events_before_event_gating * 100 
  
  
  
  
  #************************************************
  #gating on Live Dead----  Cis, 
  #************************************************
  number_of_events_before_cis_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_cis_gating)
  
  
  #update lower_gate_percent, upper_gate_percent
  cis_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Pt194Di", lower_gate_percent = 15, upper_gate_percent = 40)
 
  cis_gates$lower_gates[file_names == "FHI0016_T1_P2_01_2"] <- 3.25
  cis_gates$upper_gates[file_names == "FHI0016_T1_P2_01_2"] <- 5
  cis_gates$lower_gates[file_names == "FHI005_T1_P2_01_1"] <- 2.5
  cis_gates$lower_gates[file_names == "K_K1560_Fe_A_78_T1_P2_01_1"] <- 3.4
  cis_gates$upper_gates[file_names == "M_FHI004_Ma_B_43_T1_P2_01_2"] <- 4.4
  cis_gates$lower_gates[file_names == "M_FHI100T1_Ma_A_26_T1_P2_01_4"] <- 3.4
  cis_gates$lower_gates[file_names == "M_FHI155T1_Fe_B_20_T1_P2_01_4"] <- 3.5
  
  
  density_plots <- density_plot(data = clean_up_data, "Pt194Di", plot_title = file_names, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  
  
  tiff(fs::path(outFigDensityPath, "fig7_cis_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Pt194Di", plot_title = file_names, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
 
  tiff(fs::path(outFigSignalPath, "Signal_fig7_cis_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  
  
    events_to_keep_after_cis_gating <- events_to_keep(data = clean_up_data, channel = "Pt194Di",  lower_gate = cis_gates$lower_gate,  upper_gate = cis_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_cis_gating, file_names = file_names)
  
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_cis_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_cis_gating)
  number_of_events_after_cis_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Pt194Di"] <- number_of_events_after_cis_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Pt194Di"] <- number_of_events_after_cis_gating/number_of_events_before_cis_gating * 100 
  
  
  
  #************************************************
  #gating on DNA1, Ir191Di+ ----
  #************************************************
  number_of_events_before_Ir191Di_gating <-  number_of_events(fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_Ir191Di_gating)
  
  #update lower_gate_percent, upper_gate_percent
  #Ir191di_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Ir191Di", lower_gate_percent = NA, upper_gate_percent = NA, perc_included = 0.99995, main_top_to_left = F)
  Ir191di_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Ir191Di", lower_gate_percent = 15, upper_gate_percent = 15)
 
  Ir191di_gates$lower_gates[file_names == "FHI0016_T1_P2_01_2"] <- 4.5
  Ir191di_gates$upper_gates[file_names == "FHI0016_T1_P2_01_2"] <- 5.4
  Ir191di_gates$lower_gates[file_names == "FHI005_T1_P2_01_1"] <- 3.9
  Ir191di_gates$upper_gates[file_names == "FHI005_T1_P2_01_1"] <- 4.6
  
   # 
  density_plots <- density_plot(data = clean_up_data, "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate, maksCellsUsed = 25000)
  #density_plots
  
  
   
  tiff(fs::path(outFigDensityPath, "fig8_Ir191_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  tiff(fs::path(outFigSignalPath, "Signal_fig8_Ir191_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  
  
  events_to_keep_after_Ir191Di_gating <- events_to_keep(data = clean_up_data, channel = "Ir191Di",  lower_gate = Ir191di_gates$lower_gate,  upper_gate = Ir191di_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_Ir191Di_gating, file_names = file_names)
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_Ir191Di_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_Ir191Di_gating)
  number_of_events_after_Ir191Di_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Ir191Di"] <- number_of_events_after_Ir191Di_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Ir191Di"] <- number_of_events_after_Ir191Di_gating/number_of_events_before_Ir191Di_gating * 100 
  
  
  
  #************************************************
  #gating on DNA1, Ir193Di+ ----
  #************************************************
  number_of_events_before_Ir193Di_gating <-  number_of_events(fcs_data, file_names = file_names)
  random_events_for_plotting <- random_events(number_of_events_before_Ir193Di_gating)
  
  #update lower_gate_percent, upper_gate_percent
  Ir193di_gates <- find_gaussian_gates_highest_top(data = clean_up_data, channel = "Ir193Di", lower_gate_percent = 20, upper_gate_percent = 20)
  
  
  density_plots <- density_plot(data = clean_up_data, "Ir193Di", plot_title = file_names, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate, maksCellsUsed = 25000)
  
  
  tiff(fs::path(outFigDensityPath, "fig9_Ir193_gating.tiff"), width = 1200, height = 2000)
  print(density_plots)
  dev.off()
  time_signal_plots <- time_signal_plot(data = clean_up_data, random_events = random_events_for_plotting, channel = "Ir193Di", plot_title = file_names, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate)
  #time_signal_plots # to see all plots
  #time_signal_plots[1] # to see first plot
  
  
  tiff(fs::path(outFigSignalPath, "Signal_fig9_Ir193_gating.tiff"), width = 1800, height = 1200)
  plotSignal(plot_list = time_signal_plots)
  dev.off()
  
  events_to_keep_after_Ir193Di_gating <- events_to_keep(data = clean_up_data, channel = "Ir193Di",  lower_gate = Ir193di_gates$lower_gate,  upper_gate = Ir193di_gates$upper_gate)
  #percent_to_keep_this_gating(kept_events = events_to_keep_after_Ir193Di_gating, file_names = file_names)
  
  
  fcs_data <- update_data_based_on_events_to_keep(data = fcs_data, kept_events = events_to_keep_after_Ir193Di_gating)
  clean_up_data <- update_data_based_on_events_to_keep(data = clean_up_data, kept_events = events_to_keep_after_Ir193Di_gating)
  number_of_events_after_Ir193Di_gating <-  number_of_events(data = fcs_data, file_names = file_names)
  percent_lost_from_full_dataset[as.character(file_names),"Ir193Di"] <- number_of_events_after_Ir193Di_gating/number_of_events_raw_data * 100 #percent remaining from total
  percent_lost_each_gating[as.character(file_names),"Ir193Di"] <- number_of_events_after_Ir193Di_gating/number_of_events_before_Ir193Di_gating * 100 
  
  
  #************************************************
  #save fcs_data ----
  #************************************************
  
  flowCore::write.flowSet(fcs_data, outdir = outDataPath, filename = as.character(file_names))
  


write.csv2(percent_lost_each_gating, fs::path(out_result, "percent_kept_each_gating.csv"))
write.csv2(percent_lost_from_full_dataset, fs::path(out_result, "percent_kept_from_full_dataset.csv"))

proc.time() - ptm