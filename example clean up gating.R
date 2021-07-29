# Example
# 
# 

# path to where the fcs files are stored
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2770 - EXIMIOUS - Mapping E_", "CyTOF in vitro WP4", "FCS-filer til Anja")
outDataPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2770 - EXIMIOUS - Mapping E_", "CyTOF in vitro WP4", "GatedData")
scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "gating.R"))
source(fs::path(scriptPath, "read_data.R"))
source(fs::path(scriptPath, "transformation.R"))


# read all files in data_path into one dataset fcs_data
fcs_data_with_info <- read_data_from_folder(data_path)
fcs_data <- fcs_data_with_info$fcs_data
file_names <- fcs_data_with_info$file_names
rm(fcs_data_with_info)

# get the parameters of fcs_data and store them in params.
params <- get_params_fcs_data(fcs_data[[1]])


#***************************************************
#gating on Beads ----
#***************************************************


bead_channels <- as.character(params$name[grepl("140|151|153|165|175", params$name)])

beads_data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = bead_channels)

number_of_cells_raw_data <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_raw_data)


#find upper_gate for noise gating based on Ce140Di
upper_gate_Ce140Di <- find_gate_perc_height_upper_noise(data = beads_data, channel =  "Ce140Di", upper_perc_height = 0.001)
time_signal_plots <- time_signal_plot(data = beads_data, random_cells = random_cells_for_plotting, 
                                        channel = "Ce140Di",  plot_title = file_names, upper_gate = upper_gate_Ce140Di)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot
density_plots <- density_plot(data = beads_data, channel = "Ce140Di", plot_title = file_names, upper_gate = upper_gate_Ce140Di)
#density_plots # to see the plots

#find upper_gate for noise gating based on Eu151Di
upper_gate_Eu151Di <- find_gate_perc_height_upper_noise(data = beads_data, channel =  "Eu151Di", upper_perc_height = 0.001)
time_signal_plots <- time_signal_plot(data = beads_data, random_cells = random_cells_for_plotting, 
                                      channel = "Eu151Di",  plot_title = file_names, upper_gate = upper_gate_Eu151Di)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot
density_plots <- density_plot(data = beads_data, channel = "Eu151Di", plot_title = file_names, upper_gate = upper_gate_Eu151Di)
#density_plots # to see the plots

#find upper_gate for noise gating based on Eu153Di
upper_gate_Eu153Di <- find_gate_perc_height_upper_noise(data = beads_data, channel =  "Eu153Di", upper_perc_height = 0.001)
time_signal_plots <- time_signal_plot(data = beads_data, random_cells = random_cells_for_plotting, 
                                      channel = "Eu153Di",  plot_title = file_names, upper_gate = upper_gate_Eu153Di)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot
density_plots <- density_plot(data = beads_data, channel = "Eu153Di", plot_title = file_names, upper_gate = upper_gate_Eu153Di)
#density_plots # to see the plots

#find upper_gate for noise gating based on Ho165Di
upper_gate_Ho165Di <- find_gate_perc_height_upper_noise(data = beads_data, channel =  "Ho165Di", upper_perc_height = 0.001)
time_signal_plots <- time_signal_plot(data = beads_data, random_cells = random_cells_for_plotting, 
                                      channel = "Ho165Di",  plot_title = file_names, upper_gate = upper_gate_Ho165Di)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot
density_plots <- density_plot(data = beads_data, channel = "Ho165Di", plot_title = file_names, upper_gate = upper_gate_Ho165Di)
#density_plots # to see the plots

#find upper_gate for noise gating based on Lu175Di
upper_gate_Lu175Di <- find_gate_perc_height_upper_noise(data = beads_data, channel =  "Lu175Di", upper_perc_height = 0.001)
time_signal_plots <- time_signal_plot(data = beads_data, random_cells = random_cells_for_plotting, 
                                      channel = "Lu175Di",  plot_title = file_names, upper_gate = upper_gate_Lu175Di)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot
density_plots <- density_plot(data = beads_data, channel = "Lu175Di", plot_title = file_names, upper_gate = upper_gate_Lu175Di)
#density_plots # to see the plots


## use the upper_gates found for gating. Either on one of the beads or why not all. 
#"Ce140Di"
cells_to_keep_after_gating <- cells_to_keep(data = beads_data, channel = "Ce140Di",  upper_gate = upper_gate_Ce140Di)
#overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_gating)
beads_data <- update_data_based_on_cells_to_keep(data = beads_data, kept_cells = cells_to_keep_after_gating)
number_of_cells_after_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_after_gating)
percent_remaining_from_total <- number_of_cells_after_gating/number_of_cells_raw_data * 100 
percent_remaining_from_total

#"Eu151Di"
cells_to_keep_after_gating <- cells_to_keep(data = beads_data, channel = "Eu151Di",  upper_gate = upper_gate_Eu151Di)
#overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_gating)
beads_data <- update_data_based_on_cells_to_keep(data = beads_data, kept_cells = cells_to_keep_after_gating)
number_of_cells_after_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_after_gating)
percent_remaining_from_total <- number_of_cells_after_gating/number_of_cells_raw_data * 100 
percent_remaining_from_total

#"Eu153Di"
cells_to_keep_after_gating <- cells_to_keep(data = beads_data, channel = "Eu153Di",  upper_gate = upper_gate_Eu153Di)
#overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_gating)
beads_data <- update_data_based_on_cells_to_keep(data = beads_data, kept_cells = cells_to_keep_after_gating)
number_of_cells_after_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_after_gating)
percent_remaining_from_total <- number_of_cells_after_gating/number_of_cells_raw_data * 100 
percent_remaining_from_total

#"Ho165Di"
cells_to_keep_after_gating <- cells_to_keep(data = beads_data, channel = "Ho165Di",  upper_gate = upper_gate_Ho165Di)
#overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_gating)
beads_data <- update_data_based_on_cells_to_keep(data = beads_data, kept_cells = cells_to_keep_after_gating)
number_of_cells_after_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_after_gating)
percent_remaining_from_total <- number_of_cells_after_gating/number_of_cells_raw_data * 100 
percent_remaining_from_total

#"Lu175Di"
cells_to_keep_after_gating <- cells_to_keep(data = beads_data, channel = "Lu175Di",  upper_gate = upper_gate_Lu175Di)
#overwrite the raw  and Beads datasett (will take lot of space if we make one new each time, and do not need it)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_gating)
beads_data <- update_data_based_on_cells_to_keep(data = beads_data, kept_cells = cells_to_keep_after_gating)
number_of_cells_after_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_after_gating)
percent_remaining_from_total <- number_of_cells_after_gating/number_of_cells_raw_data * 100 
percent_remaining_from_total


number_of_cells_after_beads_gating <- number_of_cells_after_gating
#finished working with beads_data, remove to keep space in memory
rm(beads_data)




#************************************************
#clean_up_data
#************************************************
clean_up_channels <- as.character(params$name[grep("Center|Offset|Width|Residual|Event|Ir191|Ir193|Pt195Di", params$name)])  #Pt195Di tilsvarer Cis
number_of_cells_before_clean_up_gating <-  number_of_cells(data = fcs_data)
clean_up_data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = clean_up_channels)

#************************************************
#gating on Residual+ ----
#************************************************
number_of_cells_before_residual_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_before_residual_gating)


#update lower_gate_percent, upper_gate_percent
residual_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Residual", lower_gate_percent = 15, upper_gate_percent = 20)
density_plots <- density_plot(data = clean_up_data, "Residual", plot_title = file_names, lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Residual", plot_title = file_names,  lower_gate = residual_gates$lower_gate, upper_gate = residual_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_gating <- cells_to_keep(data = clean_up_data, channel = "Residual", lower_gate = residual_gates$lower_gates, 
                                                  upper_gate = residual_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_gating)
number_of_cells_after_residual_gating <-  number_of_cells(data = fcs_data)
number_of_cells_after_residual_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_residual_gating/number_of_cells_after_beads_gating * 100 #percent remaining from bead gating


#************************************************
#gating on Center+ ----
#************************************************
number_of_cells_before_center_gating <- number_of_cells_after_residual_gating
random_cells_for_plotting <- random_cells(number_of_cells_before_center_gating)

#update lower_gate_percent, upper_gate_percent
center_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Center", lower_gate_percent = 10, upper_gate_percent = 10)
density_plots <- density_plot(data = clean_up_data, "Center", plot_title = file_names, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Center", plot_title = file_names, lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_center_gating <- cells_to_keep(data = clean_up_data, channel = "Center", lower_gate = center_gates$lower_gate, upper_gate = center_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_center_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_center_gating)
number_of_cells_after_center_gating <-  number_of_cells(data = fcs_data)
number_of_cells_after_center_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_center_gating/number_of_cells_after_residual_gating * 100 #percent remaining from residual gating

#************************************************
#gating on Offset+ ----
#************************************************
number_of_cells_before_offset_gating <-  number_of_cells_after_center_gating
random_cells_for_plotting <- random_cells(number_of_cells_before_offset_gating)


#update lower_gate_percent, upper_gate_percent
offset_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Offset", lower_gate_percent = 10, upper_gate_percent = 10)
density_plots <- density_plot(data = clean_up_data, "Offset", plot_title = file_names, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Offset", plot_title = file_names, lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_offset_gating <- cells_to_keep(data = clean_up_data, channel = "Offset",  lower_gate = offset_gates$lower_gate, upper_gate = offset_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_offset_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_offset_gating)
number_of_cells_after_offset_gating <- number_of_cells(data = fcs_data)
number_of_cells_after_offset_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_offset_gating/number_of_cells_after_center_gating * 100 #percent remaining from  center gating

#************************************************
#gating on Width+ ----
#************************************************
number_of_cells_before_width_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_before_width_gating)



#update lower_gate_percent, upper_gate_percent
width_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Width", lower_gate_percent = 15, upper_gate_percent = 15)

density_plots <- density_plot(data = clean_up_data, "Width", plot_title = file_names, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Width", plot_title = file_names, lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

width_gates$lower_gate
#lower_gate failed for three datasets. Put it manually to 2.8 for all datasets.
width_gates$lower_gate <- 2.8
cells_to_keep_after_width_gating <- cells_to_keep(data = clean_up_data, channel = "Width", lower_gate = width_gates$lower_gate, upper_gate = width_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_width_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_width_gating)
number_of_cells_after_width_gating <- number_of_cells(data = fcs_data)
number_of_cells_after_width_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_width_gating/number_of_cells_after_offset_gating * 100 #percent remaining from offset gating

#************************************************
#gating on Event+ ----
#************************************************
number_of_cells_before_event_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_before_event_gating)

#update lower_gate_percent, upper_gate_percent
EventGates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Event_length", lower_gate_percent = 15, upper_gate_percent = 10)
density_plots <- density_plot(data = clean_up_data, "Event_length", plot_title = file_names, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data =  clean_up_data, random_cells = random_cells_for_plotting, channel = "Event_length", plot_title = file_names, lower_gate = EventGates$lower_gate, upper_gate = EventGates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_event_gating <- cells_to_keep(data = clean_up_data, channel = "Event_length",  upper_gate = EventGates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_event_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_event_gating)
number_of_cells_after_event_gating <-  number_of_cells(data = fcs_data)
number_of_cells_after_event_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_event_gating/number_of_cells_after_width_gating * 100 #percent remaining from width gating


#************************************************
#gating on Live Dead----  Cis, 
#************************************************
number_of_cells_before_cis_gating <-  number_of_cells(data = fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_before_cis_gating)

#NB Her har jeg overstyrt algoritmen og sier at jeg vil ha minste negativ gate på 2.
#update lower_gate_percent, upper_gate_percent
cis_gates <- findNegativeGateFCS(exprFCSdata = clean_up_data, channel = "Pt195Di", lower_gate_percent = 15, upper_gate_percent = 10, minupper_gate = 2)
density_plots <- density_plot(data = clean_up_data, "Pt195Di", plot_title = file_names, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Pt195Di", plot_title = file_names, lower_gate = cis_gates$lower_gate, upper_gate = cis_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_cis_gating <- cells_to_keep(data = clean_up_data, channel = "Pt195Di",  lower_gate = cis_gates$lower_gate,  upper_gate = cis_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_cis_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_cis_gating)
number_of_cells_after_cis_gating <-  number_of_cells(data = fcs_data)
number_of_cells_after_cis_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_cis_gating/number_of_cells_after_event_gating * 100 #percent remaining from event gating


#************************************************
#gating on DNA1, Ir191Di+ ----
#************************************************
number_of_cells_before_Ir191Di_gating <-  number_of_cells(fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_before_Ir191Di_gating)

#update lower_gate_percent, upper_gate_percent
Ir191di_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Ir191Di", lower_gate_percent = 15, upper_gate_percent = 15)

density_plots <- density_plot(data = clean_up_data, "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Ir191Di", plot_title = file_names, lower_gate = Ir191di_gates$lower_gate, upper_gate = Ir191di_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_Ir191Di_gating <- cells_to_keep(data = clean_up_data, channel = "Ir191Di",  lower_gate = Ir191di_gates$lower_gate,  upper_gate = Ir191di_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_Ir191Di_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_Ir191Di_gating)
number_of_cells_after_Ir191Di_gating <-  number_of_cells(data = fcs_data)
number_of_cells_after_Ir191Di_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_Ir191Di_gating/number_of_cells_after_cis_gating * 100 #percent remaining from cis gating



#************************************************
#gating on DNA1, Ir193Di+ ----
#************************************************
number_of_cells_before_Ir193Di_gating <-  number_of_cells(fcs_data)
random_cells_for_plotting <- random_cells(number_of_cells_before_Ir193Di_gating)

#update lower_gate_percent, upper_gate_percent
Ir193di_gates <- find_gaussian_gates_second_top(data = clean_up_data, channel = "Ir193Di", lower_gate_percent = 15, upper_gate_percent = 15)

density_plots <- density_plot(data = clean_up_data, "Ir193Di", plot_title = file_names, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = clean_up_data, random_cells = random_cells_for_plotting, channel = "Ir193Di", plot_title = file_names, lower_gate = Ir193di_gates$lower_gate, upper_gate = Ir193di_gates$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_Ir193Di_gating <- cells_to_keep(data = clean_up_data, channel = "Ir193Di",  lower_gate = Ir193di_gates$lower_gate,  upper_gate = Ir193di_gates$upper_gate)
fcs_data <- update_data_based_on_cells_to_keep(data = fcs_data, kept_cells = cells_to_keep_after_Ir193Di_gating)
clean_up_data <- update_data_based_on_cells_to_keep(data = clean_up_data, kept_cells = cells_to_keep_after_Ir193Di_gating)
number_of_cells_after_Ir193Di_gating <-  number_of_cells(data = fcs_data)
number_of_cells_after_Ir193Di_gating/number_of_cells_raw_data * 100 #percent remaining from total
number_of_cells_after_Ir193Di_gating/number_of_cells_after_cis_gating * 100 #percent remaining from cis gating

#************************************************
#save fcs_data ----
#************************************************

write.flowSet(fcs_data, outdir = outDataPath, filename = filesToOpen)

stopp_time <- Sys.time()

stopp_time - start_time

