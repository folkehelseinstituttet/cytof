# Example
# 
# 

# path to where the fcs files are stored
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2770 - EXIMIOUS - Mapping E_", "CyTOF in vitro WP4", "GatedData")
outDataPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2770 - EXIMIOUS - Mapping E_", "CyTOF in vitro WP4", "GatedData", "Tceller")
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

####
# functional cells
####

params$desc[grepl("CD3", params$desc)]
VarNavn1 <- params$name[grepl("CD3", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD19", params$desc)]
VarNavn2 <- params$name[grepl("CD19", params$desc)][1]
params$desc[grepl("CD4", params$desc)]
VarNavn3 <- params$name[grepl("CD4", params$desc)][2]  #NB pass på å velge riktig....
params$desc[grepl("CCR7", params$desc)]
VarNavn4 <- params$name[grepl("CCR7", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD45RA", params$desc)]
VarNavn5 <- params$name[grepl("CD45RA", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD8", params$desc)]
VarNavn6 <- params$name[grepl("CD8", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD25", params$desc)]
VarNavn7 <- params$name[grepl("CD25", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD127", params$desc)]
VarNavn8 <- params$name[grepl("CD127", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD56", params$desc)]
VarNavn9 <- params$name[grepl("CD56", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD16", params$desc)]
VarNavn10 <- params$name[grepl("CD16", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD14", params$desc)]
VarNavn11 <- params$name[grepl("CD14", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("CD11c", params$desc)] #denne er ikke med
VarNavn12 <- params$name[grepl("CD11c", params$desc)][1]  #NB pass på å velge riktig....
params$desc[grepl("HLA-DR", params$desc)] #denne er ikke med
VarNavn13 <- params$name[grepl("HLA-DR", params$desc)][1]  #NB pass på å velge riktig....

#params$desc[grepl("CD123", params$desc)] #denne er ikke med
#params$desc[grepl("CD38", params$desc)] #denne er ikke med
#params$desc[grepl("CD27", params$desc)] #denne er ikke med
#params$desc[grepl("CD24", params$desc)] #denne er ikke med



arcsinhexprData <- arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(VarNavn1, VarNavn2, VarNavn3, 
                                                                                 VarNavn4, VarNavn5, VarNavn6, 
                                                                                 VarNavn7, VarNavn8, VarNavn9, 
                                                                                 VarNavn10, VarNavn11, VarNavn12,
                                                                                 VarNavn13))

for(i in 1:length(arcsinhexprData)){
  colnames(arcsinhexprData[[i]]) <- c("CD3", "CD19", "CD4", "CCR7", "CD45RA", "CD8", "CD25", "CD127", "CD56", "CD16", "CD14", "CD11c", "HLA-DR", "Time")
}



number_of_cells_after_clean_up <-  number_of_cells(arcsinhexprData)
random_cells_for_plotting <- random_cells(number_of_cells_after_clean_up)


#CD3+
CD3_pos_gate <- find_gaussian_gates_second_top(data = arcsinhexprData, channel = "CD3", lower_gate_percent = 20, upper_gate_percent = 20)
density_plots <- density_plot(data = arcsinhexprData, "CD3", plot_title = file_names, lower_gate = CD3_pos_gate$lower_gate, upper_gate = CD3_pos_gate$upper_gate)
#density_plots
time_signal_plots <- time_signal_plot(data = arcsinhexprData, random_cells = random_cells_for_plotting, channel = "CD3", plot_title = file_names,  lower_gate = CD3_pos_gate$lower_gate, upper_gate = CD3_pos_gate$upper_gate)
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_CD3pos_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD3", lower_gate = CD3_pos_gate$lower_gates, 
                                            upper_gate = CD3_pos_gate$upper_gate)
# for(i in 1:length(cells_to_keep_after_CD3pos_gating)){
#   print(table(cells_to_keep_after_CD3pos_gating[[i]]))
# }
CD3_neg_gate <- CD3_pos_gate$lower_gates
density_plots <- density_plot(data = arcsinhexprData, "CD3", plot_title = file_names, lower_gate = NA, upper_gate = CD3_neg_gate)
cells_to_keep_after_CD3neg_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD3", lower_gate = NA, 
                                                   upper_gate = CD3_neg_gate)

#density_plots


#CD19-
#plot to see distrubution of data before desiding on gating strategy 
density_plots <- density_plot(data = arcsinhexprData, "CD19", plot_title = file_names, lower_gate = NA, upper_gate = NA, xlim = c(0,2))

CD19_neg_gate <- find_gate_perc_height_upper_noise(data = arcsinhexprData, channel = "CD19",  upper_perc_height = 0.001)
density_plots <- density_plot(data = arcsinhexprData, "CD19", plot_title = file_names, lower_gate = NA, upper_gate = CD19_neg_gate, xlim = c(0,2))
#density_plots
time_signal_plots <- time_signal_plot(data = arcsinhexprData, random_cells = random_cells_for_plotting, channel = "CD19", 
                                      plot_title = file_names,  lower_gate = CD19_neg_gate, upper_gate = CD19_neg_gate, ylim = c(0,1))
#time_signal_plots # to see all plots
#time_signal_plots[1] # to see first plot

cells_to_keep_after_CD19neg_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD19", lower_gate = NA, 
                                                   upper_gate = CD19_neg_gate)
# for(i in 1:length(cells_to_keep_after_CD19neg_gating)){
#   print(table(cells_to_keep_after_CD19neg_gating[[i]]))
# }
# 
# 




density_plots <- density_plot(data = arcsinhexprData, channel = "CD19", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD4", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CCR7", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD45RA", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD8", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD25", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD127", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD56", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD16", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD14", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots
density_plots <- density_plot(data = arcsinhexprData, channel = "CD11c", plot_title = file_names, lower_gate = NA, upper_gate = NA)
density_plots




CD4gates <-  find_gaussian_gates_second_top(data = arcsinhexprData, channel = "CD4", lower_gate_percent = 0.2,upper_gate_percent = 0.2)
density_plots <- density_plot(data = arcsinhexprData, "CD4", plot_title = file_names, lower_gate = CD4gates$lower_gates, upper_gate = CD4gates$upper_gates)

CD4gates2 <-  find_split_first_second_top(data = arcsinhexprData, channel = "CD4")
density_plots <- density_plot(data = arcsinhexprData, "CD4", plot_title = file_names, lower_gate = CD4gates$lower_gates, upper_gate = CD4gates2)



CCR7gates <-  find_gate_gaussian_first_top(data = arcsinhexprData, channel = "CCR7", perc_included = 0.995)
density_plots <- density_plot(data = arcsinhexprData, "CCR7", plot_title = file_names, lower_gate = CCR7gates)
density_plots

CD8gates <-  find_gate_gaussian_first_top(data = arcsinhexprData, channel = "CD8", perc_included = 0.995)
density_plots <- density_plot(data = arcsinhexprData, "CD8", plot_title = file_names, lower_gate = CD8gates)
density_plots


CD19gates <-  find_gate_gaussian_first_top(data = arcsinhexprData, channel = "CD19", perc_included = 0.995)
density_plots <- density_plot(data = arcsinhexprData, "CD19", plot_title = file_names, lower_gate = CD19gates, xlim = c(0,2))
density_plots
cells_to_keep_after_CD19neg_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD19", lower_gate = NA, 
                                                    upper_gate = CD19gates)


CD45RAgates <-  find_gate_gaussian_first_top(data = arcsinhexprData, channel = "CD45RA", perc_included = 0.995)
density_plots <- density_plot(data = arcsinhexprData, "CD45RA", plot_title = file_names, lower_gate = CD45RAgates)
density_plots
cells_to_keep_after_CD45RAneg_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD45RA", lower_gate = NA, 
                                                    upper_gate = CD45RAgates)



CD127gates <-  find_gate_gaussian_first_top(data = arcsinhexprData, channel = "CD127", perc_included = 0.995)
density_plots <- density_plot(data = arcsinhexprData, "CD127", plot_title = file_names, lower_gate = CD127gates)
density_plots
cells_to_keep_after_CD127neg_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD127", lower_gate = NA, 
                                                    upper_gate = CD127gates)



CD11cgates <-  find_gate_gaussian_first_top(data = arcsinhexprData, channel = "CD11c", perc_included = 0.995)
density_plots <- density_plot(data = arcsinhexprData, "CD11c", plot_title = file_names, lower_gate = CD11cgates, xlim = c(0,2))
density_plots
cells_to_keep_after_CD11cneg_gating <- cells_to_keep(data = arcsinhexprData, channel = "CD11c", lower_gate = NA, 
                                                     upper_gate = CD11cgates)


CCR7gates2 <-  find_split_first_second_top(data = arcsinhexprData, channel = "CCR7")
density_plots <- density_plot(data = arcsinhexprData, "CCR7", plot_title = file_names, lower_gate = CCR7gates)



CD45RAgates <-  find_gaussian_gates_second_top(data = arcsinhexprData, channel = "CD45RA", lower_gate_percent = 0.2,upper_gate_percent = 0.2)
density_plots <- density_plot(data = arcsinhexprData, "CD45RA", plot_title = file_names, lower_gate = CD45RAgates$lower_gates, upper_gate = CD45RAgates$upper_gates)

CD45RAgates2 <-  find_split_first_second_top(data = arcsinhexprData, channel = "CD45RA")
density_plots <- density_plot(data = arcsinhexprData, "CD45RA", plot_title = file_names, lower_gate = CD45RAgates$lower_gates, upper_gate = CD45RAgates2)


CD45RAgates <-  find_gaussian_gates_first_top(data = arcsinhexprData, channel = "CD45RA", lower_gate_percent = 0.2,upper_gate_percent = 0.2)
density_plots <- density_plot(data = arcsinhexprData, "CD45RA", plot_title = file_names, lower_gate = CD45RAgates$lower_gates, upper_gate = CD45RAgates$upper_gates)


plots <- signal_signal_just_plot(data = arcsinhexprData, random_cells = random_cells_for_plotting, 
                                             channel1 = "CD4", channel2 = "CD45RA", 
                                             plot_title = NA, xlim = NA, ylim = NA)
  

