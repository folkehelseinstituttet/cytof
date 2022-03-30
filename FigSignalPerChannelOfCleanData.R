data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "clean data")
scriptPath <- fs::path("H:", "git", "cytof")
outFigPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "FigSignalPerChannelOfCleanData")

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
                               plot_list[[66]],  ncol = 11, nrow = 6)
  return(g)
}


# read all files in data_path into one dataset fcs_data
fcs_data_with_info <- read_some_data_from_folder(data_path, file_number = filenumbers)
fcs_data <- fcs_data_with_info$fcs_data
file_names <- factor(fcs_data_with_info$file_names, levels = fcs_data_with_info$file_names)
rm(fcs_data_with_info)

# get the parameters of fcs_data and store them in params.
params <- get_params_fcs_data(fcs_data[[1]])

params$kanal <- c("Time", NA, "CD45", NA, NA, NA, "CD57", NA, NA, "CD19", "CD5", "CD8", "HLADR", "CD3",
                  NA, NA, NA, NA, NA, NA, "CCR6", "KLRG1", "CD127", "CD15", "CD4", "IgD", "CD11c", "CD16",
                  "CD25", "CD134_OX40", "CD123", "TCRgd", "CCR4", "TIGIT", "CD45RA", "CXCR3", "CD27", "IgG", "CD28",
                  "CD160", "CD95", "CRTH2", "CD161", "CD85j", "TCRVa7.2", "CCR7", "ICOS", "NKG2A", "CD169", "CXCR5", 
                  "CD38", "CD141",  "PD-1", "CD14", "CD56",   rep(NA, 7), "CD11b", rep(NA,4))

kanaler <- c("CD3","CD4", "CD5", "CD8", "CD19", "CD45", "CD57", "CD56", "CCR4",
             "KLRG1", "CD127", "CD15", "IgD", "CD11c", "CD16", "CD25", 
             "CD134_OX40", "CD123", "TCRgd", "TIGIT", "CD45RA", "CXCR3", "CD27",
             "IgG", "CD28", "CD160", "CD85j", "TCRVa7.2", "CD161", "CRTH2", "CD95",
             "CCR7", "ICOS", "NKG2A", "CD169", "CXCR5", "CD38", "CD141", "PD-1", "CD14", 
             "CD56", "CD11b", "CCR6", "HLADR")

kanalname <- params$name[which(!is.na(params$kanal))] #må sjekke at vi får riktig kanal.

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanalname) 


number_of_events_data <-  number_of_events(data = data, file_names = file_names)
random_events_for_plotting <- random_events(number_of_events_data, n = 10000)
file_names_in_plot <- 1:length(data)

for(i in 2:nrow(data[[1]])){
  navn <- params$kanal[grep(colnames(data[[1]])[i], params$name)]
  time_signal_plots <- time_signal_plot(data = data, random_events = random_events_for_plotting, 
                                        channel = colnames(data[[1]])[i],  plot_title = file_names_in_plot)
  tiff(fs::path(outFigPath, paste0("Signal_clean_data_", navn,".tiff")), width = 1800, height = 1200)
  print(plotSignal(plot_list = time_signal_plots))
  dev.off()
}



