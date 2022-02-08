

# Example
# 
# 
ptm <- proc.time()
# path to where the fcs files are stored F:\Forskningsprosjekter\PDB 2794 - Immune responses aga_\Forskningsfiler\JOBO\CyTOF\Datafiles\Panel 1 all files
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Gating", "Gating fra R_FINAL")
resPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS")
out_result <- fs::path(resPath, "Resultat_Panel_1")
outDataPath <- fs::path(out_result, "Data")
outFigPath <- fs::path(out_result, "FigAlle")
outFigPathSignal <- fs::path(out_result, "FigSignalSignal")

scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "gating_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))

fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))
fcs_files <- fcs_files[grep(".fcs", fcs_files)]
# fcs_files_sever <- fcs_files[grep("_FINAL/S_", fcs_files)] 
# fcs_files_sever_T1 <- fcs_files_sever[grep("T1", fcs_files_sever)] 
# fcs_files_moderat <- fcs_files[grep("_FINAL/M_", fcs_files)] 
# fcs_files_moderat_T1 <- fcs_files_moderat[grep("T1", fcs_files_moderat)] 
# fcs_files_ref <- fcs_files[grep("_FINAL/Ref1_", fcs_files)] 
# 
# fcs_files <- c(fcs_files_sever_T1, fcs_files_moderat_T1, fcs_files_ref)

files_to_open <- basename(fcs_files)
setwd(dirname(fcs_files[1]))
file_names <- gsub(".fcs", "", files_to_open)
# Read the files into a flowset


n_files <- length(file_names)
filenumber <- 1:n_files


result <- list()


kanaler <- c("CD3","CD4", "CD5", "CD8", "CD19", "CD45", "CD57", "CD56", "CCR4",
             "KLRG1", "CD127", "CD15", "IgD", "CD11c", "CD16", "CD25", 
             "CD134_OX40", "CD123", "TCRgd", "TIGIT", "CD45RA", "CXCR3", "CD27",
             "IgG", "CD28", "CD160", "CD85j", "TCRVa7.2", "CD161", "CRTH2", "CD95",
             "CCR7", "ICOS", "NKG2A", "CD169", "CXCR5", "CD38", "CD141", "PD-1", "CD14", 
             "CD56", "CD11b", "CCR6", "HLADR")


filene <- 1:n_files
# read all files in data_path into one dataset fcs_data
fcs_data_with_info <- read_specific_data_from_folder(data_path = data_path, files_to_open = files_to_open[i])
fcs_data <- fcs_data_with_info$fcs_data
file_names <- factor(fcs_data_with_info$file_names, levels = fcs_data_with_info$file_names)
rm(fcs_data_with_info)
for(j in filene){
  mat <-  as.data.frame(matrix(NA, ncol = length(kanaler), nrow =  nrow(fcs_data[[j]])))
  colnames(mat) <- kanaler
  result[[j]] <- mat
}

# get the parameters of fcs_data and store them in params.
params <- get_params_fcs_data(fcs_data[[1]])


#***************************************************
#pos/neg CD3 ---- OK
#***************************************************
CD45 <- params$name[grep("CD45", params$desc)][1] 
x <- "CD3"
params$desc[grep(x, params$desc)][1] #må sjekke at vi får riktig kanal.

kanal <- params$name[grep(x, params$desc)][1] #må sjekke at vi får riktig kanal.
kanal
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45)) 

#split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 10, upper_gate_percent = 0.001)
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


#signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = "CD3", plot_title = file_names)
#
# for(i in 1:length(signal$plotList)){
#   tiff(fs::path(outFigPathSignal, paste0("fig_", x, "_gating", i, ".tiff")), height = 600, width = 600)
#     print(signal$plotList[[i]])
#   dev.off()
# }

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA
#***************************************************
#pos/neg CD4 ---- OK
#***************************************************
x <- "CD4"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][2] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD5 ---- OK  NB FHI084?
#***************************************************
x <- "CD5"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][2] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA


#***************************************************
#pos/neg CD8 ----  kun pos neg
#***************************************************
x <- "CD8"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 1] <- mean(split$lower_gates[split$lower_gates > 1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



split <- NA
#***************************************************
#pos/neg CD19 ---- kun pos neg
#***************************************************
x <- "CD19"
y <- "CD45"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
params$desc[grep(y, params$desc)]#må sjekke at vi får riktig kanal.
kanaly <- params$name[grep(y, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, kanaly)) 
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 4] <- mean(split$lower_gates[split$lower_gates > 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



split <- NA
#***************************************************
#pos/neg CD45 ----  OK Velger 01
#***************************************************
x <- "CD45"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA
#***************************************************
#pos/neg CD57 ----
#***************************************************
#Johanna

x <- "CD57"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])
split$low_high_splits[split$low_high_splits > 5] <- mean(split$low_high_splits[split$low_high_splits < 5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()



# 
# split <- find_split_neg_low_high(data = data, channel = kanal, negProp = 0.5, minLowHig = 2)
# split$low_high_splits[split$low_high_splits > 5] <- mean(split$low_high_splits[split$low_high_splits < 5])
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate =  split$low_high_splits, main_title = x)
# density_plots # to see the plots
# 
# # density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate =  split$low_high_splits)
# density_plots_pos # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")))
print(density_plots)
dev.off()

for(j in filene){
#  result[[j]][,x] <- data[[j]][, kanal] > split$neg_splits[j]
#  result[[j]][,x][data[[j]][, kanal] > split$low_high_splits[j]] <- 2

   result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA
#***************************************************
#pos/neg CD56 ---- kun pos neg
#***************************************************
x <- "CD56"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()



for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA
#***************************************************
#pos/neg CCR4 ---- OK
#***************************************************
x <- "CCR4"


params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


split <- NA
#***************************************************
#pos/neg KLRG1 ---- OK!
#***************************************************
x <- "KLRG1"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA
#***************************************************
#pos/neg CD127 ---- OK h?y lav neg
#***************************************************

x <- "CD127"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


split <- NA


#***************************************************
#pos/neg CD15 ---- OK
#***************************************************
x <- "CD15"
#NBNBNB
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA
#***************************************************
#pos/neg IgD ---- OK
#***************************************************
#NBNBNB
x <- "IgD"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_IgD_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_IgD_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"IgD"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"IgD"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD11c ---- OK
#***************************************************
#NBNBNBNB  veldig stor og rar spredning. Vil vi gjøre noe med denne?
x <- "CD11c"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD11c_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD11c_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD11c"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD11c"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD16 ---- OK
#***************************************************
x <- "CD16"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD16_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD16_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD16"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD16"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD25 ---- her tror jeg at vi m? sette en nedre grense p? 1
#***************************************************
#NBNBNB
#Johanna se svar fast nedre grense 1

x <- "CD25"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 3)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD25c_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD25c_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD25c"] <- data[[j]][, kanal] > split$neg_splits[j]
#   #  result[[j]][,"CD25c"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD134_OX40 ---- OK 
#***************************************************
#NBNBNB
x <- "CD134_OX40"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD134_OX40_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD134_OX40_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD134_OX40"] <- data[[j]][, kanal] > split$neg_splits[j]
#   #  result[[j]][,"CD134_OX40"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD123---- her m? vi sette fiks gating fra 5 ned til 0.1?
#***************************************************
#NBNBNB
#Johanna tror p? fiks gate her

x <- "CD123"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4])
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD123_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD123_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD123"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD123"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg TCRgd---- her ?nsker vi den lille toppen til h?yre, kanskje m? gate mot CD45?
#***************************************************
#NBNBNB

x <- "TCRgd"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_TCRgd_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_TCRgd_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"TCRgd"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"TCRgd"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg TIGIT---- her kan vi sette en nedre cut off p? 1, alt annet er positivt
#***************************************************
#NBNBNB
x <- "TIGIT"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 1.8] <- 1.8 #Johanna

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 1.8] <- 1.8 #Johanna
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_TIGIT_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_TIGIT_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"TIGIT"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"TIGIT"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD45RA----OK  high, low, negativ
#***************************************************
#NBNBNB
x <- "CD45RA"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 6] <- mean(split$low_high_splits[split$low_high_splits < 6]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD45RA_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD45RA_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD45RA"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD45RA"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CXCR3---- OK
#***************************************************
x <- "CXCR3"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 3] <- mean(split$low_high_splits[split$low_high_splits < 3]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CXCR3_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CXCR3_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CXCR3"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CXCR3"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA
#***************************************************
#pos/neg CD27---- OK
#***************************************************
x <- "CD27"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 6] <- mean(split$low_high_splits[split$low_high_splits < 6]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD27_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD27_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD27"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD27"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg IgG---- OK
#***************************************************
x <- "IgG"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_IgG_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_IgG_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"IgG"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"IgG"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD28---- OK
#***************************************************
x <- "CD28"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 6] <- mean(split$low_high_splits[split$low_high_splits < 6]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD28_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD28_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD28"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD28"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD160---- OK
#***************************************************
x <- "CD160"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD160_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD160_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD160"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD160"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD85j---- OK 
#***************************************************
x <- "CD85j"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# 
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD85j_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD85j_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD85j"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD85j"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg TCRVa7.2---- OK
#***************************************************
x <- "TCRVa7.2"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# #split$low_high_splits[split$low_high_splits > 6] <- mean(split$low_high_splits[split$low_high_splits < 6]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_TCRVa7.2_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_TCRVa7.2_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"TCRVa7.2"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"TCRVa7.2"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD161---- OK
#***************************************************
x <- "CD161"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}
# 
# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD161_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD161_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD161"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD161"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA

#***************************************************
#pos/neg CRTH2---- OK
#***************************************************
x <- "CRTH2"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates < 1.5] <- mean(split$lower_gates[split$lower_gates > 1.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits < 1.5] <- mean(split$low_high_splits[split$low_high_splits > 1.5]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CRTH2_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CRTH2_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CRTH2"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CRTH2"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA

#***************************************************
#pos/neg CD95---- OK
#***************************************************
x <- "CD95"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD95_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD95_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD95"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD95"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA

#***************************************************
#pos/neg CCR7---- OK
#***************************************************
x <- "CCR7"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CCR7_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CCR7_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CCR7"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CCR7"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA

#***************************************************
#pos/neg ICOS---- OK
#***************************************************
x <- "ICOS"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 2.5] <- mean(split$low_high_splits[split$low_high_splits < 2.5]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_ICOS_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_ICOS_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"ICOS"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"ICOS"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA

#***************************************************
#pos/neg NKG2A---- OK
#***************************************************
x <- "NKG2A"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_NKG2A_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_NKG2A_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"NKG2A"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"NKG2A"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA


#***************************************************
#pos/neg CD169---- OK?
#***************************************************
x <- "CD169"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD169_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD169_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD169"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD169"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA

#***************************************************
#pos/neg CD38---- OK
#***************************************************
x <- "CD38"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# split$low_high_splits <- rep(2, length(split$low_high_splits))
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD38_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD38_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD38"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD38"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA

#***************************************************
#pos/neg CD141---- OK
#***************************************************
x <- "CD141"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 2.5] <- mean(split$low_high_splits[split$low_high_splits < 2.5]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD141_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD141_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD141"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD141"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA

#***************************************************
#pos/neg PD-1----OK
#***************************************************
x <- "PD-1"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_PD-1_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_PD-1_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"PD-1"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"PD-1"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA


#***************************************************
#pos/neg CD14---- OK
#***************************************************
x <- "CD14"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4] <- mean(split$low_high_splits[split$low_high_splits < 4]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD14_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD14_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD14"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD14"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD56---- OK high low
#***************************************************
x <- "CD56"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
#split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD56_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD56_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD56"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD56"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD11b---- her velger vi den lille toppen helt til h?yre ca 4, resten negativt
#***************************************************
x <- "CD11b"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 4 | split$low_high_splits < 2] <- mean(split$low_high_splits[split$low_high_splits < 4 & split$low_high_splits > 2 ]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD11b_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD11b_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD11b"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD11b"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CCR6---- dette er ok, high, low, negativ
#***************************************************
x <- "CCR6"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 2.5] <- mean(split$low_high_splits[split$low_high_splits < 2.5]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CCR6_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CCR6_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CCR6"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CCR6"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg CD 160 OK positiv=topp til h?yre
#***************************************************
x <- "CD160"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 3.5] <- mean(split$lower_gates[split$lower_gates < 3.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}
# 
# 
# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 3.5] <- mean(split$low_high_splits[split$low_high_splits < 3.5]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD160_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CD160_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CD160"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CD160"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }

split <- NA
#***************************************************
#pos/neg HLADR---- OK positiv/negativ topp til h?yre
#***************************************************
x <- "HLADR"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 3] <- mean(split$low_high_splits[split$low_high_splits < 3]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_HLADR_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_HLADR_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"HLADR"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"HLADR"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
split <- NA

#***************************************************
#pos/neg ICOS OK
#***************************************************
x <- "ICOS"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



split <- NA



#***************************************************
#pos/neg CXCR5 OK
#***************************************************
x <- "CXCR5"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[split$lower_gates > 2] <- mean(split$lower_gates[split$lower_gates < 2])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


# 
# split <- find_split_neg_low_high(data = data, channel = kanal, neg = 0.05, minLowHig = 1)
# split$low_high_splits[split$low_high_splits > 2] <- mean(split$low_high_splits[split$low_high_splits < 2]) #Johanna
# 
# density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots # to see the plots
# 
# 
# density_plots_pos <- density_plot_without_neg(data = data, channel = kanal, plot_title = file_names, lower_gate = split$neg_splits, upper_gate = split$low_high_splits)
# density_plots_pos # to see the plots
# 
# 
# tiff(fs::path(outFigPath, paste0("fig1_CXCR5_gating", ".tiff")))
# print(density_plots)
# dev.off()
# 
# tiff(fs::path(outFigPath, paste0("fig1_CXCR5_gating_pos", ".tiff")))
# print(density_plots_pos)
# dev.off()
# 
# if(length(split$neg_splits) == 1){
#   split$neg_splits <- rep(split$neg_splits, length(data))
# }
# 
# for(j in filene){
#   result[[j]][,"CXCR5"] <- data[[j]][, kanal] > split$neg_splits[j]
#   result[[j]][,"CXCR5"][data[[j]][, kanal] > split$low_high_splits[j]] <- 2
# }
# 


saveRDS(result, fs::path(outDataPath, "posNeg.rds"))
saveRDS(file_names, fs::path(outDataPath, "posNegFilnavn.rds"))

proc.time() - ptm
