

# Example
#  
# 
ptm <- proc.time()
# path to where the fcs files are stored F:\Forskningsprosjekter\PDB 2794 - Immune responses aga_\Forskningsfiler\JOBO\CyTOF\Datafiles\Panel 1 all files
#data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Gating", "Gating fra R Panel 2")
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "clean data")
resPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg")
outDataPath <- fs::path(resPath, "Data")
outFigPath <- fs::path(resPath, "FigDensity")
outFigPathSignal <- fs::path(resPath, "FigSignalSignal")

scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "gating_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))


#denne funksjonen må tilpasset riktig antall figurer. 

plotTiff <- function(signal, filnavn){
  tiff(fs::path(outFigPathSignal, filnavn), height = 1200, width = 1800)
  gridExtra::grid.arrange(signal$plotList[[1]], signal$plotList[[2]], signal$plotList[[3]], signal$plotList[[4]], signal$plotList[[5]],
                          signal$plotList[[6]],  signal$plotList[[7]], signal$plotList[[8]], signal$plotList[[9]], signal$plotList[[10]],
                          signal$plotList[[11]], signal$plotList[[12]], signal$plotList[[13]], signal$plotList[[14]], signal$plotList[[15]],
                          signal$plotList[[16]],  signal$plotList[[17]], signal$plotList[[18]], signal$plotList[[19]], signal$plotList[[20]],
                          signal$plotList[[21]], signal$plotList[[22]], signal$plotList[[23]], signal$plotList[[24]], signal$plotList[[25]],
                          signal$plotList[[26]],  signal$plotList[[27]], signal$plotList[[28]], signal$plotList[[29]], signal$plotList[[30]],
                          signal$plotList[[31]], signal$plotList[[32]], signal$plotList[[33]], signal$plotList[[34]], signal$plotList[[35]],
                          signal$plotList[[36]],  signal$plotList[[37]], signal$plotList[[38]], signal$plotList[[39]], signal$plotList[[40]],
                          signal$plotList[[41]], signal$plotList[[42]], signal$plotList[[43]], signal$plotList[[44]], signal$plotList[[45]],
                          signal$plotList[[46]],  signal$plotList[[47]], signal$plotList[[48]], signal$plotList[[49]], signal$plotList[[50]],
                          signal$plotList[[51]], signal$plotList[[52]], signal$plotList[[53]], signal$plotList[[54]], signal$plotList[[55]],
                          signal$plotList[[56]],  signal$plotList[[57]], signal$plotList[[58]], signal$plotList[[59]], signal$plotList[[60]],
                          signal$plotList[[61]], signal$plotList[[62]], signal$plotList[[63]], signal$plotList[[64]], signal$plotList[[65]],
                          signal$plotList[[66]],  ncol = 11, nrow = 6)
  #print(signal$plotList[[i]])
  dev.off()
}


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






filene <- 1:n_files
# read all files in data_path into one dataset fcs_data
fcs_data_with_info <- read_specific_data_from_folder(data_path = data_path, files_to_open = files_to_open[filene])
fcs_data <- fcs_data_with_info$fcs_data
file_names <- factor(fcs_data_with_info$file_names, levels = fcs_data_with_info$file_names)
rm(fcs_data_with_info)
# get the parameters of fcs_data and store them in params.
params <- get_params_fcs_data(fcs_data[[1]])

kanaler <- c("CD45", "CD57",  "CD19", "CD8", "HLADR", "CD3", "CD4", "TCRgd", "CXCR5",
             "CD45RA",  "CD27", "CD28", "CCR7", "CD25", "CD38", "PD1", "CD14", "CD56",
              "CD16",  #de før dette var også på Plate 1.
             "CD107a", "CD44", "CD223", "IL-1b", "CD127", "IL-2", "TNFa", "TIM-3",
             "PD-L1", "IL-12p70", "MIP-1b", "CD137", "CD272", "IL-6", "GM-CSF", 
             "IL-17A", "FoxP3", "CD33", "Perforin", "IFNg", "IL-10",  "CD154", "CTLA-4",
             "GranzymeB", "PD-L2") 
             


result <- list()


for(j in filene){
  mat <-  as.data.frame(matrix(NA, ncol = length(kanaler), nrow =  nrow(fcs_data[[j]])))
  colnames(mat) <- kanaler
  result[[j]] <- mat
}


#***************************************************
#pos/neg CD45 ---- 14 feb: fast grense på 0.5
#***************************************************
CD107a <- params$name[grep( "CD107a", params$desc)][1] 

x <- "CD45"
print(x)

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates <-rep(0.5, length(split$lower_gates))

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD107a))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD107a, channel2 = kanal, ylow = split$lower_gates, xname = "CD107a", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)

plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))


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
#pos/neg CD57 ----  14 feb fast på 2.5
#***************************************************
CD45 <- params$name[grep("CD45", params$desc)][1] 

x <- "CD57"
print(x)

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 5, upper_gate_percent = 0.001, minimum = 1) #NB var minimum 3.75 for panel 1
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates <- rep(2.5, length(split$lower_gates))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA


#***************************************************
#pos/neg CD19 ---- 
#***************************************************
x <- "CD19"
print(x)


params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates < 4] <- mean(split$lower_gates[split$lower_gates > 4])


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


split <- NA

#***************************************************
#pos/neg CD8 ----  
#***************************************************
x <- "CD8"
print(x)

kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 2)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA




saveRDS(result, fs::path(outDataPath, "posNeg.rds"))
saveRDS(file_names, fs::path(outDataPath, "posNegFilnavn.rds"))

proc.time() - ptm


#***************************************************
#pos/neg HLADR ----  
#***************************************************
x <- "HLADR"
print(x)

kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))


splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates[is.na(splitHigh$lower_gates)] <- mean(splitHigh$lower_gates[!is.na(splitHigh$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA


#***************************************************
#pos/neg CD3 ----  
#***************************************************
x <- "CD3"
print(x)

kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 0.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA



#***************************************************
#pos/neg CD4 ----  
#***************************************************
x <- "CD4"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][3] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 2] <- mean(split$lower_gates[split$lower_gates < 2])
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)

plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA


#***************************************************
#pos/neg CD25 ----  
#***************************************************
x <- "CD25"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))

 
splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates[is.na(splitHigh$lower_gates)] <- mean(splitHigh$lower_gates[!is.na(splitHigh$lower_gates)])
splitHigh$lower_gates[splitHigh$lower_gates < 0.1] <- mean(splitHigh$lower_gates[splitHigh$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA


#***************************************************
#pos/neg TCRgd ----  
#***************************************************
x <- "TCRgd"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
# split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CXCR5 ----  14 feb Velger fast på 2
#***************************************************
x <- "CXCR5"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
# split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])
split$lower_gates <- rep(2, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA


#***************************************************
#pos/neg CD45RA ----  
#***************************************************
x <- "CD45RA"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
# split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])



density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA



#***************************************************
#pos/neg CD27 ----  
#***************************************************
x <- "CD27"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][3] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))



splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates[is.na(splitHigh$lower_gates)] <- mean(splitHigh$lower_gates[!is.na(splitHigh$lower_gates)])



density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
# split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])



kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA


#***************************************************
#pos/neg CD28 ----  
#***************************************************
x <- "CD28"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = NA)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
# split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CCR7 ----  14 feb inkl >3.5
#***************************************************
x <- "CCR7"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 3.5] <- mean(split$lower_gates[split$lower_gates < 3.5])



density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA



#***************************************************
#pos/neg CD38 ----  
#***************************************************
x <- "CD38"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))


splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates <- rep(0.8, length(splitHigh$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA


#***************************************************
#pos/neg PD1 ----  15 feb fix 1
#***************************************************
x <- "PD1"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
#split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD14 ----  
#***************************************************
x <- "CD14"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))


splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates <- rep(1.5, length(splitHigh$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA


#***************************************************
#pos/neg CD56 ----  
#***************************************************
x <- "CD56"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))


splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates[is.na(splitHigh$lower_gates)] <- mean(splitHigh$lower_gates[!is.na(splitHigh$lower_gates)])


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA


#***************************************************
#pos/neg CD16 ----  
#***************************************************
x <- "CD16"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))


splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates <- rep(1.5, length(splitHigh$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA

#***************************************************
#pos/neg CD107a ----  
#***************************************************
x <- "CD107a"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
#split$lower_gates <- rep(1.5, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD44 ----  
#***************************************************
x <- "CD44"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
#split$lower_gates <- rep(1.5, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD223 ----  
#***************************************************
x <- "CD223"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1.5, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg IL-1b ----  #14 feb fast grense på 3
#***************************************************
x <- "IL-1b"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(3, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

kanal_max <- 6

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD127 ----  #legg inn upper ... 15 feb
#***************************************************
x <- "CD127"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 0.25)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 1] <- mean(split$lower_gates[split$lower_gates < 1])
split$upper_gates[is.na(split$lupper_gates)] <- mean(split$upper_gates[!is.na(split$upper_gates)])
split$upper_gates[split$upper_gates > 3] <- mean(split$upper_gates[split$upper_gates < 3])
split$lower_gates <- rep(0.3, length(split$lower_gates))
split$upper_gates <- rep(2.5, length(split$upper_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, upper_gate = split$upper_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates,  yhigh = split$upper_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
  result[[j]][data[[j]][, kanal] > split$upper_gates[j],x] <- 2
}

split <- NA

#***************************************************
#pos/neg  IL-2 ----  14 feb, fast grense på 2.5
#***************************************************
x <- "IL-2"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 1] <- mean(split$lower_gates[split$lower_gates < 1])
split$lower_gates <- rep(2.5, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

kanal_max <- 6
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg TNFa ----  15 feb fix 2
#***************************************************
x <- "TNFa"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 1] <- mean(split$lower_gates[split$lower_gates < 1])
split$lower_gates <- rep(2, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}


signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg TIM-3 ----  15 feb fix 0.25
#***************************************************
x <- "TIM-3"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(0.25, length(split$lower_gates))


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}


kanal_max <- 2

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg PD-L1 ----  15 feb fix 1
#***************************************************
x <- "PD-L1"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1, length(split$lower_gates))

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
kanal_max <- 4

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg IL-12p70 ----  15 feb fix grense på 2.4
#***************************************************
x <- "IL-12p70"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(2.4, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg MIP-1b ----  
#***************************************************
x <- "MIP-1b"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.8)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
#split$lower_gates <- rep(1.5, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD137 ----  fast grense på 0.5
#***************************************************
x <- "CD137"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(0.5, length(split$lower_gates))



kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

kanal_max <- 2

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD272 ----  14 feb fast grense på 1
#***************************************************
x <- "CD272"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1, length(split$lower_gates))



kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
kanal_max <- 5

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg IL-6 ----  # 14 feb fast grense på 1.5
#***************************************************
x <- "IL-6"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1.5, length(split$lower_gates))



kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
kanal_max <- 4

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg GM-CSF ----  14 feb, fast grense på 1
#***************************************************
x <- "GM-CSF"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 0.9)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg IL-17A ----  15 feb fix på 1-6
#***************************************************
x <- "IL-17A"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 0.9)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1.7, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg FoxP3 ----  14 feb fast grense på 0.3
#***************************************************
x <- "FoxP3"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))

splitLow <- find_gaussian_gates_first_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 15)

splitHigh <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001, minimum = 1)
splitHigh$lower_gates <- rep(0.3, length(splitHigh$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = splitLow$upper_gates, upper_gate = splitHigh$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = splitLow$upper_gates, yhigh = splitHigh$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > splitLow$upper_gates[j]
  result[[j]][data[[j]][, kanal] > splitHigh$lower_gates[j],x] <- 2
}



splitHigh <- NA
splitLow <- NA

#***************************************************
#pos/neg CD33 ----  fast grense 1
#***************************************************
x <- "CD33"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(1, length(split$lower_gates))

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

kanal_max <- 4

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg Perforin ----  15 feb fix 2.3
#***************************************************
x <- "Perforin"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.8)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(2.3, length(split$lower_gates))


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg IFNg ----  14 feb fast 2.5
#***************************************************
x <- "IFNg"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 2.1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(2.5, length(split$lower_gates))



kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg IL-10 ----  15 feb fix grense på 2
#***************************************************
x <- "IL-10"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates <- rep(2, length(split$lower_gates))


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
kanal_max <- 5

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CD154 ----  
#***************************************************
x <- "CD154"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
#split$lower_gates <- rep(1.5, length(split$lower_gates))


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg CTLA-4 ----  fast på 1.4
#***************************************************
x <- "CTLA-4"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])
split$lower_gates <- rep(1.4, length(split$lower_gates))



kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg GranzymeB ----  14 feb øker minimum til 3
#***************************************************
x <- "GranzymeB"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 3)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 5] <- mean(split$lower_gates[split$lower_gates < 5])
#split$lower_gates <- rep(1.5, length(split$lower_gates))


density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA

#***************************************************
#pos/neg PD-L2 ----  
#***************************************************
x <- "PD-L2"
print(x)

params$desc[grep(x, params$desc)]
kanal <- params$name[grep(x, params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)#, minimum = 1.4)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
#split$lower_gates <- rep(1.5, length(split$lower_gates))


kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
kanal_max <- 4

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x, xlim = c(0, kanal_max))
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}

split <- NA



saveRDS(result, fs::path(outDataPath, "posNeg.rds"))
saveRDS(file_names, fs::path(outDataPath, "posNegFilnavn.rds"))

proc.time() - ptm
