

# Example
#  
# 
ptm <- proc.time()
# path to where the fcs files are stored F:\Forskningsprosjekter\PDB 2794 - Immune responses aga_\Forskningsfiler\JOBO\CyTOF\Datafiles\Panel 1 all files
data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Gating", "Gating fra R_FINAL")
resPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS")
out_result <- fs::path(resPath, "Resultat_Panel_1")
outDataPath <- fs::path(out_result, "Data")
outFigPath <- fs::path(out_result, "FigAlle_9feb")
outFigPathSignal <- fs::path(out_result, "FigSignalSignal_9feb")

scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "gating_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))


#denne funksjonen må tilpasset riktig antall figurer. 

plotTiff <- function(signal, filnavn){
  tiff(fs::path(outFigPathSignal, filnavn), height = 1000, width = 1800)
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


result <- list()


kanaler <- c("CD3","CD4", "CD5", "CD8", "CD19", "CD45", "CD57", "CD56", "CCR4",
             "KLRG1", "CD127", "CD15", "IgD", "CD11c", "CD16", "CD25", 
             "CD134_OX40", "CD123", "TCRgd", "TIGIT", "CD45RA", "CXCR3", "CD27",
             "IgG", "CD28", "CD160", "CD85j", "TCRVa7.2", "CD161", "CRTH2", "CD95",
             "CCR7", "ICOS", "NKG2A", "CD169", "CXCR5", "CD38", "CD141", "PD-1", "CD14", 
             "CD56", "CD11b", "CCR6", "HLADR")


filene <- 1:n_files
# read all files in data_path into one dataset fcs_data
fcs_data_with_info <- read_specific_data_from_folder(data_path = data_path, files_to_open = files_to_open[filene])
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
#pos/neg CD3 ---- OK  8 feb ok
#***************************************************
CD45 <- params$name[grep("CD45", params$desc)][1] 
x <- "CD3"
params$desc[grep(x, params$desc)][1] #må sjekke at vi får riktig kanal.

kanal <- params$name[grep(x, params$desc)][1] #må sjekke at vi får riktig kanal.
kanal
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45)) 

#split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 10, upper_gate_percent = 0.001)
split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates < 1] <- mean(split$lower_gates[split$lower_gates > 1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD4 ---- OK
#***************************************************
x <- "CD4"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][2] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
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
#pos/neg CD5 ---- OK    8 feb ok.
#***************************************************
x <- "CD5"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][2] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

split <- NA


#***************************************************
#pos/neg CD8 ----  kun pos neg
#***************************************************
x <- "CD8"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates < 2.5] <- mean(split$lower_gates[split$lower_gates > 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))


split <- NA
#***************************************************
#pos/neg CD19 ---- kun pos neg
#***************************************************
x <- "CD19"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates < 4] <- mean(split$lower_gates[split$lower_gates > 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}




data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))




split <- NA
#***************************************************
#pos/neg CD45 ----  OK Velger 01
#***************************************************
x <- "CD45"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 2, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD57 ----  8 feb   SJEKK øker til 3.75 og minker til 5% på toppen
#***************************************************
#Johanna

x <- "CD57"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 5, upper_gate_percent = 0.001, minimum = 3.75)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))

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
#pos/neg CD56 ---- kun pos neg   8 feb ok
#***************************************************
x <- "CD56"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()



data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CCR4 ---- 8 feb ok 
#***************************************************
x <- "CCR4"


params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 2)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

split <- NA
#***************************************************
#pos/neg KLRG1 ---- OK!
#***************************************************
x <- "KLRG1"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}




data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

split <- NA
#***************************************************
#pos/neg CD127 ---- OK h?y lav neg
#***************************************************

x <- "CD127"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 0.8)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))



split <- NA


#***************************************************
#pos/neg CD15 ---- OK    8 feb. lavt signal overalt. Ta ut av datasettet? SJEKK. Har endret til 0.5
#***************************************************
x <- "CD15"
#NBNBNB
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 0.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg IgD ---- OK
#***************************************************
#NBNBNB
x <- "IgD"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}




data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))


split <- NA
#***************************************************
#pos/neg CD11c ---- OK  8 feb. bytt til minimum 3 sjekk! sjekk FHI004, FHI111???
#***************************************************
#NBNBNBNB  veldig stor og rar spredning. Vil vi gjøre noe med denne?
x <- "CD11c"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 3)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)] )
density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD16 ---- OK  8 feb- endr til minimum 1.5 SJEKK
#*************************************************** 
x <- "CD16"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD25 ---- her tror jeg at vi m? sette en nedre grense p? 1  8 feb endret minimum til 1. SJEKK om man trenger grense?
#***************************************************

x <- "CD25"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates < 0.1] <- mean(split$lower_gates[split$lower_gates > 0.1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD134_OX40 ---- OK 
#***************************************************
#NBNBNB
x <- "CD134_OX40"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD123---- her m? vi sette fiks gating fra 5 ned til 0.1?
#***************************************************

x <- "CD123"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 

split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots


tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg TCRgd---- her ?nsker vi den lille toppen til h?yre, kanskje m? gate mot CD45?
#***************************************************
#NBNBNB

x <- "TCRgd"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg TIGIT---- her kan vi sette en nedre cut off p? 1, alt annet er positivt
#***************************************************
#NBNBNB
x <- "TIGIT"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

split$lower_gates[split$lower_gates > 1.8] <- 1.8 #Johanna

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD45RA----OK  high, low, negativ  8 feb:  FHI103???? SJEKK, bytt til minimum 1.5?  TENK på klar bunn ved minimum??? tester også 5% istedenfor 15%
#***************************************************  #tenk på de med klar andre topp. 
#NBNBNB
x <- "CD45RA"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 5, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CXCR3---- OK
#***************************************************
x <- "CXCR3"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 20, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])
split$lower_gates[split$lower_gates < 1] <- mean(split$lower_gates[split$lower_gates > 1])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD27---- OK  8 feb beholder min grense på 1
#***************************************************
x <- "CD27"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg IgG---- OK
#***************************************************
x <- "IgG"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 


split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
#split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD28---- OK  8 feb ok
#***************************************************
x <- "CD28"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD160---- OK
#***************************************************
x <- "CD160"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD85j---- OK 
#***************************************************
x <- "CD85j"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg TCRVa7.2---- OK
#***************************************************
x <- "TCRVa7.2"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

#split$lower_gates[split$lower_gates > 6] <- mean(split$lower_gates[split$lower_gates < 6])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD161---- OK
#***************************************************
x <- "CD161"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CRTH2---- OK
#***************************************************
x <- "CRTH2"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD95---- OK
#***************************************************
x <- "CD95"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CCR7---- 8 feb, bruker fix grense på 2
#***************************************************
x <- "CCR7"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

#innfører fix grense på 2
split$lower_gates <- rep(2, length(split$lower_gates))

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg ICOS---- OK
#***************************************************
x <- "ICOS"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg NKG2A---- OK
#***************************************************
x <- "NKG2A"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD169---- OK?
#***************************************************
x <- "CD169"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD38----   8 feb må ha høyere minimum,  prøv 2.. SJEKK
#***************************************************
x <- "CD38"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 2)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD141---- OK   
#***************************************************
x <- "CD141"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg PD-1----OK
#***************************************************
x <- "PD-1"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}



data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))



split <- NA


#***************************************************
#pos/neg CD14---- OK  8 feb.. litt usikker men beholder minimum 1...
#***************************************************
x <- "CD14"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][2] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD11b---- her velger vi den lille toppen helt til h?yre ca 4, resten negativt, feb 8 prøv minim 3, sjekk... MÅ JOBBES MED!
#***************************************************
x <- "CD11b"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 3)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 4] <- mean(split$lower_gates[split$lower_gates < 4])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CCR6---- 8 feb ok
#***************************************************
x <- "CCR6"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 2.5] <- mean(split$lower_gates[split$lower_gates < 2.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CD 160 OK positiv=topp til h?yre
#***************************************************
x <- "CD160"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 3.5] <- mean(split$lower_gates[split$lower_gates < 3.5])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg HLADR---- OK positiv/negativ topp til h?yre
#***************************************************
x <- "HLADR"

params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 3] <- mean(split$lower_gates[split$lower_gates < 3])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
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
#pos/neg CXCR5 OK
#***************************************************
x <- "CXCR5"
params$desc[grep(x, params$desc)]#må sjekke at vi får riktig kanal.
kanal <- params$name[grep(x, params$desc)][1] 
data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = kanal) 



split <- find_gaussian_gates_second_top(data = data, channel = kanal, lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates > 2] <- mean(split$lower_gates[split$lower_gates < 2])

density_plots <- density_plot(data = data, channel = kanal, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)
#density_plots # to see the plots

tiff(fs::path(outFigPath, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()


data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(kanal, CD45))
kanal_max <- max(data[[1]][,kanal])
for(i in 1:n_files){
  kanal_max <- max(kanal_max, max(data[[i]][,kanal]))
}
signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD45, channel2 = kanal, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(0, kanal_max), title_size = 10)
plotTiff(signal = signal, filnavn = paste0("fig_", x, "_gating", ".tiff"))

for(j in filene){
  result[[j]][,x] <- data[[j]][, kanal] > split$lower_gates[j]
}




saveRDS(result, fs::path(outDataPath, "posNeg.rds"))
saveRDS(file_names, fs::path(outDataPath, "posNegFilnavn.rds"))

proc.time() - ptm
