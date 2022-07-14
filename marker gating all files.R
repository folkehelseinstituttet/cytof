

# Example
#  
# 
ptm <- proc.time()



result <- list()

for(j in 1:length(fcs_data)){
  mat <-  as.data.frame(matrix(NA, ncol = nrow(marker_info), nrow =  nrow(fcs_data[[j]])))
  colnames(mat) <- marker_info$markers_short_name
  result[[j]] <- mat
}

mean_gates <- as.data.frame(matrix(NA, nrow = nrow(marker_info), ncol = 2))
rownames(mean_gates) <- marker_info$markers_short_name
colnames(mean_gates) <- c("low", "high")


#***************************************************
#pos/neg CD3 ---- 
#***************************************************
CD45 <- marker_info$markers_name[marker_info$markers_short_name == "CD45"]
x <- "CD3"

marker_i <-  marker_info$markers_name[marker_info$markers_short_name == x]
data <-  transform_selected_markers(fcs_data = fcs_data, markers = c(marker_i, CD45), method = "arc_sinh") 

split <- find_gaussian_gates(data = data, marker = marker_i, method = "second top", lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, marker = marker_i, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)

tiff(fs::path(paths$marker_gating_fig_density_path, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

ymax <- max_marker(data = data, marker = marker_i)
ymin <- min_marker(data = data, marker = marker_i)

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), marker1 = CD45, marker2 = marker_i, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(ymin, ymax), title_size = 10)


tiff(fs::path(paths$marker_gating_fig_signal_path, paste0("fig_", x, "_gating", ".tiff")), width = 1800, height = 1200)
plotSignal(plot_list = signal$plotList)
dev.off()


for(j in 1:length(data)){
  result[[j]][,x] <- data[[j]][, marker_i] > split$lower_gates[j]
}


mean_gates[x,1] <- mean(split$lower_gates)

split <- NA


#***************************************************
#pos/neg CD4 ---- 
#***************************************************
x <- "CD4"
marker_i <-  marker_info$markers_name[marker_info$markers_short_name == x]
data <-  transform_selected_markers(fcs_data = fcs_data, markers = c(marker_i, CD45), method = "arc_sinh") 

split <- find_gaussian_gates(data = data, marker = marker_i, method = "second top", lower_gate_percent = 15, upper_gate_percent = 0.001, minimum = 1.5)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, marker = marker_i, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)

tiff(fs::path(paths$marker_gating_fig_density_path, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

ymax <- max_marker(data = data, marker = marker_i)
ymin <- min_marker(data = data, marker = marker_i)

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), marker1 = CD45, marker2 = marker_i, ylow = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(ymin, ymax), title_size = 10)


tiff(fs::path(paths$marker_gating_fig_signal_path, paste0("fig_", x, "_gating", ".tiff")), width = 1800, height = 1200)
plotSignal(plot_list = signal$plotList)
dev.off()


for(j in 1:length(data)){
  result[[j]][,x] <- data[[j]][, marker_i] > split$lower_gates[j]
}


mean_gates[x,1] <- mean(split$lower_gates)

split <- NA

#***************************************************
#pos/neg CD8 ----  example with two split one for noise...
#***************************************************
x <- "CD8"
marker_i <-  marker_info$markers_name[marker_info$markers_short_name == x]
data <-  transform_selected_markers(fcs_data = fcs_data, markers = c(marker_i, CD45), method = "arc_sinh") 

split <- find_gaussian_gates(data = data, marker = marker_i, method = "second top", lower_gate_percent = 15, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])
split$lower_gates[split$lower_gates < 2.5] <- mean(split$lower_gates[split$lower_gates > 2.5])

splitLow <- find_gaussian_gates(data = data, marker = marker_i, method = "lower noise")

density_plots <- density_plot(data = data, marker = marker_i, plot_title = file_names, lower_gate = splitLow, upper_gate =  split$lower_gates, main_title = x)

tiff(fs::path(paths$marker_gating_fig_density_path, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

ymax <- max_marker(data = data, marker = marker_i)
ymin <- min_marker(data = data, marker = marker_i)

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), marker1 = CD45, marker2 = marker_i, ylow = splitLow,  yhigh = split$lower_gates, xname = "CD45", yname = x, plot_title = file_names, ylim = c(ymin, ymax), title_size = 10)


tiff(fs::path(paths$marker_gating_fig_signal_path, paste0("fig_", x, "_gating", ".tiff")), width = 1800, height = 1200)
plotSignal(plot_list = signal$plotList)
dev.off()


for(j in 1:length(data)){
  result[[j]][,x] <- data[[j]][, marker_i] > splitLow[j]
  result[[j]][split$lower_gates[j],x] <- 2
}


mean_gates[x,1] <- mean(splitLow)
mean_gates[x,2] <- mean(split$lower_gates)

splitLow <- NA
split <- NA


split <- NA
#***************************************************
#pos/neg CD45 ---- 
#***************************************************
CD3 <- marker_info$markers_name[marker_info$markers_short_name == "CD3"]
x <- "CD45"
marker_i <-  marker_info$markers_name[marker_info$markers_short_name == x]
data <-  transform_selected_markers(fcs_data = fcs_data, markers = c(marker_i, CD3), method = "arc_sinh") 

split <- find_gaussian_gates(data = data, marker = marker_i, method = "second top", lower_gate_percent = 2, upper_gate_percent = 0.001)
split$lower_gates[is.na(split$lower_gates)] <- mean(split$lower_gates[!is.na(split$lower_gates)])

density_plots <- density_plot(data = data, marker = marker_i, plot_title = file_names, lower_gate = split$lower_gates, main_title = x)

tiff(fs::path(paths$marker_gating_fig_density_path, paste0("fig_", x, "_gating", ".tiff")), height = 1800, width = 600)
print(density_plots)
dev.off()

ymax <- max_marker(data = data, marker = marker_i)
ymin <- min_marker(data = data, marker = marker_i)

signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), marker1 = CD3, marker2 = marker_i, ylow = split$lower_gates, xname = "CD3", yname = x, plot_title = file_names, ylim = c(ymin, ymax), title_size = 10)


tiff(fs::path(paths$marker_gating_fig_signal_path, paste0("fig_", x, "_gating", ".tiff")), width = 1800, height = 1200)
plotSignal(plot_list = signal$plotList)
dev.off()


for(j in 1:length(data)){
  result[[j]][,x] <- data[[j]][, marker_i] > split$lower_gates[j]
}


mean_gates[x,1] <- mean(split$lower_gates)

split <- NA



saveRDS(result, fs::path(paths$marker_gating_results_path, "posNeg.rds"))
saveRDS(file_names, fs::path(paths$marker_gating_results_path, "posNegFilnavn.rds"))
write.csv2(mean_gates, fs::path(paths$marker_gating_results_path, "mean_gates.csv"))

proc.time() - ptm
