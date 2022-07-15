
# this code is to be able to install packages from https://folkehelseinstituttet.github.io/drat/

#not sure if this is needed:
#usethis::edit_r_profile()

# write this into the file .Rprofile that comes up when runing the command "usethis::edit_r_profile()", 
# options(repos=structure(c(
#   "FHI=https://folkehelseinstituttet.github.io/drat/",
#   "CRAN=https://cran.rstudio.com"
# )))
# 
# save the .rprofile file, close R and reopen, then the next line should work. 

#install.packages("org")
#install.packages("flowCore")

library(ComplexHeatmap)
library(gridExtra)
library(ggplot2)
library(grid)
library(flowCore)
library(FlowSOM)
library(org)
library(data.table)
library(fs)
library(betareg)


#start with making a folder where you will analyse your project
path0 <- fs::path("C:", "Cytof data", "test")
#make a folder script in this folder and place all the files from GitHub in this folder

source(fs::path(path0, "script", "functions.R"))

#then to make the rest of the folderes within your project folder just run:
paths <- make_project_folders(path0)
rm(path0)
#now you have a directory system, place all your raw data in the folder "raw data"


#THIS HAS TO BE CHANGED ACCORDING TO NUMBER OF FILES IN YOUR DATASET.
#the plot_list have to include as many memberes as there are files in your dataset, also remember to update ncol and nrow.
plotSignal <- function(plot_list){
  g <- gridExtra::grid.arrange(plot_list[[1]], plot_list[[2]],
                               plot_list[[3]], plot_list[[4]], plot_list[[5]],
                               plot_list[[6]],  plot_list[[7]], plot_list[[8]],
                               plot_list[[9]], plot_list[[10]], plot_list[[11]],
                               plot_list[[12]], plot_list[[13]], plot_list[[14]],
                               plot_list[[15]], #plot_list[[16]], plot_list[[17]],
                               # plot_list[[18]], plot_list[[19]], plot_list[[20]],
                               # plot_list[[21]], plot_list[[22]], plot_list[[23]],
                               # plot_list[[24]], plot_list[[25]], plot_list[[26]],
                               # plot_list[[27]], plot_list[[28]], plot_list[[29]],
                               # plot_list[[30]], plot_list[[31]], plot_list[[32]],
                               # plot_list[[33]], plot_list[[34]], plot_list[[35]],
                               # plot_list[[36]], plot_list[[37]], plot_list[[38]],
                               # plot_list[[39]], plot_list[[40]], plot_list[[41]],
                               # plot_list[[42]], plot_list[[43]], plot_list[[44]],
                               # plot_list[[45]], plot_list[[46]], plot_list[[47]],
                               # plot_list[[48]], plot_list[[49]], plot_list[[50]],
                               # plot_list[[51]], plot_list[[52]], plot_list[[53]],
                               # plot_list[[54]], plot_list[[55]], plot_list[[56]],
                               # plot_list[[57]], plot_list[[58]], plot_list[[59]],
                               # plot_list[[60]], plot_list[[61]], plot_list[[62]],
                               # plot_list[[63]], plot_list[[64]], plot_list[[65]],
                               # plot_list[[66]],  ncol = 11, nrow = 6)
                               ncol = 5,nrow = 3)
  return(g)
}
###################
# CLEAN Up ----
###################
# The file "clean up gating all files.R" has to be updated according to existing data, 
# The function plotSignal make a plot per file, and have to be adjusted to include the correct number of files. 
# NB do not have  to be run each time. The if sentence ensure that the cleanup only are run when the folder is empty.
# if you still want to rerun skip the if sentence and only run source(..)
# Remember to look at the figures in paths$clean_data_fig_density_path and paths$clean_data_fig_signal_path

if(nrow(file.info(list.files(paths$clean_data_path))) == 0){
  source(fs::path(paths$script_path, "clean up gating all files.R"))
}
# when running line 114 you get many comments from R:
#Coordinate system already present. Adding new coordinate system, which will replace the existing one.
# Sabin: I have not understood how to not get this comments. Please feal free to think about it. 


# marker names are not "nice" in R, manually change to shorter names in ...\clean_up\clean_data_info\marker_names_included_manual_shortnames.csv. (save also with different name if you want to rerun cleanup)
# make sure that only markers are left, if not something has to be done in the end of the file "clean up gating all files.R"

###################
# Read clean data ----
###################

files_to_open <- "all" # this will include all files in the clean data folder, can be changed to a vector of filenames to use.
#rownames(file.info(list.files(paths$clean_data_path))) give you a list of all files in the clean data folder, only .fcs files are read

fcs_data_with_info <- read_data_from_folder(data_path = paths$clean_data_path, 
                                            files_to_open = files_to_open)
fcs_data <- fcs_data_with_info$fcs_data
file_names <- fcs_data_with_info$file_names
rm(fcs_data_with_info)

params_fcs <- get_params_fcs_data(fcs_data[[1]])


marker_info <- read.csv2(fs::path(paths$clean_data_info_path, "marker_names_included_manual_shortnames.csv"), stringsAsFactors = FALSE)

###################
# Marker gating ----
###################
# this section can be skipped if manual gates for the markers are not of interest. 
# but is needed for some markers if selected events should be used in further analysis. 
# this is not necessary to rerun. Hence the if sentence. If you still want to rerun just run the line within the if sentence. 
# this example only include gating for CD3, CD4, CD8 and CD45.
# all gates can be adjusted inside the file and manually displayed
if(nrow(file.info(list.files(paths$marker_gating_results_path)) ) == 0){
  source(fs::path(paths$script_path, "marker gating all files.R"))
}

# when running line 114 you get many comments from R:
#Coordinate system already present. Adding new coordinate system, which will replace the existing one.
# Sabin: I have not understood how to not get this comments. Please feal free to think about it. 



###################
# Run FlowSOM ----
###################

run_flowSOM(fcs_data = fcs_data, # clean data, not transformed
            file_names = file_names, # file_names in fcs_data in same order as the data
            included_files = file_names, # could be changed to only some files
            n_per_file = 15000, # number of events to include per file
            included_markers = marker_info$marker_name, # use all markers, can be manually changed with a list of markers to use in analysis.
            transformation = "arc_sinh", # by now only option. 
            scaling_flowSOM = TRUE, # scaling is recommanded, but will mean that all markers count as much regardless of size
            k_s = c(10, 15, 20), # number of clusters, k_s has to be a number or a vector of numbers.
            xdim = 10, # xdim * ydim gives the number of nodes that FlowSOM will combined to k_s clusters.
            ydim = 10, 
            resultpath = paths$clean_data_flowSOM_results_path,  # can be changed if new folder is wanted a new for different analysis, the folder have to exists before running function. recommand one folder per selected event analysis.
            seed = 1234, # the seed ensure that the same events are chosen, and hence the same result are uptained next time the exactly same analysis are done
            selectedEvents = "all", # if not "all" has to be the same name as a coloum-name in the matrices in posNeg
            posNeg = NULL, # has to have same structure as fcs_data, creat by....
            make_heatmap = TRUE, # need library ComplexHeatmap to make a heatmap. 
            heatmap_cluster_column = FALSE # if FALSE the column are not clusters and will appear in the same order.
            )
  
  
# this should be repeated for different seeds, and if wanted different selectedEvents (e.g. only CD3, Cd45, CD8 positive)


posNeg <- readRDS( fs::path(paths$marker_gating_results_path, "posNeg.rds"))

posNeg <- extra_column_posneg(posNeg = posNeg, column_name ="CD45CD3CD8", markers = c("CD45", "CD3", "CD8"), markers_values = c(1,1,12)) 

number_of_positive_events(posNeg = posNeg, marker = "CD45CD3CD8")

run_flowSOM(fcs_data = fcs_data, # clean data, not transformed
            file_names = file_names, # file_names in fcs_data in same order as the data
            included_files = "all", # could be changed to only some files
            n_per_file = 5000, # number of events to include per file
            included_markers = "all", # use all markers, can be manually changed with a list of markers to use in analysis.
            transformation = "arc_sinh", # by now only option. 
            scaling_flowSOM = TRUE, # scaling is recommanded, but will mean that all markers count as much regardless of size
            k_s = c(10, 15, 20), # number of clusters, k_s has to be a number or a vector of numbers.
            xdim = 10, # xdim * ydim gives the number of nodes that FlowSOM will combined to k_s clusters.
            ydim = 10, 
            resultpath = fs::path(paths$clean_data_flowSOM_results_path,  "CD45CD3CD8"),   
            seed = 2134, # the seed ensure that the same events are chosen, and hence the same result are uptained next time the exactly same analysis are done
            selectedEvents = "CD45CD3CD8", # if not "all" has to be the same name as a coloum-name in the matrices in posNeg
            posNeg = posNeg, # has to have same structure as fcs_data, creat by marker gating....
            make_heatmap = TRUE, # need library ComplexHeatmap to make a heatmap. 
            heatmap_cluster_column = FALSE # if FALSE the column are not clusters and will appear in the same order.
)

###################
# Markerplot ----
###################

resultpath <- paths$clean_data_flowSOM_results_path  # path to where the q5, q10, ... files from run_flowSOM is saved.
# if gates included in markerplots the csv file has to be made, colnames "X"	"low"	"high"
# X has to be equal column marker_short_name in ...clean_up\clean_data_info\marker_names_included_manual_shortnames.csv
# low is first gate, high second gate, both can be NA
# gates <- read.csv2(fs::path(paths$clean_data_posNeg_path, "mean_gates.csv"), stringsAsFactors = FALSE)

k <- 10
seed <- 1234
selectedEvents <- "all"
highlight_cluster <- 4
order_marker_shortname <- marker_info$marker_short_name
gates <- NULL

tiff(fs::path(resultpath, paste0("markerplot_k_", k,"_seed", seed, selectedEvents, ".tiff")), width = 1150, height = 900)
marker_plot(path = resultpath, 
            k = k, 
            seed = seed, 
            selectedEvents = selectedEvents, 
            highlight_cluster = highlight_cluster, 
            gates = gates, 
            order_marker_shortname = order_marker_shortname)
dev.off()




k <- 15
seed <- 1234
selectedEvents <- "all"
highlight_cluster <- NA
order_marker_shortname <- marker_info$marker_short_name
gates <- NULL
tiff(fs::path(resultpath, paste0("markerplot_k_", k,"_seed", seed, selectedEvents, ".tiff")), width = 1150, height = 900)
marker_plot(path = resultpath, 
            k = k, 
            seed = seed, 
            selectedEvents = selectedEvents, 
            highlight_cluster = highlight_cluster, 
            gates = gates, 
            order_marker_shortname = order_marker_shortname)
dev.off()



k <- 20
seed <- 1234
selectedEvents <- "all"
highlight_cluster <- 1:20
order_marker_shortname <- marker_info$marker_short_name
gates <- NULL
tiff(fs::path(resultpath, paste0("markerplot_k_", k,"_seed", seed, selectedEvents, ".tiff")), width = 1150, height = 900)
marker_plot(path = resultpath, 
            k = k, 
            seed = seed, 
            selectedEvents = selectedEvents, 
            highlight_cluster = highlight_cluster, 
            gates = gates, 
            order_marker_shortname = rker_info$marker_short_name)
dev.off()


#
k <- 10
seed <- 2134
selectedEvents <- "CD45CD3CD8"
highlight_cluster <- 4
order_marker_shortname <- marker_info$marker_short_name
gates <- read.csv2(fs::path(paths$marker_gating_results_path, "mean_gates.csv"))
resultpath <- fs::path(paths$clean_data_flowSOM_results_path,  "CD45CD3CD8")

tiff(fs::path(resultpath, paste0("markerplot_k_", k,"_seed", seed, selectedEvents, ".tiff")), width = 1150, height = 900)
marker_plot(path = resultpath, 
            k = k, 
            seed = seed, 
            selectedEvents = selectedEvents, 
            highlight_cluster = highlight_cluster, 
            gates = gates, 
            order_marker_shortname = order_marker_shortname)
dev.off()










# 
# 
# #calculate percentage based on cluster markers for all events. Can also be used to manually set lower and upper limites. lower or upper limites equal to NA is allowed.
# data <-  transform_selected_markers(fcs_data, markers = marker_name, method = "arc_sinh") 
# 
# 
# marker_info <- read.csv2(fs::path(paths$clean_data_info_path, "marker_names_included_manual_shortnames.csv"))
# k <- 15
# seed <- 1234
# cluster <- 8
# temp <- t(read.csv2(fs::path(resultpath, paste0("q5_per_cluster_k_", k, "_seed", seed, selectedEvents, ".csv")))[cluster, -1])
# lower <- data.frame(marker_name = rownames(temp), lower = temp)
# colnames(lower) <- c("marker_name", "lower")
# temp <- t(read.csv2(fs::path(resultpath, paste0("q95_per_cluster_k_", k, "_seed", seed, selectedEvents, ".csv")))[cluster, -1])
# upper <- data.frame(marker_name = rownames(temp), upper = temp)
# colnames(upper) <- c("marker_name", "upper")
# temp <- merge(marker_info, lower)
# lower_upper_limites <- merge(temp, upper)
# files_to_open <- "all"
# selected_markers = "all"
# 
# temp <- get_cluster_from_quantiles(lower_upper_limites = lower_upper_limites, data = data, selected_markers = selected_markers)
# temp
# 
# 
# temp <- get_cluster_from_quantiles(lower_upper_limites = lower_upper_limites, data = data, selected_markers = temp$high_cor_markers, files= files_to_open, correlation = 0.05)
# temp
# 
# #examples on analysing proportion data
# #if testing differences between two groups use wilcox.test (t-test assume normal distribution which proportion data is not)
# 
# #regression use beta regression.
# 
# #if metadata consist of 2 factors: x_Factor1 and x_Factor2 and one continues variable x_continiuos
# set.seed(100)
# y <- temp$res_prop
# x_continiuos <- runif(n = length(y), min = 15, max = 79)
# x_factor1 <- sample(c("F", "M"), length(y), replace = T)
# x_factor1 <- factor(x_factor1, levels = c("F", "M"))
# x_factor2 <- rep("low", length(y))
# x_factor2[y > mean(y)] <- "high"
# x_factor2 <- factor(x_factor2, levels = c("low", "high"))
# d <- data.frame(y = y, x_continiuos = x_continiuos, x_factor1 = x_factor1, x_factor2 = x_factor2)
# 
# fit0 <- betareg(y ~ 1, data = d)
# fit1 <- betareg(y ~ x_factor1, data = d)
# fit2 <- betareg(y ~ x_factor2, data = d)
# fit3 <- betareg(y ~ x_continiuos, data = d)
# fit4 <- betareg(y ~ x_factor1 + x_factor2, data = d)
# fit5 <- betareg(y ~ x_factor1 + x_continiuos, data = d)
# fit6 <- betareg(y ~ x_factor2 + x_continiuos, data = d)
# fit7 <- betareg(y ~ x_factor1 + x_factor2 + x_continiuos, data = d)
# AIC(fit0, fit1, fit2, fit3, fit4, fit5, fit6, fit7)
# 
# 
# fit <- fit2 # enkleste modell med AIC < minste(AIC) + 2
# res <- summary(fit)
# print(res)
# #predData <- expand.grid(x_continiuos= min(d$x_continiuos):max(d$x_continiuos), x_factor1 = levels(d$x_factor1), x_factor2 = levels(d$x_factor2))
# predData <- expand.grid(x_factor2 = levels(d$x_factor2))
# 
# # ONLY INCLUDE DATA WHERE YOU HAVE DATA...
# # for(n in levels(d$x_factor2)){
# #   min_n <- min(d$x_continiuos[d$x_factor2 == n])
# #   max_n <- max(d$x_continiuos[d$x_factor2 == n])
# #   predData$x_continiuos[predData$x_factor2 == n & predData$x_continiuos > max_n] <- NA
# #   predData$x_continiuos[predData$x_factor2 == n & predData$x_continiuos < min_n] <- NA
# # }
# pred_intervall <- 50
# nedre <- ((100 - pred_intervall)/2)/100
# ovre <- 1 - nedre
# 
# 
# predData$predict <- predict(fit, newdata = predData) #usikker på hvordan jeg får konfidensintervall. Vent til det skal brukes
# predData$navn <- predData$predict
# predData$pred_lower <- predict(fit, newdata = predData, type = "quantile", at = nedre)
# predData$pred_upper <- predict(fit, newdata = predData, type = "quantile", at = ovre)
# predData <- data.frame(predData)
# #i_x_continiuos <- which(colnames(d) == "x_continiuos")
# i_x_factor2 <- which(colnames(d) == "x_factor2")
# g <- ggplot(data = d, aes(x = x_continiuos, y = navn, col= x_factor2)) +
#   geom_point(size = 2) + 
#   geom_line(data = predData, aes(x=x_continiuos, y = navn, col = x_factor2), size = 2) +
#   geom_ribbon(data = predData, aes(ymin = pred_lower, ymax = pred_upper, fill = x_factor2, color = NULL), alpha = 0.1)
# }
# 
