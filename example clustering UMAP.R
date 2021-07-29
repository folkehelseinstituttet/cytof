data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2770 - EXIMIOUS - Mapping E_", "CyTOF in vitro WP4", "GatedData")
scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "read_data.R"))
source(fs::path(scriptPath, "transformation.R"))
source(fs::path(scriptPath, "ploting.R"))

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


arcsinhexprDataMatrix <- rbind(arcsinhexprData[[1]], arcsinhexprData[[2]],arcsinhexprData[[3]],
                               arcsinhexprData[[4]], arcsinhexprData[[5]],arcsinhexprData[[6]],
                               arcsinhexprData[[7]], arcsinhexprData[[8]])

arcsinhexprDataMatrix$dataset <- c(rep(file_names[1], nrow(arcsinhexprData[[1]])), 
                                   rep(file_names[2], nrow(arcsinhexprData[[2]])),
                                   rep(file_names[3], nrow(arcsinhexprData[[3]])),
                                   rep(file_names[4], nrow(arcsinhexprData[[4]])),
                                   rep(file_names[5], nrow(arcsinhexprData[[5]])),
                                   rep(file_names[6], nrow(arcsinhexprData[[6]])),
                                   rep(file_names[7], nrow(arcsinhexprData[[7]])),
                                   rep(file_names[8], nrow(arcsinhexprData[[8]])))


data <- as.matrix(arcsinhexprDataMatrix[,1:(ncol(arcsinhexprDataMatrix) -2)])
out <- FlowSOM::ReadInput(as.matrix(data), transform = F, scale = F)
out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(data)))
out <- FlowSOM::BuildMST(out)
labels_pre <- out$map$mapping[, 1]
#Modified k for each data set,
outk14 <-  FlowSOM::metaClustering_consensus(out$map$codes, k = 14, seed = 1) 
labels_14 <-outk14[labels_pre]

#Modified k for each data set,
outk28 <-  FlowSOM::metaClustering_consensus(out$map$codes, k = 28, seed = 1) 
labels_28 <-outk28[labels_pre]



#umapAshinfac5 <- umap::umap(data) # take to much time if all data are plottet.
#Choose to plot only N observations per data set, and repeat 3 times to see differences.

set.seed(100) # to ensure same plot every time
random_cells_for_plotting <- random_cells_vector(arcsinhexprDataMatrix$dataset)
umapAshinfac5 <- umap::umap(data[random_cells_for_plotting,])

plot(umapAshinfac5$layout, col = col25[labels_14[random_cells_for_plotting]])
plot(umapAshinfac5$layout, col = col25[labels_14[random_cells_for_plotting]], xlim = c(-10,10), ylim = c(-15,10))
plot(umapAshinfac5$layout, col = col25[labels_28[random_cells_for_plotting]], xlim = c(-10,10), ylim = c(-15,10))


