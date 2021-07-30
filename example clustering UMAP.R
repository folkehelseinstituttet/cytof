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
params$desc[grepl("CD15", params$desc)]
VarNavn11 <- params$name[grepl("CD15", params$desc)][1]  #NB pass på å velge riktig....
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
  colnames(arcsinhexprData[[i]]) <- c("CD3", "CD19", "CD4", "CCR7", "CD45RA", "CD8", "CD25", "CD127", "CD56", "CD16", "CD15", "CD11c", "HLA-DR", "Time")
}



arcsinhexprDataMatrix <- list_to_matrix(data = arcsinhexprData, file_names = file_names)

# data here only include those column used for clustering (not time and dataset, could also exclude more)
data <- as.matrix(arcsinhexprDataMatrix[,1:(ncol(arcsinhexprDataMatrix) -2)])
#clustering using FlowSOM.
out <- FlowSOM::ReadInput(as.matrix(data), transform = F, scale = F)
out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(data)))
out <- FlowSOM::BuildMST(out)
cluster_FlowSOM_pre <- out$map$mapping[, 1]
#Modified k for each data set, try differnt k?
outk15 <-  FlowSOM::metaClustering_consensus(out$map$codes, k = 15, seed = 1) 
cluster_FlowSOM_15 <-outk15[cluster_FlowSOM_pre]

#Modified k for each data set,
outk28 <-  FlowSOM::metaClustering_consensus(out$map$codes, k = 28, seed = 1) 
cluster_FlowSOM_28 <-outk28[cluster_FlowSOM_pre]
# FlowSOM do one clustering of 100 clusters, then it is grouped down to wished number of clusters k.
# by tabulate two of these we see which clusters that are grouped further to get fewer clusters
table(cluster_FlowSOM_15, cluster_FlowSOM_28)




#clustering with k-means
k_mean_clustering_different_ks(data, k =c(1,3,5,10,15,20))
# choose the number of clusters where the prop_within fail to decrease further. 
Cluster_kmean15 <- kmeans(data, centers = 15)$cluster

#compare how two clustering methods differs
table(Cluster_kmean15, cluster_FlowSOM_15)
table(Cluster_kmean15, cluster_FlowSOM_28)


#umapAshinfac5 <- umap::umap(data) # take to much time if all data are plottet.
#Choose to plot only N observations per data set, and repeat 3 times with different seeding to see differences.

set.seed(100) # to ensure same plot every time
random_cells_for_plotting <- random_cells_vector(arcsinhexprDataMatrix$dataset)
umapAshinfac5 <- umap::umap(data[random_cells_for_plotting,])

plot(umapAshinfac5$layout, col = col25[cluster_FlowSOM_15[random_cells_for_plotting]])
plot(umapAshinfac5$layout, col = col25[cluster_FlowSOM_15[random_cells_for_plotting]], xlim = c(-10,10), ylim = c(-14,10))
plot(umapAshinfac5$layout, col = col25[cluster_FlowSOM_28[random_cells_for_plotting]], xlim = c(-10,10), ylim = c(-14,10))



density_plot_per_cluster(data = data, cluster_per_cell = cluster_FlowSOM_15, rand_cells = random_cells_for_plotting)
#density_plot_per_cluster(data = data, cluster_per_cell = cluster_FlowSOM_15, rand_cells = random_cells_for_plotting, plot_cluster = c(1,2))

barplot_per_sample(file_names_per_cell = arcsinhexprDataMatrix$dataset,  cluster_per_cell = cluster_FlowSOM_15, rand_cells = random_cells_for_plotting)
