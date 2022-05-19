
if(selectedEvents == T){
  set.seed(seed)
  random_events_data <- random_events_from_selected_events(posNeg = posNeg, marker = channel, n = n_per_file)
  ext_name <- channel
} else {
  set.seed(seed)
  number_of_events_data <-  number_of_events(data = fcs_data, file_names = file_names)
  random_events_data <- random_events(number_of_events_data, n = n_per_file)
  ext_name <- "all"
}


arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data, 
                                                          kept_events = random_events_data, 
                                                          file_names = file_names, 
                                                          channels = kanalnavn,
                                                          archSin = T, 
                                                          cofactor = 5, 
                                                          scale = F
)
colnames(arcSindataMatrix)[1:length(kanaler)] <- kanaler
print("arcsin")

if(column_cluster == TRUE){
  o <- kanaler
}

params <- list(seed = seed,
               data = arcSindataMatrix[arcSindataMatrix$dataset %in% tamed, ],
               n_random_for_plotting_per_fil = 5000,
               kanaler = kanaler,
               scaling = TRUE,
               column_cluster = column_cluster, 
               o <- o,
               ydim = xdim,  
               xdim = ydim
)

set.seed(params$seed) # to ensure same plot every time
out <- FlowSOM::ReadInput(as.matrix(params$data[,params$kanaler]), transform = F, scale = params$scaling)
print(1)
out <- FlowSOM::BuildSOM(out, colsToUse = 1:(ncol(params$data[,params$kanaler])), xdim = params$xdim, ydim = params$ydim)
print(2)
out <- FlowSOM::BuildMST(out)
print(3)
cluster_FlowSOM_pre <- out$map$mapping[, 1]
set.seed(params$seed)

for(k in ks){
  print(k)
  out_k <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = params$seed)
  cluster_FlowSOM_k <- out_k[cluster_FlowSOM_pre]
  cluster_FlowSOM_k_factor <- factor(cluster_FlowSOM_k, levels = 1:k)
print("cluster")
# ved å inkludere disse ble ikke datasettet likt. 
  q5_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.05)
  write.csv2(q5_k, fs::path(utSti, paste0("q5_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
  q10_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.1)
  write.csv2(q10_k, fs::path(utSti, paste0("q10_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
  q25_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.25)
  write.csv2(q25_k, fs::path(utSti, paste0("q25_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
  q75_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.75)
  write.csv2(q75_k, fs::path(utSti, paste0("q75_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
  q90_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.9)
  write.csv2(q90_k, fs::path(utSti, paste0("q90_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
  q95_k <- q_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k, probs = 0.95)
  write.csv2(q95_k, fs::path(utSti, paste0("q95_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
  
print("q")

  medians_k <- median_per_cluster_marker(data = params$data[,params$kanaler], kluster = cluster_FlowSOM_k)
  write.csv2(medians_k, fs::path(utSti, paste0("medians_per_kluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
print("median")



  tiff(fs::path(utSti, paste0("heatmap_median_k_", k, "_cluster_seed", seed, ext_name, ".tiff")), width = 1000, height = 800)
    print(Heatmap(as.matrix(medians_k[,o], cluster_columns = params$column_cluster)))
  dev.off()
print("heatmap")


  totalt_antall_i_kluster <- table(cluster_FlowSOM_k_factor)
  

  antallPerKluster <- cbind(totalt_antall_i_kluster, totalt_antall_i_kluster/nrow(params$data))
  colnames(antallPerKluster) <- c("antall", "prosent")
  write.csv2(antallPerKluster, fs::path(utSti, paste0("antallPerCluster_k_", k, "_seed", params$seed, ext_name, ".csv")))
print("totalt antall")

  
  perSample <- table(params$data$dataset, cluster_FlowSOM_k_factor)
  colnames(perSample) <- paste0("kluster_", 1:k)

  write.csv2(perSample, fs::path(utSti, paste0("antallPerClusterOgPat_k_", k, "_seed", params$seed, ext_name, ".csv")))
print("antall per fil")



#  randomEvents <- sample(nrow(params$data), 10000)  #might have to rethink to ensure given cluster is represented among random cells. 
#   pdf(fs::path(utSti, paste0("density_plot_per_cluster_k_", k, "_cluster_seed", params$seed, ext_name, ".pdf")), width = 10, height = 7)
#   density_plot_per_cluster(data = params$data[randomEvents,params$kanaler], cluster_per_cell = cluster_FlowSOM_k_factor[randomEvents],  nrow_plot = 5, strip_text_size = 5)
#   dev.off()
# print("density")
}

  