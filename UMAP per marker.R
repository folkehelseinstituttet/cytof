data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "clean data")
posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data")
fig_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "FigUMAP")
scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "clustering_functions.R"))


fcs_data_with_info <- read_data_from_folder(data_path)
fcs_data <- fcs_data_with_info$fcs_data
file_names <- fcs_data_with_info$file_names
rm(fcs_data_with_info)

params_fcs <- get_params_fcs_data(fcs_data[[1]])

posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))



kanaler <- params_fcs$desc[grepl("_",params_fcs$desc)]
kanaler <- kanaler[!grepl("_DNA", kanaler)]
kanaler <- kanaler[!grepl("_Cisp", kanaler)]
kanalnavn <- params_fcs$name[params_fcs$desc %in% kanaler]
print("antall markører med")
length(kanalnavn)
set.seed(100) # to ensure same plot every time
random_events_data <- random_events(number_of_events(fcs_data, file_names = file_names), n = 1000)

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data, 
                                                        kept_events = random_events_data, 
                                                        file_names = file_names, 
                                                        channels = kanalnavn,
                                                        archSin = T, 
                                                        cofactor = 5, 
                                                        scale = F)


d2 <- as.matrix(arcSindataMatrix[,1:(ncol(arcSindataMatrix) -2)])
posNegMat <- list_to_matrix_selected_events(data = posNeg, 
                                            kept_events = random_events_data, 
                                            file_names = posNegFilnavn, 
  #                                          channels = kanalnavn,
                                            archSin = F, cytofData = F
)


set.seed(100) # to ensure same plot every time
umapAshinfac5 <- umap::umap(d2)


par(mfrow=c(5,9))
for(i in 1:(ncol(posNegMat)-1)){
  tiff(fs::path(fig_path, paste0("UMAP_", colnames(posNegMat)[i], ".tiff")))
      print(plot(umapAshinfac5$layout, col = posNegMat[, i] + 1, main = colnames(posNegMat)[i]))
  dev.off()
}


### same for different datasett. 

data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Gating", "Gating fra R_FINAL")
posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_1", "Data")
fig_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_1", "FigUMAP")
scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "clustering_functions.R"))


fcs_data_with_info <- read_data_from_folder(data_path)
fcs_data <- fcs_data_with_info$fcs_data
file_names <- fcs_data_with_info$file_names
rm(fcs_data_with_info)

params_fcs <- get_params_fcs_data(fcs_data[[1]])

posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))



kanaler <- params_fcs$desc[grepl("_",params_fcs$desc)]
kanaler <- kanaler[!grepl("_DNA", kanaler)]
kanaler <- kanaler[!grepl("_Cisp", kanaler)]
kanalnavn <- params_fcs$name[params_fcs$desc %in% kanaler]
print("antall markører med")
length(kanalnavn)
set.seed(100) # to ensure same plot every time
random_events_data <- random_events(number_of_events(fcs_data, file_names = file_names), n = 1000)

arcSindataMatrix <- list_to_matrix_selected_events(data = fcs_data, 
                                                   kept_events = random_events_data, 
                                                   file_names = file_names, 
                                                   channels = kanalnavn,
                                                   archSin = T, 
                                                   cofactor = 5, 
                                                   scale = F)


d2 <- as.matrix(arcSindataMatrix[,1:(ncol(arcSindataMatrix) -2)])
posNegMat <- list_to_matrix_selected_events(data = posNeg, 
                                            kept_events = random_events_data, 
                                            file_names = posNegFilnavn, 
                                            #                                          channels = kanalnavn,
                                            archSin = F, cytofData = F
)


set.seed(100) # to ensure same plot every time
umapAshinfac5 <- umap::umap(d2)


par(mfrow=c(5,9))
for(i in 1:(ncol(posNegMat)-1)){
  tiff(fs::path(fig_path, paste0("UMAP_", colnames(posNegMat)[i], ".tiff")))
  print(plot(umapAshinfac5$layout, col = posNegMat[, i] + 1, main = colnames(posNegMat)[i]))
  dev.off()
}



