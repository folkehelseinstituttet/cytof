#posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1", "Data")
posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data")

resultposNeg <- readRDS(fs::path(posNegPath, "posNeg.rds"))
posNegFileName <- readRDS(fs::path(posNegPath, "posNegFilnavn.rds"))

#outPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1")
outPath <-fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "surrface")

scriptPath <- fs::path("H:", "git", "cytof")

source(fs::path(scriptPath, "analysis_functions.R"))

surrface <- as.data.frame(readxl::read_excel(fs::path(outPath, "Gating surface i R.xlsx")))
surrface <- surrface[!is.na(surrface$Population),]
if(any(duplicated(surrface))){
  surrface <- surrface[-duplicated(surrface),]
}

if(any(duplicated(surrface$Population))){
  print("Need unique population names")
  print(paste0("delete last entery of: ", surrface$Population[duplicated(surrface$Population)]))
  surrface <- surrface[!duplicated(surrface),]
}

surrface <- as.data.frame(surrface)
rownames(surrface) <- surrface$Population
surrface <- surrface[,!colnames(surrface) == "Population"]

#
colnames(surrface)[!colnames(surrface) %in% colnames(posNeg[[1]])]
colnames(posNeg[[1]])[!colnames(posNeg[[1]]) %in% colnames(surrface)]

#change colnames so that they are equal between datasets
colnames(surrface)[colnames(surrface) == "PD1"] <- "PD-1"
colnames(surrface)[colnames(surrface) == "OX40"] <- "CD134_OX40"


extra_col <- colnames(surrface)[!colnames(surrface) %in% colnames(posNeg[[1]])]
if(length(extra_col) > 0){
  print(paste0("excel file includes other columns than dataset"))
  print(paste0("column: ", extra_col, ", are deleted"))
  surrface <- surrface[,!colnames(surrface) %in% extra_col]
}

if(any(rownames(surrface) == "possible combinations")){
  surrface <- surrface[-which(rownames(surrface) == "possible combinations"), ]
}

result <- matrix(NA, nrow = length(posNeg), ncol = nrow(surrface))
rownames(result) <- posNegFileName
colnames(result) <- rownames(surrface)
for(i in 1:nrow(surrface)){
  pos_i <- colnames(surrface)[which(surrface[i,] == 1)]
  neg_i <- colnames(surrface)[which(surrface[i,] == 0)]
  pos12_i <- colnames(surrface)[which(surrface[i,] == 12)] 
  neg_01_i <- colnames(surrface)[which(surrface[i,] == 10)] 
  either_or_i <- colnames(surrface)[which(surrface[i,] == 99)]
  markers_i <- c(pos_i, neg_i, pos12_i, neg_01_i)
  markersPosNeg_i <- c(rep(1, length(pos_i)), rep(0, length(neg_i)), rep(12, length(pos12_i)), rep(10, length(neg_01_i)))
  if(length(markers_i) > 0){
    if(length(either_or_i) == 0){
      result[,i] <- prosent_senario(posneg = posNeg, channels = markers_i, values = markersPosNeg_i)
    } else {
      result[,i] <- prosent_senario_with_atleat_one_of_some_channels(posneg = posNeg, channels = markers_i, 
                                                                     values = markersPosNeg_i, atleast_one_of = either_or_i, value_atleast_one_of = 12)		
    }
  }
  
}


write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_Surrface_Panel1.csv"))

