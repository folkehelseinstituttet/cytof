#posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1", "Data")
posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Data")

posNeg <- readRDS(fs::path(posNegPath, "posNeg.rds"))
posNegFileName <- readRDS(fs::path(posNegPath, "posNegFilnavn.rds"))

#outPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1")
outPath <-fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Stim")

scriptPath <- fs::path("H:", "git", "cytof")

source(fs::path(scriptPath, "analysis_functions.R"))

#surrface <- as.data.frame(readxl::read_excel(fs::path(outPath, "Gating stimulated panel i R.xlsx")))
surrface <- read.csv2(fs::path(outPath, "Gating stimulated panel i R.csv"))
surrface <- surrface[!is.na(surrface$ï..Population),]
surrface <- surrface[!surrface$ï..Population == "",]
if(any(duplicated(surrface))){
  surrface <- surrface[-duplicated(surrface),]
}

if(any(duplicated(surrface$Population))){
  print("Need unique population names")
  print(paste0("delete last entery of: ", surrface$Population[duplicated(surrface$Population)]))
  surrface <- surrface[!duplicated(surrface),]
}

surrface <- as.data.frame(surrface)
rownames(surrface) <- surrface$ï..Population
surrface <- surrface[,!colnames(surrface) == "Population"]

#
colnames(surrface)[!colnames(surrface) %in% colnames(posNeg[[1]])]
colnames(posNeg[[1]])[!colnames(posNeg[[1]]) %in% colnames(surrface)]

#change colnames so that they are equal between datasets
colnames(surrface)[colnames(surrface) == "TIM3"] <- "TIM-3"
colnames(surrface)[colnames(surrface) == "CD137/4IBB"] <- "CD137"
colnames(surrface)[colnames(surrface) == "IL12p70"] <- "IL-12p70"
colnames(surrface)[colnames(surrface) == "MIP1b"] <- "MIP-1b"
colnames(surrface)[colnames(surrface) == "BTLA"] <- "CD272"
colnames(surrface)[colnames(surrface) == "Granzyme.B"] <- "GranzymeB"
colnames(surrface)[colnames(surrface) == "LAG3"] <- "CD223"
colnames(surrface)[colnames(surrface) == "CD40L"] <- "CD154"
colnames(surrface)[colnames(surrface) == "IL.2"] <- "IL-2"
colnames(surrface)[colnames(surrface) == "PD.L1"] <- "PD-L1"
colnames(surrface)[colnames(surrface) == "CD137.4IBB"] <- "CD137"
colnames(surrface)[colnames(surrface) == "CTLA.4"] <- "CTLA-4"
colnames(surrface)[colnames(surrface) == "PD.L2"] <- "PD-L2"
colnames(surrface)[colnames(surrface) == "IL.17A"] <- "IL-17A"
colnames(surrface)[colnames(surrface) == "IL.1b"] <- "IL-1b"
colnames(surrface)[colnames(surrface) == "GM.CSF"] <- "GM-CSF"
colnames(surrface)[colnames(surrface) == "IL.6"] <- "IL-6"
colnames(surrface)[colnames(surrface) == "IL.10"] <- "IL-10"






extra_col <- colnames(surrface)[!colnames(surrface) %in% colnames(posNeg[[1]])]
if(length(extra_col) > 0){
  print(paste0("excel file includes other columns than dataset"))
  print(paste0("column: ", extra_col, ", are deleted"))
  surrface <- surrface[,!colnames(surrface) %in% extra_col]
}

if(any(rownames(surrface) == "possible values")){
  surrface <- surrface[-which(rownames(surrface) == "possible values"), ]
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
  if(length(either_or_i) ==1) {
    either_or_i <- NULL
  }
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

summary(result)

write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_Stim_panel2.csv"))

