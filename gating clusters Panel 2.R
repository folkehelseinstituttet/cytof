#posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1", "Data")
posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Data")

posNeg <- readRDS(fs::path(posNegPath, "posNeg.rds"))
posNegFileName <- readRDS(fs::path(posNegPath, "posNegFilnavn.rds"))

#outPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1")
outPath <-fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")

scriptPath <- fs::path("H:", "git", "cytof")

source(fs::path(scriptPath, "analysis_functions.R"))

#klustre <- as.data.frame(read.csv2(fs::path(outPath, "sjekk disse.csv")))
#klustre <- as.data.frame(read.csv2(fs::path(outPath, "Kopi av sjekk disse JB.csv")))
#klustre <- as.data.frame(read.csv2(fs::path(outPath, "sjekk disse CD8.csv")))
#klustre <- as.data.frame(read.csv2(fs::path(outPath, "sjekk disse CD8_1og2.csv")))
#klustre <- as.data.frame(read.csv2(fs::path(outPath, "sjekk disse fra Summary.csv")))
#klustre <- as.data.frame(read.csv2(fs::path(outPath, "gjennskapt.csv")))
klustre <- as.data.frame(read.csv2(fs::path(outPath, "gjennskapt kun less markers.csv")))
klustre <- klustre[-1,]
#klustre <- klustre[!is.na(klustre$Population),]
##klustre <- klustre[!klustre$CD45 == "0, 1",]
alleNA <- function(x){
  all(is.na(x))
}
alleTom <- function(x){
  all(x == "" | x == " ")
}

klustre <- klustre[!apply(klustre, 1, alleNA),]
klustre <- klustre[!apply(klustre, 1, alleTom),]
colnames(klustre)  
info <- klustre[,c(1, 2, 3, ncol(klustre)), ]  #nb sjekk
rownames(klustre) <- klustre$Poplulation
klustre <- klustre[, - c(1, 2,3, ncol(klustre))]

if(any(duplicated(klustre))){
  print(paste("duplisert kluster blir slettet:", as.character(info$Population[which(duplicated(klustre))])))
  info <- info[!duplicated(klustre), ]
  klustre <- klustre[!duplicated(klustre),]
}


#
colnames(klustre)[!colnames(klustre) %in% colnames(posNeg[[1]])]
colnames(posNeg[[1]])[!colnames(posNeg[[1]]) %in% colnames(klustre)]

#change colnames so that they are equal between datasets
colnames(klustre)[colnames(klustre) == "TIM3"] <- "TIM-3"
colnames(klustre)[colnames(klustre) == "CD137.4IBB"] <- "CD137"
colnames(klustre)[colnames(klustre) == "IL12p70"] <- "IL-12p70"
colnames(klustre)[colnames(klustre) == "MIP1b"] <- "MIP-1b"
colnames(klustre)[colnames(klustre) == "BTLA"] <- "CD272"
colnames(klustre)[colnames(klustre) == "Granzyme.B"] <- "GranzymeB"
colnames(klustre)[colnames(klustre) == "LAG3"] <- "CD223"
colnames(klustre)[colnames(klustre) == "IL.1b"] <- "IL-1b"
colnames(klustre)[colnames(klustre) == "IL.2"] <- "IL-2"
colnames(klustre)[colnames(klustre) == "PD.L1"] <- "PD-L1"
colnames(klustre)[colnames(klustre) == "IL.6"] <- "IL-6"
colnames(klustre)[colnames(klustre) == "GM.CSF"] <- "GM-CSF"
colnames(klustre)[colnames(klustre) == "IL.17A"] <- "IL-17A"
colnames(klustre)[colnames(klustre) == "IL.10"] <- "IL-10"
colnames(klustre)[colnames(klustre) == "CTLA.4"] <- "CTLA-4"
colnames(klustre)[colnames(klustre) == "PD.L2"] <- "PD-L2"
colnames(klustre)[colnames(klustre) == "CD40L"] <- "CD154"



extra_col <- colnames(klustre)[!colnames(klustre) %in% colnames(posNeg[[1]])]
if(length(extra_col) > 0){
  print(paste0("excel file includes other columns than dataset"))
  print(paste0("column: ", extra_col, ", are deleted"))
  klustre <- klustre[,!colnames(klustre) %in% extra_col]
}

if(any(rownames(klustre) == "possible values")){
  klustre <- klustre[-which(rownames(klustre) == "possible values"), ]
}

result <- matrix(NA, nrow = length(posNeg), ncol = nrow(klustre))
rownames(result) <- posNegFileName
colnames(result) <- paste(info$Population, info$navn)
for(i in 1:nrow(klustre)){
  pos_i <- colnames(klustre)[which(klustre[i,] == 1)]
  neg_i <- colnames(klustre)[which(klustre[i,] == 0)]
  pos12_i <- colnames(klustre)[which(klustre[i,] == 12)] 
  neg_01_i <- colnames(klustre)[which(klustre[i,] == 10)] 
  either_or_i <- colnames(klustre)[which(klustre[i,] == 99)]
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


#write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_.csv"))
#write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_Kopi av sjekk disse JB.csv"))
#write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_CD8.csv"))
#write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_CD8_1og2.csv"))
#write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_fra_Summary.csv"))
#write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_Gjenskapt.csv"))
write.csv2(as.data.frame(result), fs::path(outPath, "Prosent_Gating_klustre_Panel2_Gjenskapt_less_markers.csv"))
