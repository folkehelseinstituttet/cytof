resPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1")
d <- read.csv2( fs::path(resPath, "Prosent_Gating_Surrface_Panel1.csv"))
analyseSti <- "F:/Forskningsprosjekter/PDB 2794 - Immune responses aga_/Forskningsfiler/JOBO/CyTOF/Analyse i R OUS/"
meta <- readxl::read_xlsx(fs::path(analyseSti, "STATAfil OUS pasienter FINAL 270121.xlsx"))

#tar ut FHI005
d <- d[!grepl("FHI005", d$X),]
rownames(d) <- d$X
d <- d[,-1]

dInfo <- data.frame(name = rownames(d))
dInfo$name <- as.character(dInfo$name)

dInfo$Pat <- NA
dInfo$Tid <- NA
dInfo$Age <- NA
dInfo$status <- NA
for(i in 1:nrow(dInfo)){
  if(grepl("FHI", dInfo$name[i])){
    dInfo$Pat[i] <- paste("FHI", substr(strsplit(dInfo$name[i], "FHI")[[1]][2], 1,3))
    dInfo$Tid[i] <- paste0("T", substr(strsplit(dInfo$name[i], "T")[[1]][[2]], 1,1))
    dInfo$Age[i] <- substr(dInfo$name[i], 8,9)
    dInfo$status[i] <- substr(dInfo$name[i],1,1)
  }
}

dInfo$Pat[dInfo$Pat %in% "FHI 95_"] <- "FHI 095"
dInfo$Pat[dInfo$Pat %in% "FHI 81_"] <- "FHI 081"
dInfo$Age[dInfo$Age %in% "_6"] <- 62
dInfo$Age[dInfo$Age %in% "_2"] <- 22
dInfo$Age[dInfo$Age %in% "_1"] <- 102
dInfo$Age[dInfo$Age %in% "_8"] <- 82
dInfo$Age[dInfo$Age %in% "2_"] <- 82


dInfo$Pat[grepl("Ref1", dInfo$name)] <- "Ref1"
dInfo$Tid[grepl("Ref1", dInfo$name)] <- "Ref1"
dInfo$status[grepl("Ref1", dInfo$name)] <- "Ref1"

dInfo$Pat[grepl("Ko", dInfo$name)] <- paste0("Ko ", substr(dInfo$name[grepl("Ko", dInfo$name)], 12,14))
dInfo$Pat[dInfo$Pat == "Ko _55"] <- "Ko 550"
dInfo$Tid[grepl("Ko", dInfo$name)] <- "Kontroll"
dInfo$status[grepl("Ko", dInfo$name)] <- "Kontroll"


dInfo$Sex <- NA
dInfo$Sex[grepl("Fe", dInfo$name)] <- "Female"
dInfo$Sex[grepl("Ma", dInfo$name)] <- "Male"
dInfo$Sex[dInfo$name %in% "S_M_A_82_T2_FHI134_Panel1"] <- "Male"






#test between M and S, T1

dT1 <- d[dInfo$Tid == "T1",]
dInfoT1 <- dInfo[dInfo$Tid == "T1",]

resT1 <- as.data.frame(matrix(NA, ncol = ncol(d), nrow = 4))
rownames(resT1) <- c("median severe, all", "median moderat, all", "p all", "p adj. all")
colnames(resT1) <- colnames(dT1)


resT1[1,] <- apply(dT1[dInfoT1$status == "S",], 2, median)
resT1[2,] <- apply(dT1[dInfoT1$status == "M",], 2, median)
for(i in 1:ncol(resT1)){
  resT1[3,i] <- wilcox.test(dT1[dInfoT1$status == "S",i], dT1[dInfoT1$status == "M",i])$p.value
}
resT1[4,] <- p.adjust(resT1[3,])

tom <- rep(NA, ncol(resT1))

#  tid 1 og tid 2

Pat2 <- dInfo$Pat[duplicated(dInfo$Pat)]
Pat2 <- Pat2[!Pat2 %in% "Ref1"]

dBegge <- d[dInfo$Pat %in% Pat2, ]
dInfoBegge <- dInfo[dInfo$Pat %in% Pat2,]


resBegge <- as.data.frame(matrix(NA, ncol = ncol(d), nrow = 4))
rownames(resBegge) <- c("median T1, all", "median T2, all", "p all", "p adj. all")
colnames(resBegge) <- colnames(dBegge)


resBegge[1,] <- apply(dBegge[dInfoBegge$Tid == "T1",], 2, median)
resBegge[2,] <- apply(dBegge[dInfoBegge$Tid == "T2",], 2, median)
dBeggeT1 <- dBegge[dInfoBegge$Tid == "T1",]
dBeggeT2 <- dBegge[dInfoBegge$Tid == "T2",]
dInfoBeggeT1 <- dInfoBegge[dInfoBegge$Tid == "T1",]
dInfoBeggeT2 <- dInfoBegge[dInfoBegge$Tid == "T2",]
dInfoBeggeT1$Pat
dInfoBeggeT2$Pat  #sjekket at pasientrekkefølge er lik
for(i in 1:ncol(resBegge)){
  resBegge[3,i] <- wilcox.test(dBeggeT1[,i], dBeggeT2[,i], paired = T)$p.value
}
resBegge[4,] <- p.adjust(resBegge[3,])





