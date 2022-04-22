posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Data")
fig_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022")

scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "clustering_functions.R"))
source(fs::path(scriptPath, "analysis_functions.R"))

posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))


resultTighter <- matrix(NA, nrow = length(posNeg), ncol = ncol(posNeg[[1]]))
rownames(resultTighter) <- posNegFilnavn
colnames(resultTighter) <- colnames(posNeg[[1]])


for(i in colnames(resultTighter)){
  maks <- max(posNeg[[1]][,i])
  resultTighter[, i] <- prosent_senario(posneg = posNeg, channels = i, values = maks)
}

rownames(resultTighter)

resultTighter <- as.data.frame(resultTighter)

resultTighter$filnavn <- rownames(resultTighter)

resultTighterLong <- data.table::melt(data.table::setDT(resultTighter),id.vars = "filnavn", variable.names = "marker")
colnames(resultTighterLong)[3] <- "Tighter"

resultTighterLong$filnavn_marker <- paste0(resultTighterLong$filnavn, resultTighterLong$variable)


posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_1", "Data")

scriptPath <- fs::path("H:", "git", "cytof")


posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))


resultWider <- matrix(NA, nrow = length(posNeg), ncol = ncol(posNeg[[1]]))
rownames(resultWider) <- posNegFilnavn
colnames(resultWider) <- colnames(posNeg[[1]])


for(i in colnames(resultWider)){
  maks <- max(posNeg[[1]][,i])
  resultWider[, i] <- prosent_senario(posneg = posNeg, channels = i, values = maks)
}

resultWider <- as.data.frame(resultWider)

resultWider$filnavn <- rownames(resultWider)

resultWiderLong <- data.table::melt(data.table::setDT(resultWider),id.vars = "filnavn", variable.names = "marker")
colnames(resultWiderLong)[3] <- "Wider"

resultWiderLong$filnavn_marker <- paste0(resultWiderLong$filnavn, resultWiderLong$variable)


result <- merge(resultTighterLong, resultWiderLong, by = "filnavn_marker")
result <- result[,-c(1, 5,6)]
colnames(result)[1] <- "filename"
colnames(result)[2] <- "marker"



library(ggplot2)
g <- ggplot(result, aes(x = Wider, y = Tighter)) +
  geom_point() +
  facet_wrap(vars(marker), ncol = 11, scale = "free") + 
  geom_abline(slope = 1)
  

tiff(fs::path(fig_path, "thighterWider.tiff"), width = 1150, height = 900)
g
dev.off()



## panel 2
posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Data")
fig_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022")

scriptPath <- fs::path("H:", "git", "cytof")
source(fs::path(scriptPath, "read_data_functions.R"))
source(fs::path(scriptPath, "transformation_functions.R"))
source(fs::path(scriptPath, "ploting_functions.R"))
source(fs::path(scriptPath, "clustering_functions.R"))
source(fs::path(scriptPath, "analysis_functions.R"))

posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))


resultTighter <- matrix(NA, nrow = length(posNeg), ncol = ncol(posNeg[[1]]))
rownames(resultTighter) <- posNegFilnavn
colnames(resultTighter) <- colnames(posNeg[[1]])


for(i in colnames(resultTighter)){
  maks <- max(posNeg[[1]][,i])
  resultTighter[, i] <- prosent_senario(posneg = posNeg, channels = i, values = maks)
}

rownames(resultTighter)

resultTighter <- as.data.frame(resultTighter)

resultTighter$filnavn <- rownames(resultTighter)

resultTighterLong <- data.table::melt(data.table::setDT(resultTighter),id.vars = "filnavn", variable.names = "marker")
colnames(resultTighterLong)[3] <- "Tighter"

resultTighterLong$filnavn_marker <- paste0(resultTighterLong$filnavn, resultTighterLong$variable)


posNeg_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2", "Data")

scriptPath <- fs::path("H:", "git", "cytof")


posNegFilnavn <- as.character(readRDS(fs::path(posNeg_path,  "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(posNeg_path,  "posNeg.rds"))


resultWider <- matrix(NA, nrow = length(posNeg), ncol = ncol(posNeg[[1]]))
rownames(resultWider) <- posNegFilnavn
colnames(resultWider) <- colnames(posNeg[[1]])


for(i in colnames(resultWider)){
  maks <- max(posNeg[[1]][,i])
  if(!is.na(maks))
    resultWider[, i] <- prosent_senario(posneg = posNeg, channels = i, values = maks)
}

resultWider <- as.data.frame(resultWider)

resultWider$filnavn <- rownames(resultWider)

resultWiderLong <- data.table::melt(data.table::setDT(resultWider),id.vars = "filnavn", variable.names = "marker")
colnames(resultWiderLong)[3] <- "Wider"

resultWiderLong$filnavn_marker <- paste0(resultWiderLong$filnavn, resultWiderLong$variable)


result <- merge(resultTighterLong, resultWiderLong, by = "filnavn_marker")
result <- result[,-c(1, 5,6)]
colnames(result)[1] <- "filename"
colnames(result)[2] <- "marker"



library(ggplot2)
g <- ggplot(result, aes(x = Wider, y = Tighter)) +
  geom_point() +
  facet_wrap(vars(marker), ncol = 11, scale = "free") + 
  geom_abline(slope = 1)


tiff(fs::path(fig_path, "thighterWider.tiff"), width = 1150, height = 900)
g
dev.off()

