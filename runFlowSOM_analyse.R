rm(list = ls())

library(ComplexHeatmap)
library(gridExtra)
library(ggplot2)
library(grid)


scriptPath <- fs::path("H:", "git", "cytof")
#scriptPath <- fs::path("C:", "CyToF data", "fra github")

source(fs::path(scriptPath, "readDataToAnalysePanel2.R"))

tamed <- as.character(dInfo$filenames)
tamed <- tamed[!grepl("FHI005", tamed)] #denne oppfører seg helt rart og tas ut av analysen.
tamed <- tamed[!grepl("Ref1", tamed)] #denne oppfører seg helt rart og tas ut av analysen.
tamed

ks <- c(40, 50, 60)
n_per_file <- 25000
xdim <- 14
ydim <- 14


utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2_ALLE")
seed <- 12345 #nb må endres vil man vil gjøre et annet uttrekk
selectedEvents <- FALSE
channel <- NULL
source(fs::path(scriptPath, "FlowSOM_analyse.R"))



for(j in 1:length(posNeg)){
  posNeg[[j]]$CD3CD45CD4 <- posNeg[[j]]$CD3 & posNeg[[j]]$CD45 & posNeg[[j]]$CD4
}

utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD4_test")
seed <- 45903#nb må endres vil man vil gjøre et annet uttrekk
selectedEvents <- TRUE
channel <- "CD3CD45CD4"
n_per_file <- 4250
ks <- c(15, 20, 25, 30, 35, 40)
xdim <- 10
ydim <- 10

source(fs::path(scriptPath, "FlowSOM_analyse.R"))





for(j in 1:length(posNeg)){
  posNeg[[j]]$CD3CD45CD8 <- posNeg[[j]]$CD3 & posNeg[[j]]$CD45 & posNeg[[j]]$CD8
}

utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD8")
seed <- 45902#nb må endres vil man vil gjøre et annet uttrekk
selectedEvents <- TRUE
channel <- "CD3CD45CD8"
n_per_file <- 4000
ks <- c(15, 20, 25, 30, 35, 40)
xdim <- 10
ydim <- 10


source(fs::path(scriptPath, "FlowSOM_analyse.R"))

#panel 1

source(fs::path(scriptPath, "readDataToAnalyse.R"))

tamed <- as.character(dInfo$filenames)
tamed <- tamed[!grepl("FHI005", tamed)] #denne oppfører seg helt rart og tas ut av analysen.
tamed <- tamed[!grepl("Ref1", tamed)] #denne oppfører seg helt rart og tas ut av analysen.
tamed

ks <- c(40, 50, 60)
n_per_file <- 25000
xdim <- 14
ydim <- 14


utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE")
seed <- 1245 #nb må endres vil man vil gjøre et annet uttrekk
selectedEvents <- FALSE
channel <- NULL
source(fs::path(scriptPath, "FlowSOM_analyse.R"))



for(j in 1:length(posNeg)){
  posNeg[[j]]$CD3CD45CD4 <- posNeg[[j]]$CD3 & posNeg[[j]]$CD45 & posNeg[[j]]$CD4
}

utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD4")
seed <- 4503#nb må endres vil man vil gjøre et annet uttrekk
selectedEvents <- TRUE
channel <- "CD3CD45CD4"
n_per_file <- 4250
ks <- c(15, 20, 25, 30, 35, 40)
xdim <- 10
ydim <- 10

source(fs::path(scriptPath, "FlowSOM_analyse.R"))





for(j in 1:length(posNeg)){
  posNeg[[j]]$CD3CD45CD8 <- posNeg[[j]]$CD3 & posNeg[[j]]$CD45 & posNeg[[j]]$CD8
}

utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD8")
seed <- 4502#nb må endres vil man vil gjøre et annet uttrekk
selectedEvents <- TRUE
channel <- "CD3CD45CD8"
n_per_file <- 4000
ks <- c(15, 20, 25, 30, 35, 40)
xdim <- 10
ydim <- 10


source(fs::path(scriptPath, "FlowSOM_analyse.R"))

