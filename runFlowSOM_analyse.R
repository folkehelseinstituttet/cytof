rm(list = ls())

library(ComplexHeatmap)
library(gridExtra)
library(ggplot2)
library(grid)


orderPanel1 <- c("X89Y_CD45", "X116Cd_CD3", "X145Nd_CD4", "X113Cd_CD8", "X111Cd_CD19", "X152Sm_TCRgd", "X166Er_TCRVa7.2", "X149Sm_CD25", 
       "X158Gd_CD27", "X160Gd_CD28", "X143Nd_CD127", "X172Yb_CD38", "X167Er_CCR7", "X155Gd_CD45RA", "X114Cd_HLADR", "X146Nd_IgD", "X159Tb_IgG",
       "X156Gd_CXCR3", "X153Eu_CCR4", "X141Pr_CCR6", "X171Yb_CXCR5", "X168Er_ICOS", "X142Nd_KLRG1", "X150Nd_CD134_OX40", "X154Sm_TIGIT",
       "X161Dy_CD160", "X164Dy_CD161", "X162Dy_CD95", "X163Dy_CRTH2", "X165Ho_CD85j", "X169Tm_NKG2A", "X173Yb_CD141", "X174Yb_CD279_PD.1",
       "X175Lu_CD14", "X148Nd_CD16", "X176Yb_CD56", "X106Cd_CD57", "X209Bi_CD11b", "X147Sm_CD11c", "X112Cd_CD5", "X144Nd_CD15", "X170Er_CD169", "X151Eu_CD123")

orderPanel2 <- c("X89Y_CD45", "X116Cd_CD3", "X113Cd_CD8", "X145Nd_CD4",  "X111Cd_CD19",  "X153Eu_CXCR5.CD185", 
       "X114Cd_HLADR", "X175Lu_CD14", "X209Bi_CD16",  "X176Yb_CD56", "X106Cd_CD57", "X158Gd_CD27",
       "X160Gd_CD28", "X169Tm_CD25", "X172Yb_CD38", "X112Cd_CD44", "X143Nd_CD127.IL7Ra",  "X167Er_CCR7.CD197",
       "X155Gd_CD45RA", "X162Dy_FoxP3", "X170Er_CTLA.4.CD152", "X163Dy_CD33",  "X154Sm_CD272.BTLA", "X110Cd_CD107a",
       "X161Dy_IL.17A", "X166Er_IL.10", "X164Dy_Perforin", "X171Yb_GranzymeB", "X148Nd_CD274.PD.L1",  
       "X173Yb_CD273.PD.L2", "X159Tb_GM.CSF",   "X168Er_CD154.CD40L", "X174Yb_CD279.PD1", "X141Pr_CD223.LAG3", 
       "X147Sm_TIM.3",  "X151Eu_CD137.4.1BB",  "X165Ho_IFNg", "X142Nd_IL.1b", "X156Gd_IL.6", "X166Er_IL.10", 
       "X152Sm_TCRgd", "X151Eu_CD137.4.1BB", "X149Sm_IL.12p70", "X150Nd_MIP.1b")



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
o <- orderPanel2
column_cluster <- FALSE

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
o <- orderPanel2
column_cluster <- FALSE

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
o <- orderPanel2
column_cluster <- FALSE


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
o <- orderPanel1
column_cluster <- FALSE

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
o <- orderPanel1
column_cluster <- FALSE


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
o <- orderPanel1
column_cluster <- FALSE


source(fs::path(scriptPath, "FlowSOM_analyse.R"))

