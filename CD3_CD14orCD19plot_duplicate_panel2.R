rm(list = ls())
scriptPath <- fs::path("H:", "git", "cytof")
#scriptPath <- fs::path("C:", "CyToF data", "fra github")

source(fs::path(scriptPath, "readDataToAnalysePanel2.R"))
utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2", "Duplet check")
outFigPath <- utSti 
params <- get_params_fcs_data(fcs_data[[1]])
n_files <- length(fcs_data)

# #posNegPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Analyse i R OUS", "Resultat_Panel_1", "Data")
# posNegPath <- fs::path("C:", "CyToF data", "immun_aga", "posNeg panel2")
#
# posNeg <- readRDS(fs::path(posNegPath, "posNeg.rds"))
# file_names <- readRDS(fs::path(posNegPath, "posNegFilnavn.rds"))
# 


plot_multiple_signal <- function(signal){
  gridExtra::grid.arrange(signal$plotList[[1]], signal$plotList[[2]], signal$plotList[[3]], signal$plotList[[4]], signal$plotList[[5]],
                          signal$plotList[[6]],  signal$plotList[[7]], signal$plotList[[8]], signal$plotList[[9]], signal$plotList[[10]],
                          signal$plotList[[11]], signal$plotList[[12]], signal$plotList[[13]], signal$plotList[[14]], signal$plotList[[15]],
                          signal$plotList[[16]],  signal$plotList[[17]], signal$plotList[[18]], signal$plotList[[19]], signal$plotList[[20]],
                          signal$plotList[[21]], signal$plotList[[22]], signal$plotList[[23]], signal$plotList[[24]], signal$plotList[[25]],
                          signal$plotList[[26]],  signal$plotList[[27]], signal$plotList[[28]], signal$plotList[[29]], signal$plotList[[30]],
                          signal$plotList[[31]], signal$plotList[[32]], signal$plotList[[33]], signal$plotList[[34]], signal$plotList[[35]],
                          signal$plotList[[36]],  signal$plotList[[37]], signal$plotList[[38]], signal$plotList[[39]], signal$plotList[[40]],
                          signal$plotList[[41]], signal$plotList[[42]], signal$plotList[[43]], signal$plotList[[44]], signal$plotList[[45]],
                          signal$plotList[[46]],  signal$plotList[[47]], signal$plotList[[48]], signal$plotList[[49]], signal$plotList[[50]],
                          signal$plotList[[51]], signal$plotList[[52]], signal$plotList[[53]], signal$plotList[[54]], signal$plotList[[55]],
                          signal$plotList[[56]],  signal$plotList[[57]], signal$plotList[[58]], signal$plotList[[59]], signal$plotList[[60]],
                          signal$plotList[[61]], signal$plotList[[62]], signal$plotList[[63]], signal$plotList[[64]], signal$plotList[[65]],
                          signal$plotList[[66]],  ncol = 11, nrow = 6)
}





#gating av CD127

#***************************************************
#pos/neg CD127 ----  #legg inn upper ... 15 feb
#***************************************************
CD3 <- params$name[grep("CD3", params$desc)][1] 
CD14 <- params$name[grep("CD14", params$desc)][1] 
CD19 <- params$name[grep("CD19", params$desc)][1] 

data <-  arc_sinh_transform_selected_channels(fcs_data = fcs_data, channels = c(CD3, CD14, CD19))


signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD3, channel2 = CD14, xname = "CD3",yname = "CD14", plot_title = file_names, title_size = 10)

tiff(fs::path(utSti, "CD3_CD14.tiff"), width = 1800, height = 1200)
plot_multiple_signal(signal)
dev.off()



signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD3, channel2 = CD19, xname = "CD3",yname = "CD19", plot_title = file_names, title_size = 10)

tiff(fs::path(utSti, "CD3_CD19.tiff"), width = 1800, height = 1200)
plot_multiple_signal(signal)
dev.off()




signal <- signal_signal_plot(data = data, random_events = random_events(number_of_events(data)), channel1 = CD14, channel2 = CD19, xname = "CD14",yname = "CD19", plot_title = file_names, title_size = 10)

tiff(fs::path(utSti, "CD14_CD19_panel2.tiff"), width = 1800, height = 1200)
plot_multiple_signal(signal)
dev.off()
