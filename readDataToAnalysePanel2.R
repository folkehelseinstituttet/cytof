# les inn stier
#data_path <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF","Gating", "Gating fra R_FINAL")
scriptPath <- fs::path("H:", "git", "cytof")
data_path <- fs::path("C:", "CyToF data", "immun_aga", "panel2")
#scriptPath <- fs::path("C:", "CyToF data", "fra github")

utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2")
#utSti <- fs::path("C:", "CyToF data", "immun_aga", "posNeg panel2")
metaDataPath <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS")

# les inn funksjoner 
source(fs::path(scriptPath,  "read_data_functions.R"))
source(fs::path(scriptPath,  "transformation_functions.R"))
source(fs::path(scriptPath,  "ploting_functions.R"))
source(fs::path(scriptPath,  "clustering_functions.R"))
source(fs::path(scriptPath,  "gating_functions.R"))
source(fs::path(scriptPath,  "analysis_functions.R"))



#lager dInfo fra filnavnene
posNegFilnavn <- as.character(readRDS(fs::path(utSti,  "Data", "posNegFilnavn.rds")))
posNeg <- readRDS(fs::path(utSti,  "Data", "posNeg.rds"))
# posNegFilnavn <- as.character(readRDS(fs::path(utSti,  "posNegFilnavn.rds")))
# posNeg <- readRDS(fs::path(utSti,   "posNeg.rds"))

status <- rep("Moderat", length(posNegFilnavn))
status[grepl("S_", posNegFilnavn)] <- "Severe"
status[grepl("Ref", posNegFilnavn)] <- "Ref"

age <- rep(NA, length(posNegFilnavn))
for(i in 1:length(posNegFilnavn)){
  if(!substr(posNegFilnavn[i],1,2) %in% "Re")
  age[i] <- strsplit(as.character(posNegFilnavn[i]), "_")[[1]][4]
  if(age[i] %in% c("01","02","1","2","A", "B")){
    age[i] <- strsplit(as.character(posNegFilnavn[i]), "_")[[1]][5]
  }
}
age <- as.numeric(age)


sex <- rep("Male", length(posNegFilnavn))
sex[grepl("Fe", posNegFilnavn)] <- "Female"
sex[grepl("Ref1", posNegFilnavn)] <- "Ref"


dInfo <- data.frame(status = status, age = age, sex = sex, filenames= posNegFilnavn)
rownames(dInfo) <- dInfo$filenames
#dInfo$x <- jitter(as.numeric(dInfo$status))


# tid <- read.csv(fs::path(metaDataPath, "tid.csv"), sep = ";")
# colnames(tid) <- c("fil", "tid")
# tid$pasient <- gsub("T1", "", tid$fil)
# tid$pasient <- gsub("T2", "", tid$pasient)
# 
# 
# dInfo$tid <- NA
# for(i in 1:nrow(tid)){
#   dInfo$tid[grep(as.character(tid$pasient[i]), dInfo$filenames)] <- as.character(tid$tid[i])
# }
# 
# dInfo$tid[grep(as.character("FHI81"), dInfo$filenames)] <- "28.06.2021"
# dInfo$tid[grep(as.character("FHI95"), dInfo$filenames)] <- "28.06.2021"


fcs_files <- fs::path(data_path, rownames(file.info(list.files(data_path))))


# read all files in data_path into one dataset fcs_data
fcs_data_with_info <- read_specific_data_from_folder(data_path = data_path, files_to_open = paste0(as.character(posNegFilnavn), ".fcs"))
fcs_data <- fcs_data_with_info$fcs_data
file_names <- fcs_data_with_info$file_names
rm(fcs_data_with_info)

params_fcs <- get_params_fcs_data(fcs_data[[1]])

kanaler <- params_fcs$desc[grepl("_",params_fcs$desc)]
kanaler <- kanaler[!grepl("_DNA", kanaler)]
kanaler <- kanaler[!grepl("_Cisp", kanaler)]
kanalnavn <- params_fcs$name[params_fcs$desc %in% kanaler]
print("antall markører med")
length(kanalnavn)

