scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE", "TigtherCleanUp")
params$seed <- 6892 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <-  "bonferroni" # "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))


params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE", "TigtherCleanUp")
params$seed <- 7856 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))



params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD4", "TigtherCleanUp")
params$seed <- 3803 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD4", "TigtherCleanUp")
params$seed <- 4903 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))



params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD8", "TigtherCleanUp")
params$seed <- 4303 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD8", "TigtherCleanUp")
params$seed <- 9003 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))


#panel 2

scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2_ALLE", "TigtherCleanUp")
params$seed <- 4392 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))


params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2_ALLE", "TigtherCleanUp")
params$seed <- 5678 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))



params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD4", "TigtherCleanUp")
params$seed <- 3049 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD4", "TigtherCleanUp")
params$seed <- 4038 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))




params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD8", "TigtherCleanUp")
params$seed <- 2390 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD8", "TigtherCleanUp")
params$seed <- 9143 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))




###

# 
# scriptPath <- fs::path("H:", "git", "cytof")
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "surrface")
# #params$seed <- 6892 #nb må endres vil man vil gjøre et annet uttrekk
# #params$muligeK <- c(40, 50 , 60)
# params$fil <- "Prosent_Gating_Surrface_Panel1.csv"
# params$selcected_name <- "surrface"
# params$panel <- 1
# params$adj_p <- 0.25
# params$adj_p_methods <- "bonferroni" #"fdr"  # see help(p.adjust) for other methods 
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_manuel.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))
# 
# 
# 
# 
# scriptPath <- fs::path("H:", "git", "cytof")
# source(fs::path(scriptPath, "gating stim from excelfile Panel 2.R"))
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Stim")
# params$fil <- "Prosent_Gating_Stim_Panel2.csv"
# params$selcected_name <- "stim"
# params$panel <- 2
# params$adj_p <- 0.25 #HM???
# params$adj_p_methods <- "fdr"  # see help(p.adjust) for other methods"bonferroni" # 
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_manuel.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))
# 




# 
# scriptPath <- fs::path("H:", "git", "cytof")
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
# params$fil <- "Prosent_Gating_klustre_Panel2.csv"
# params$selcected_name <- "surrface"
# params$panel <- 2
# params$intervall <- 50 #dvs 50% prediksjonsintervall. 
# params$seed <- "resultat fra flowsom, kluster 2"
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil), ".docx")))




# scriptPath <- fs::path("H:", "git", "cytof")
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
# params$fil <- "Prosent_Gating_klustre_Panel2_Kopi av sjekk disse JB.csv"
# params$selcected_name <- "surrface"
# params$panel <- 2
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", "kopi_JB", ".docx")))



# scriptPath <- fs::path("H:", "git", "cytof")
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
# params$fil <- "Prosent_Gating_klustre_Panel2_CD8.csv"
# params$selcected_name <- "surrface"
# params$panel <- 2
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil), ".docx")))






# scriptPath <- fs::path("H:", "git", "cytof")
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
# params$fil <- "Prosent_Gating_klustre_Panel2_CD8_1og2.csv"
#                
# params$panel <- 2
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil), ".docx")))
# 





# scriptPath <- fs::path("H:", "git", "cytof")
# params <- list()
# params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
# params$fil <- "Prosent_Gating_klustre_Panel2_fra_Summary.csv"
# params$fil2 <- "Panel2_fra_Summary.csv"
# 
# params$panel <- 2
# 
# rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))








scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel2_Gjenskapt.csv"
params$fil2 <- "Panel2_Gjenskapt.csv"
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$seed <- "resultat fra flowsom, kluster 2"

params$panel <- 2

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel2_Gjenskapt_less_markers.csv"
params$fil2 <- "Panel2_Gjenskapt_less_markers.csv"
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$seed <- "resultat fra flowsom, kluster 2"

params$panel <- 2

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "res FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel1_Gjenskapt.csv"
params$fil2 <- "Panel1_Gjenskapt.csv"
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$seed <- "resultat fra flowsom, kluster 1"

params$panel <- 1

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "res FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel1_Gjenskapt_Johanna_less_markers.csv"
params$fil2 <- "Panel1_Gjenskapt_Johanna_less_Markesr.csv"
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$seed <- "resultat fra flowsom, kluster 1, less markers"

params$panel <- 1

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))












#panel 2 1 og 2




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel2_ALLE_CD8_1og2", "TigtherCleanUp")

params$seed <- 9143 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))



params$seed <- 2390 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel1_ALLE_CD8_1og2", "TigtherCleanUp")

params$seed <- 7003 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))



params$seed <- 7733 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
params$adj_p_methods <- "bonferroni" # "fdr"  # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_Bon.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_Bon.docx")))


## regresjon citrus panel 1


scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Resultater fra citrus P1")
params$fil <- "Citrus no 1 panel 1.xlsx"
params$fil2 <- "Citrus no 1 panel 1.csv"
params$seed <- "citrus no 1"
params$adj_p_methods <- "bonferroni"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Resultater fra citrus P1")
params$fil <- "Citrus no 2 panel 1.xlsx"
params$fil2 <- "Citrus no 2 panel 1.csv"
params$seed <- "citrus no 2"
params$adj_p_methods <- "bonferroni"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), "test.docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Resultater fra citrus P1")
params$fil <- "Citrus no 3 panel 1.xlsx"
params$fil2 <- "Citrus no 3 panel 1.csv"
params$seed <- "citrus no 3"
params$adj_p_methods <- "bonferroni"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))





scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Resultater fra citrus P1")
params$fil <- "Run1 selection.xlsx"
params$fil2 <- "Run1 selection.csv"
params$seed <- "Run1 selection"
params$adj_p_methods <- "bonferroni"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))





scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Resultater fra citrus P1")
params$fil <- "Run2 selection.xlsx"
params$fil2 <- "Run2 selection.csv"
params$seed <- "Run2 selection"
params$adj_p_methods <- "bonferroni"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "Resultater fra citrus P1")
params$fil <- "Run3 selection.xlsx"
params$fil2 <- "Run3 selection.csv"
params$seed <- "Run3 selection"
params$adj_p_methods <- "bonferroni"
params$panel <- 1
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))




## regresjon citrus panel 2


scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Resultater fra citrus")
params$fil <- "Citrus no 1 panel 2.xlsx"
params$fil2 <- "Citrus no 1 panel 2.csv"
params$seed <- "citrus no 1"
params$adj_p_methods <- "bonferroni"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Resultater fra citrus")
params$fil <- "Citrus no 2 panel 2.xlsx"
params$fil2 <- "Citrus no 2 panel 2.csv"
params$seed <- "citrus no 2"
params$adj_p_methods <- "bonferroni"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Resultater fra citrus")
params$fil <- "Citrus no 3 panel 2.xlsx"
params$fil2 <- "Citrus no 3 panel 2.csv"
params$seed <- "citrus no 3"
params$adj_p_methods <- "bonferroni"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))





scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Resultater fra citrus")
params$fil <- "Run1 selection panel 2.xlsx"
params$fil2 <- "Run1 selection panel 2.csv"
params$seed <- "Run1 selection"
params$adj_p_methods <- "bonferroni"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))





scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Resultater fra citrus")
params$fil <- "Run2 selection panel 2.xlsx"
params$fil2 <- "Run2 selection panel 2.csv"
params$seed <- "Run2 selection"
params$adj_p_methods <- "bonferroni"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Resultater fra citrus")
params$fil <- "Run3 selection panel 2.xlsx"
params$fil2 <- "Run3 selection panel 2.csv"
params$seed <- "Run3 selection"
params$adj_p_methods <- "bonferroni"
params$panel <- 2
params$adj_p <- 0.1
params$intervall <- 50 #dvs 50% prediksjonsintervall. 
rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fra_citrus.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", gsub(".csv", "",  params$fil2), ".docx")))



