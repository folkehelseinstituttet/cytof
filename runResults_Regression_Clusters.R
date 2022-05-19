scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE", "TigtherCleanUp")
params$seed <- 6892 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))


params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE", "TigtherCleanUp")
params$seed <- 7856 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))



params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD4", "TigtherCleanUp")
params$seed <- 3803 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD4", "TigtherCleanUp")
params$seed <- 4903 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))



params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD8", "TigtherCleanUp")
params$seed <- 4303 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_1_ALLE_CD8", "TigtherCleanUp")
params$seed <- 9003 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))


#panel 2

scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2_ALLE", "TigtherCleanUp")
params$seed <- 4392 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))


params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel_2_ALLE", "TigtherCleanUp")
params$seed <- 5678 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(40, 50 , 60)
params$selcected_name <- "all"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))



params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD4", "TigtherCleanUp")
params$seed <- 3049 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD4", "TigtherCleanUp")
params$seed <- 4038 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD4"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))




params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD8", "TigtherCleanUp")
params$seed <- 2390 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))





params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resulta_Panel_2_ALLE_CD8", "TigtherCleanUp")
params$seed <- 9143 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel1_mars2022", "posNeg", "surrface")
#params$seed <- 6892 #nb må endres vil man vil gjøre et annet uttrekk
#params$muligeK <- c(40, 50 , 60)
params$fil <- "Prosent_Gating_Surrface_Panel1.csv"
params$selcected_name <- "surrface"
params$panel <- 1
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_manuel.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "Stim")
params$fil <- "Prosent_Gating_Stim_Panel2.csv"
params$selcected_name <- "stim"
params$panel <- 2
params$adj_p <- 0.25
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_manuel.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel2.csv"
params$selcected_name <- "surrface"
params$panel <- 2

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel2.csv"
params$selcected_name <- "surrface"
params$panel <- 2

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel2_Kopi av sjekk disse JB.csv"
params$selcected_name <- "surrface"
params$panel <- 2

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))



scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "CleanUpGatingMarch2022", "gating_results_Panel2_mars2022", "posNeg", "result FlowSOM")
params$fil <- "Prosent_Gating_klustre_Panel2_CD8.csv"
params$selcected_name <- "surrface"
params$panel <- 2

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters_fraFlowSOM.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$fil, ".docx")))










#panel 2 1 og 2




scriptPath <- fs::path("H:", "git", "cytof")
params <- list()
params$utSti <- fs::path("F:", "Forskningsprosjekter", "PDB 2794 - Immune responses aga_", "Forskningsfiler", "JOBO", "CyTOF", "Analyse i R OUS", "Resultat_Panel2_ALLE_CD8_1og2", "TigtherCleanUp")

params$seed <- 9143 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))



params$seed <- 2390 #nb må endres vil man vil gjøre et annet uttrekk
params$muligeK <- c(15, 20, 25, 30, 35, 40)
params$selcected_name <- "CD3CD45CD8"
params$panel <- 2
params$adj_p <- 0.1
params$adj_p_methods <- "fdr" # see help(p.adjust) for other methods

rmarkdown::render(fs:::path(scriptPath, "resultater_Regression_clusters.Rmd"), output_file = fs::path(params$utSti, paste0("Result_", params$seed, "_V2.docx")))

